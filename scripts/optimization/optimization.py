# Optimization methods and other parameters of the scipy.optimize.minimize function are documented in the reference:
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html
#
# TODO Nelder-Mead Simplex Algoritm doesn't support defining constraints
# TODO Automate setting up the mesh prior to optimization with parametrization
# TODO Disable all TBPP to avoid failures
# TODO Ensure that the class always parses the correct results
# TODO Test saving the optimal parameters when exiting (EXIT GRACEFULLY)

import os

import numpy as np
from CSTOptimizer import CSTOptimizer
from scipy.optimize import minimize  # noqa

CST_PROJECT_PATH = (
    r"C:\Users\marti\Repositories\Workspace\IDE projects\CST\diploma_thesis"
    r"\coax-to-waveguide-adapter\CoaxToWaveguideAdapter.cst"
)
OPTIMIZATION_CONFIG = {
    "method": "nelder-mead",
    "frequencyRange": [5, 6],
    "goals": [
        {"result": "1D Results\S-Parameters\S1,1", "target": -12, "weight": 1, "active": True},
        {"result": "1D Results\S-Parameters\S2,2", "target": -12, "weight": 1, "active": False},
        {"result": "1D Results\S-Parameters\S2,1", "target": -40, "weight": 0.7, "active": False},
    ],
    "variables": [
        {"name": "Probe1Distance", "initialValue": 15.01, "bounds": [10, 16.25], "optimize": True},
        {"name": "Probe1Length", "initialValue": 12.22, "bounds": [5, 15], "optimize": True},
        {"name": "Probe2Distance", "initialValue": 15.01, "bounds": [10, 16.25], "optimize": False},
        {"name": "Probe2Length", "initialValue": 12.22, "bounds": [5, 15], "optimize": False},
        {"name": "GratingDistance", "initialValue": 48.76, "bounds": [30, 50], "optimize": False},
        {"name": "GratingGap", "initialValue": 3, "bounds": [1, 10], "optimize": False},
        {"name": "WireDiameter", "initialValue": 1, "bounds": [0.5, 3], "optimize": False},
    ],
}
VBA_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "gratingAdjustment.vba")


def main() -> None:
    """
    The main function that initializes the optimizer and runs the CST optimization using the Nelder-Mead method.
    This function configures logging, initializes the CSTOptimizer with necessary parameters, and runs the optimization.
    """
    # Initialize the CSTOptimizer with the project, variable, goal, and frequency range
    optimizer = CSTOptimizer(
        project_path=CST_PROJECT_PATH,
        optimization_config=OPTIMIZATION_CONFIG,
        # vba_path=VBA_PATH,
    )

    # # Test run
    # optimizer.objective_function(np.array([round(var["initialValue"] + 3) for var in PARAMETERS if var["optimize"]]))

    # Run optimization using the minimize function from scipy.optimize
    optimization_variables = [var for var in OPTIMIZATION_CONFIG["variables"] if var["optimize"]]
    result = minimize(
        optimizer.objective_function,  # Objective function to minimize
        x0=np.array([var["initialValue"] for var in optimization_variables]),  # Initial guess
        method=OPTIMIZATION_CONFIG["method"],  # Optimization method (Nelder-Mead)
        bounds=[var["bounds"] for var in optimization_variables],  # Bounds for each variable
        constraints=optimizer.get_constraints(),  # Constraints for the optimization
        options={
            "disp": True,  # Display optimization progress
            "initial_simplex": None,  # Initial simplex for Nelder-Mead
            "maxiter": None,  # Max iterations (default is 200 * len(variables))
            "xatol": 1e-2,  # Tolerance for convergence in terms of x (parameter change)
        },
    )
    print(f"Optimization finished, result: {result}")


if __name__ == "__main__":
    main()
