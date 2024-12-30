# Optimization methods and other parameters of the scipy.optimize.minimize function are documented in the reference:
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html
#
# TODO Constraints are definable only for COBYLA, COBYQA, SLSQP and trust-constr methods.
# TODO Automate setting up the mesh prior to optimization with parametrization
# TODO Ensure saving the optimal parameters when exiting
# TODO Implement automated graceful exit from program

import os

import numpy as np
from CSTOptimizer import CSTOptimizer
from scipy.optimize import minimize

CST_PROJECT_PATH = r"C:\Users\marti\Repositories\diploma-thesis\cst\dual-feed\dual_feed.cst"
OPTIMIZATION_CONFIG = {
    "method": "nelder-mead",
    "frequencyRange": [4.7, 5.7],
    "goals": [
        {"result": "1D Results\S-Parameters\S1,1", "target": -15, "weight": 1, "norm": "SoD", "active": True},
        {"result": "1D Results\S-Parameters\S2,2", "target": -15, "weight": 1, "norm": "SoD", "active": True},
        {"result": "1D Results\S-Parameters\S2,1", "target": -40, "weight": 1, "norm": "SoD", "active": True},
    ],
    "variables": [
        {"name": "probe1Distance", "initialValue": 14.53, "bounds": [10, 18], "optimize": True},
        {"name": "probe1Length", "initialValue": 11.67, "bounds": [5, 15], "optimize": True},
        {"name": "gratingDistance", "initialValue": 32.88, "bounds": [20, 50], "optimize": True},
        {"name": "probe2Distance", "initialValue": 14.53, "bounds": [10, 18], "optimize": True},
        {"name": "probe2Length", "initialValue": 11.67, "bounds": [5, 15], "optimize": True},
        {"name": "gratingGap", "initialValue": 3, "bounds": [1, 10], "optimize": True},
        {"name": "wireDiameter", "initialValue": 1, "bounds": [0.5, 3], "optimize": True},
    ],
}
VBA_PATH = os.path.join(os.path.dirname(__file__), "gratingAdjustment.vba")


def main() -> None:
    # Initialize the CSTOptimizer
    optimizer = CSTOptimizer(
        project_path=CST_PROJECT_PATH,
        optimization_config=OPTIMIZATION_CONFIG,
        vba_path=VBA_PATH,
    )

    # Run optimization using the minimize function from scipy.optimize
    result = minimize(
        optimizer.objective_function,
        x0=np.array(optimizer.initial_values),
        method=optimizer.method,
        bounds=optimizer.bounds,
        constraints=optimizer.get_constraints(),
        options=optimizer.get_options(),
    )
    print(f"Optimization finished, result: {result}")


if __name__ == "__main__":
    main()
