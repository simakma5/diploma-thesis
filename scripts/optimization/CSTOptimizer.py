# From the official documentation (https://space.mit.edu/RADIO/CST_online/Python/main.html):
# The CST Studio Suite installation comes with Python 3.6, which requires no further setup to start using it with the
# CST Python Libraries. However, this is just a plain Python interpreter with only the standard libraries available.
# You may want to use a more feature rich Python 3.6 distribution like Anaconda with the provided CST Python Libraries.
# To start working with the CST Python Libraries from an external Python distribution, you will have to include the
# directory in Python’s system path.
#
# Possible solutions:
#   - add CST_LIBRARY_PATH to the PYTHONPATH system environment variable
#   - start the script by running `sys.path.append(CST_LIBRARY_PATH)`
#   - solutions specific to IDEs, e.g., VS Code supports editing the PYTHONPATH in the /.env file
#     "PYTHONPATH=C:\Program Files (x86)\CST Studio Suite 2023\AMD64\python_cst_libraries"
#
# While the last option edits the PYTHONPATH variable only when launching from the given IDE, it registers the packages
# to the IDE's language server, allowing it to source the packages for code completion, documentation etc.

import logging
import os
import sys
import time
from typing import Any, Optional

import numpy as np

# Import CST libraries
sys.path.append(r"C:\Program Files (x86)\CST Studio Suite 2024\AMD64\python_cst_libraries")
from cst.interface import DesignEnvironment
from cst.results import ProjectFile, ResultItem

LOG_LEVEL = logging.DEBUG


class CSTOptimizer:
    EVALUATION_TIMEOUT = 10

    def __init__(
        self,
        project_path: str,
        optimization_config: dict[str, Any],
        vba_path: Optional[str] = None,
    ):
        """
        Initializes the CSTOptimizer instance with the necessary CST application environment, project path,
        variable configurations, and optimization goals. The optimizer reads the specified VBA script for
        geometry adjustments and registers a cleanup function to close the CST application upon program exit.

        Args:
            project_path (str): Path to the CST project file.
            optimization_config (dict): Dictionary containing optimization parameters.
            vba_path (Optional[str]): Path to the VBA script used for model geometry adjustment.
        """
        self.logger = logging.getLogger("CSTOptimizer")
        self.configure_logging(log_level=LOG_LEVEL, log_file="optimizer.log")

        # Unpack arguments
        self.project_path = project_path
        self.method = optimization_config["method"]
        self.frequency_min, self.frequency_max = optimization_config["frequencyRange"]
        self.goals = [goal for goal in optimization_config["goals"] if goal["active"]]

        # Runtime variables
        self.step_counter = 0
        self.variable_names, self.initial_values, self.bounds = zip(
            *(
                (var["name"], var["initialValue"], var["bounds"])
                for var in optimization_config["variables"]
                if var["optimize"]
            )
        )
        self.frequency_mask = None
        self.best_parameters = {"objectiveValue": float("inf")}

        # Load VBA script
        if vba_path:
            try:
                with open(vba_path, "r") as vba_file:
                    self.check_geometry = vba_file.read()
            except Exception as e:
                self.logger.error(f"Failed to load VBA script: {e}")
                self.check_geometry = None
        else:
            self.check_geometry = None
            self.logger.warning("No VBA script provided. Geometry adjustments will be skipped.")

        # Initialize CST application
        self.cst = DesignEnvironment().connect_to_any_or_new()
        self.project = self.cst.open_project(self.project_path)
        self.result_module = ProjectFile(self.project_path, allow_interactive=True).get_3d()
        self.last_run_id = self.result_module.get_all_run_ids()[-1]

    def configure_logging(self, log_level=logging.DEBUG, log_file="optimizer.log"):
        """
        Configures the logger for the class with the specified log level and log file.

        Args:
            log_level (int): The logging level (e.g., logging.DEBUG, logging.INFO).
            log_file (str): The file to which logs will be written.
        """
        self.logger.setLevel(log_level)
        if not self.logger.handlers:
            formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")

            # Console handler
            stream_handler = logging.StreamHandler()
            stream_handler.setFormatter(formatter)
            self.logger.addHandler(stream_handler)

            # File handler
            log_file = os.path.join(os.path.dirname(__file__), log_file)
            file_handler = logging.FileHandler(log_file, mode="w")
            file_handler.setFormatter(formatter)
            self.logger.addHandler(file_handler)

    def get_constraints(self) -> list[dict[str, str]]:
        """
        Defines constraints for the optimization process, specifically for the 'GratingGap' and 'WireDiameter'
        parameters if they are part of the optimization.

        Returns:
            list[dict[str, str]]: A list of constraints to be applied during optimization.
        """
        if {"gratingGap", "wireDiameter"}.issubset(self.variable_names):
            gap_index = self.variable_names.index("gratingGap")
            diameter_index = self.variable_names.index("wireDiameter")
            return [{"type": "ineq", "fun": lambda x: x[gap_index] > x[diameter_index]}]
        return []

    def get_options(self) -> dict[str, Any]:
        if self.method == "nelder-mead":
            return {
                "disp": True,  # Display optimization progress
                "initial_simplex": None,  # Initial simplex for Nelder-Mead
                "maxiter": None,  # Max iterations (default is 200 * len(variables))
                "xatol": 1e-2,  # Tolerance for convergence in terms of x (parameter change)
            }
        else:
            return {}

    def update_parameters(self, names: list[str], values: np.ndarray) -> tuple[str, dict[str, Any]]:
        """
        Generates VBA code to update the parameters of the CST model.

        Args:
            x (np.ndarray): Array of parameter values.

        Returns:
            str: VBA script for updating the parameters.
        """
        self.logger.debug("Generating VBA script for parameter updates")
        update_parameters = "Sub Main()\n"
        current_step = {}
        for name, value in zip(names, values):
            self.logger.info(f"{name} = {value}")
            update_parameters += f'StoreParameter("{name}", {value})\n'
            current_step[name] = value
        update_parameters += "RebuildOnParametricChange(True, False)\nEnd Sub\n"
        return update_parameters, current_step

    def objective_function(self, x: np.ndarray) -> float:
        """
        Evaluates cost of the current optimization step.

        Args:
            x (np.ndarray): Array of parameter values.

        Returns:
            float: The computed objective value.
        """
        self.step_counter += 1
        self.logger.info(f"Optimization step {self.step_counter}")

        # Update parameters
        self.logger.debug("Updating parameters")
        update_parameters, current_step = self.update_parameters(self.variable_names, x.round(5))
        self.project.schematic.execute_vba_code(update_parameters, timeout=None)

        # Execute geometry adjustment VBA script, if provided
        if self.check_geometry:
            self.logger.debug("Executing geometry adjustment")
            self.project.schematic.execute_vba_code(self.check_geometry, timeout=None)

        # Run simulation
        self.logger.debug("Running CST simulation")
        try:
            self.project.modeler.run_solver()
            self.last_run_id += 1
        except RuntimeError:
            # TODO Analyse: tolerate failed simulations?
            raise
        finally:
            self.project.save(include_results=True)

        return self.evaluate_step(current_step)

    def evaluate_step(self, current_step: dict[str, Any]) -> float:
        """
        Evaluates the results of the current optimization step.

        Args:
            current_step (dict): The current parameter values.

        Returns:
            float: The computed objective value for this step.
        """
        # Get results
        self.logger.debug(f"Processing results for run ID {self.last_run_id}")
        results = None
        start = time.time()
        while time.time() - start < self.EVALUATION_TIMEOUT:
            try:
                results = [
                    self.parse_results(
                        self.result_module.get_result_item(treepath=goal["result"], run_id=self.last_run_id)
                    )
                    for goal in self.goals
                ]
                break
            except Exception as e:
                self.logger.warning(e)
        if not results:
            self.logger.error(
                f"Run ID {self.last_run_id} not found! Step count: {self.step_counter}, run ID: {self.last_run_id}"
            )
            raise RuntimeError(f"Run ID {self.last_run_id} not found!")

        # Compute objective value
        objective_value = 0
        for result, goal in zip(results, self.goals):
            # Maximum difference
            if goal["norm"] == "MD":
                max_difference = np.max(result - goal["target"])
                self.logger.info(
                    f"Result: {goal['result']} | Weight: {goal['weight']} | "
                    f"Target: {goal['target']} | Maximum difference: {max_difference :.2f}"
                )
                objective_value += goal["weight"] * max_difference
            # Sum of differences
            elif goal["norm"] == "SoD":
                average_difference = sum(result - goal["target"]) / len(result)
                self.logger.info(
                    f"Result: {goal['result']} | Weight: {goal['weight']} | "
                    f"Target: {goal['target']} | Average difference: {average_difference :.2f}"
                )
                objective_value += goal["weight"] * average_difference
            # Neither
            else:
                raise ValueError("Undefined goal norm!")

        # Update best parameters if necessary
        self.logger.info(f"Objective value: {objective_value}\n")
        if objective_value < self.best_parameters["objectiveValue"]:
            current_step["objectiveValue"] = objective_value
            self.best_parameters = current_step

        return objective_value

    def parse_results(self, result_item: ResultItem) -> np.ndarray:
        """
        Parses the S-parameter results from the CST simulation.

        Args:
            result_item (ResultItem): The CST result item containing frequency and complex S-parameter data.

        Returns:
            np.ndarray: The magnitudes of S-parameters in decibels, filtered by the defined frequency range.
        """
        if not self.frequency_mask:
            frequency = np.array(result_item.get_xdata())
            self.frequency_mask = (self.frequency_min <= frequency) & (frequency <= self.frequency_max)
        magnitude = np.abs(result_item.get_ydata())[self.frequency_mask]
        return 10 * np.log10(np.maximum(magnitude, 1e-12))  # Avoid log of zero

    def exit(self):
        self.logger.debug("Program exited - entering the exit method.")
        try:
            if self.best_parameters["objectiveValue"] == float("inf"):
                self.logger.info("No optimization steps completed - reverting to initial values")
                update_parameters, _ = self.update_parameters(self.variable_names, self.initial_values)
            else:
                self.logger.info("Setting the best parameter values found during optimization")
                update_parameters, _ = self.update_parameters(
                    *zip(*((key, value) for key, value in self.best_parameters.items() if key != "objectiveValue"))
                )
            self.project.schematic.execute_vba_code(update_parameters, timeout=None)
            self.logger.debug("Parameters saved")
        except Exception as e:
            self.logger.error(f"Error saving parameters: {e}")
        finally:
            self.cst.close()
