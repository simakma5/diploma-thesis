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

import atexit
import logging
import os
import signal
import sys
import threading
from typing import Any, Optional

import numpy as np

# Import CST libraries
sys.path.append(r"C:\Program Files (x86)\CST Studio Suite 2023\AMD64\python_cst_libraries")
import cst.interface
import cst.results

LOG_LEVEL = logging.DEBUG


class CSTOptimizer:
    interrupt_flag = threading.Event()  # Shared flag to signal interruption

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

        self.project_path = project_path
        self.method = optimization_config["method"]
        self.frequency_min, self.frequency_max = optimization_config["frequencyRange"]
        self.goals = [goal for goal in optimization_config["goals"] if goal["active"]]
        optimized_variables = (
            (var["name"], var["initialValue"]) for var in optimization_config["variables"] if var["optimize"]
        )
        self.variable_names, self.initial_values = zip(*optimized_variables)
        self.step_counter = 0
        self.best_parameters = {"objectiveValue": float("inf")}

        # Load VBA script
        if vba_path:
            try:
                with open(vba_path, "r") as vba_file:
                    self.vba_script = vba_file.read()
            except Exception as e:
                self.logger.error(f"Failed to load VBA script: {e}")
                self.vba_script = None
        else:
            self.vba_script = None
            self.logger.warning("No VBA script provided. Geometry adjustments will be skipped.")

        # Initialize CST application
        self.cst_app = cst.interface.DesignEnvironment()
        self.project = self.cst_app.open_project(self.project_path)
        self.result_module = cst.results.ProjectFile(self.project_path, allow_interactive=True).get_3d()

        # Register SIGINT handler and graceful shutdown
        signal.signal(signal.SIGINT, self.handle_interrupt)
        atexit.register(self.atexit)

    def configure_logging(self, log_level=logging.DEBUG, log_file="optimizer.log"):
        """
        Configures the logger for the class with the specified log level and log file.

        Args:
            log_level (int): The logging level (e.g., logging.DEBUG, logging.INFO).
            log_file (str): The file to which logs will be written.
        """
        self.logger.setLevel(log_level)
        formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")

        # Console handler
        stream_handler = logging.StreamHandler()
        stream_handler.setFormatter(formatter)
        self.logger.addHandler(stream_handler)

        # File handler
        log_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), log_file)
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
        if {"GratingGap", "WireDiameter"}.issubset(self.variable_names):
            gap_index = self.variable_names.index("GratingGap")
            diameter_index = self.variable_names.index("WireDiameter")
            return [{"type": "ineq", "fun": lambda x: x[gap_index] > x[diameter_index]}]
        return []

    def handle_interrupt(self, signum, frame):
        """
        Signal handler for SIGINT to set the interrupt flag but allow solver to finish.
        """
        self.logger.warning("SIGINT received. The program will exit after the solver finishes.")
        self.interrupt_flag.set()

    def run_solver_in_thread(self) -> bool:
        """
        Runs the CST solver in a separate thread, monitoring for interruption.
        Ensures proper thread cleanup and respects the interrupt flag.
        """
        result = {"complete": None}
        exception = {"error": None}

        def solver_task():
            try:
                result["complete"] = self.project.modeler.run_solver()
            except Exception as e:
                exception["error"] = e

        thread = threading.Thread(target=solver_task)
        thread.start()
        while thread.is_alive():
            if self.interrupt_flag.is_set():
                self.logger.warning("Interrupt flag set. Waiting for CST to finish the current simulation...")
                break
            thread.join(timeout=0.1)
        thread.join()
        if exception["error"]:
            raise exception["error"]
        return result["complete"]

    def objective_function(self, x: np.ndarray) -> float:
        """
        Synchronous objective function that leverages the threaded solver.

        Args:
            x (np.ndarray): Array of parameter values to be evaluated, given by the optimization algorithm.

        Returns:
            float: The computed objective function value.
        """
        self.step_counter += 1
        current_step = dict()
        simulation_complete = False
        self.logger.info(f"Optimization step {self.step_counter}")

        # if self.step_counter == 1:
        #     results = cst.results.ProjectFile(self.project_path).get_3d()
        #     run_ids = results.get_all_run_ids()
        #     last_run_id = results.get_all_run_ids()[-1]
        #     self.logger.debug(f"Processing run ID {last_run_id}")

        # Open the CST project and apply parameter changes
        try:
            # Update parameters
            self.logger.debug("Updating parameters")
            parameter_update = """Sub Main ()\n"""
            for name, value in zip(self.variable_names, x):
                self.logger.info(f"{name} = {value}")
                current_step[f"{name}"] = value
                parameter_update += f"""StoreParameter("{name}", {value})\n"""
            parameter_update += """RebuildOnParametricChange(True, False)\nEnd Sub\n"""
            self.project.schematic.execute_vba_code(parameter_update, timeout=None)

            # Execute geometry adjustment VBA script
            if self.vba_script:
                self.logger.debug("Executing geometry adjustment VBA script")
                self.project.schematic.execute_vba_code(self.vba_script.read(), timeout=None)

            # Run the simulation in a background thread
            self.logger.debug("Running simulation in a background thread")
            simulation_complete = self.run_solver_in_thread()

        except Exception as e:
            self.logger.error(f"Error during project handling: {e}")
        finally:
            self.project.save(include_results=True)

        if not simulation_complete:
            self.logger.error("Simulation did not complete successfully.")
            return float("inf")

        return self.evaluate_step(current_step)

    def evaluate_step(self, current_step):
        # Retrieve and parse results
        last_run_id = self.result_module.get_all_run_ids()[-1]
        self.logger.debug(f"Processing run ID {last_run_id}")
        results = [
            self.parse_results(self.result_module.get_result_item(treepath=goal["result"], run_id=last_run_id))
            for goal in self.goals
        ]

        # Calculate the objective value
        objective_value = 0
        for result, goal in zip(results, self.goals):
            maximum_difference = np.max(np.abs(result - goal["target"]))
            self.logger.info(
                f"Processed result: {goal['result']} | Weight: {goal['weight']} | "
                f"Target: {goal['target']} | Maximum difference: {maximum_difference:.2f}"
            )
            objective_value += goal["weight"] * maximum_difference

        # Compare with previous results and return
        self.logger.info(f"Objective value (weighted sum): {objective_value}\n")
        if objective_value < self.best_parameters["objectiveValue"]:
            current_step["objectiveValue"] = objective_value
            self.best_parameters = current_step
        return objective_value

    def parse_results(self, result_item: cst.results.ResultItem) -> np.ndarray:
        """
        Parses the S-parameter results from the CST simulation.

        Args:
            result_item (cst.results.ResultItem): The CST result item containing frequency and complex S-parameter data.

        Returns:
            np.ndarray: The magnitudes of S-parameters in decibels, filtered by the defined frequency range.
        """
        frequency = np.array(result_item.get_xdata())
        magnitude = np.abs(result_item.get_ydata())[(frequency > self.frequency_min) & (frequency < self.frequency_max)]
        return 10 * np.log10(np.maximum(magnitude, 1e-12))  # Avoid log of zero

    def atexit(self):
        """
        Saves the best parameters and closes CST gracefully at program exit.
        """
        self.logger.debug("Program exited - entering the atexit method.")
        try:
            parameter_update = """Sub Main ()\n"""
            if self.best_parameters["objectiveValue"] == float("inf"):
                self.logger.info("No optimization steps completed - reverting to initial values")
                for name, value in zip(self.variable_names, self.initial_values):
                    self.logger.info(f"{name} = {value}")
                    parameter_update += f"""StoreParameter("{name}", {value})\n"""
            else:
                self.logger.info("Setting the best parameter values found during optimization")
                for name, value in self.best_parameters.items():
                    if name != "objectiveValue":
                        self.logger.info(f"{name} = {value}")
                        parameter_update += f"""StoreParameter("{name}", {value})\n"""
            parameter_update += """RebuildOnParametricChange(True, False)\nEnd Sub\n"""
            self.project.schematic.execute_vba_code(parameter_update, timeout=None)
            self.logger.debug("Parameters saved")
        except Exception as e:
            self.logger.error(f"Error during atexit: {e}")
        finally:
            # Ensure CST closes
            self.logger.debug("Closing CST application.")
            self.cst_app.close()
