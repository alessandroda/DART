# pipelines/base_pipeline.py

"""
Base Assimilation Pipeline
==========================

This module defines the abstract base class for all model–data assimilation
pipelines used within the MIMESI framework.

The class implements a **Template Method pattern**, providing a common
time-stepping orchestration while delegating model-specific and
assimilation-specific logic to concrete subclasses (e.g. FARM, CHIMERE).

Key design principles
---------------------
- The base class is **model-agnostic**
- The base class is **observation-agnostic**
- Time advancement is centralized and consistent
- Concrete pipelines decide *if, when, and how* data assimilation is applied
"""

from abc import ABC, abstractmethod
from datetime import timedelta
import logging
from mimesi_orch.orchestrator_utils import set_date_gregorian

logger = logging.getLogger(__name__)


class BaseAssimilationPipeline(ABC):
    """
    Abstract base class for model–data assimilation pipelines.

    This class controls the **temporal orchestration** of a simulation–
    assimilation workflow, while leaving model execution, observation
    handling, and assimilation logic to subclasses.

    The main execution flow is implemented in :meth:`run_pipeline`.
    """

    def __init__(self, time_manager):
        """
        Initialize the pipeline.

        Parameters
        ----------
        time_manager : TimeManager
            Object responsible for handling start time, end time,
            time stepping, and simulated timestamps.
        """
        self.time_manager = time_manager

    def run_pipeline(self):
        """
        Run the full time loop of the pipeline.

        For each time step, the following sequence is executed:

        1. `before_step`        – optional pre-processing hook
        2. `run_model`          – execute the forward model
        3. Define model time    – set simulated_time for assimilation
        4. `after_model`        – prepare model outputs for DA
        5. `run_assimilation_if_needed`
                               – optional data assimilation step
        6. `after_assimilation` – post-processing of assimilation results
        7. `finalize_step`      – cleanup, logging, bookkeeping

        The loop continues until `time_manager.current_time`
        exceeds `time_manager.end_time`.
        """

        logger.info("[PIPELINE] ---- TIME LOOP START ----")

        while self.time_manager.current_time <= self.time_manager.end_time:
            # Optional hook: e.g. emission perturbations, cleanup
            self.before_step()

            # Run the forward model (mandatory)
            self.run_model()

            # Define the time associated with model outputs
            # (typically current_time + forecast step)
            self.time_manager.simulated_time = (
                self.time_manager.current_time + timedelta(hours=1)
            )

            # Convert simulated_time to (days, seconds) for DA systems
            self.set_days_seconds_model()

            # Optional hook: prepare model outputs for assimilation
            self.after_model()

            logger.info(f"Increment time")
            self.time_manager.increment_time()

            # Perform data assimilation if applicable
            self.run_assimilation_if_needed()

            # Optional hook: map analysis back to the model
            self.after_assimilation()

            # Optional hook: cleanup, logging, archiving
            self.finalize_step()

        logger.info("[PIPELINE] ---- TIME LOOP END ----")

    def set_days_seconds_model(self):
        """
        Convert the simulated model time to Gregorian day/second format.

        This is typically required by data assimilation systems such as DART,
        which represent time using a (days, seconds) convention.
        """
        self.seconds_model, self.days_model = set_date_gregorian(
            self.time_manager.simulated_time.year,
            self.time_manager.simulated_time.month,
            self.time_manager.simulated_time.day,
            self.time_manager.simulated_time.hour,
            self.time_manager.simulated_time.minute,
            self.time_manager.simulated_time.second,
        )

    # ------------------------------------------------------------------
    # Optional hooks (no-op by default)
    # ------------------------------------------------------------------

    def before_step(self):
        """
        Optional hook executed at the beginning of each time step.

        Typical use cases:
        - emission perturbation updates
        - cleanup of old files
        - preparation of working directories
        """
        pass

    def after_model(self):
        """
        Optional hook executed after the forward model run.

        Typical use cases:
        - conversion of model outputs to DA-ready format
        - ensemble post-processing
        """
        pass

    def after_assimilation(self):
        """
        Optional hook executed after the assimilation step.

        Typical use cases:
        - mapping analysis fields back to the model format
        - updating boundary or initial conditions
        """
        pass

    def finalize_step(self):
        """
        Optional hook executed at the end of the time step.

        Typical use cases:
        - cleanup of temporary files
        - logging and diagnostics
        - preparation for next time step
        """
        pass

    # ------------------------------------------------------------------
    # Mandatory methods (must be implemented by subclasses)
    # ------------------------------------------------------------------

    @abstractmethod
    def run_model(self):
        """
        Execute the forward model for the current time step.

        This method **must** be implemented by subclasses and should
        encapsulate all logic required to run the numerical model
        (e.g. FARM, CHIMERE).
        """
        ...

    @abstractmethod
    def run_assimilation_if_needed(self):
        """
        Execute data assimilation if applicable.

        Subclasses decide:
        - whether observations are available
        - whether assimilation should be performed
        - which DA system to use (e.g. DART)

        This method may be a no-op for pipelines that do not
        assimilate data.
        """
        ...
