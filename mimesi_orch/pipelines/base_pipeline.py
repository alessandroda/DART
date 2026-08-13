# pipelines/base_pipeline.py

"""
Base Assimilation Pipeline
==========================

This module defines the abstract base class for all model-data assimilation
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
import logging
import pipeline_errors
from pipeline_time import AssimWindow
from time_utils import set_date_gregorian

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
        self.current_window: AssimWindow | None = None

    @staticmethod
    def _fmt_timestamp(value) -> str:
        if value is None:
            return "None"
        return value.strftime("%Y-%m-%d %H:%M:%S")

    def _log_time_context(self, phase: str) -> None:
        logger.info(
            "[TIME] %s current_time=%s simulated_time=%s dt=%s end_time=%s",
            phase,
            self._fmt_timestamp(self.time_manager.current_time),
            self._fmt_timestamp(self.time_manager.simulated_time),
            self.time_manager.dt,
            self._fmt_timestamp(self.time_manager.end_time),
        )

    def build_assim_window(self) -> AssimWindow:
        """
        Build the timing window for the current cycle.

        The default implementation preserves the existing hourly behavior.
        Pipelines can override this to model longer free-forecast windows.
        """
        return AssimWindow(
            start_time=self.time_manager.current_time,
            end_time=self.time_manager.current_time + self.time_manager.dt,
            run_hours=1,
            has_assimilation=False,
        )

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

        logger.info("[STEP] ---- TIME LOOP START ----")

        while self.time_manager.current_time <= self.time_manager.end_time:
            try:
                self._log_time_context("step_start")
                self.current_window = self.build_assim_window()
                logger.info(
                    "[TIME] window start=%s end=%s run_hours=%s has_assimilation=%s",
                    self._fmt_timestamp(self.current_window.start_time),
                    self._fmt_timestamp(self.current_window.end_time),
                    self.current_window.run_hours,
                    self.current_window.has_assimilation,
                )
                self.time_manager.slot_time = self.current_window.end_time
                self.time_manager.run_hours = self.current_window.run_hours

                # --- computing next slot time and hours ---
                if self.current_window.end_time < self.time_manager.end_time:
                    original_current_time = self.time_manager.current_time
                    try:
                        self.time_manager.current_time = self.current_window.end_time
                        next_window = self.build_assim_window()
                        self.time_manager.next_slot_time = next_window.end_time
                        self.time_manager.next_run_hours = next_window.run_hours
                    except Exception as e:
                        logger.debug("Next window computation failed (es. final time reached): %s", e)
                    finally:
                        self.time_manager.current_time = original_current_time
                # ---------------------------------------------------

                # Optional hook: e.g. emission perturbations, cleanup
                self.before_step()

                # Run the forward model (mandatory)
                self.run_model()

                # Define the time associated with model outputs
                # (typically current_time + forecast step)
                self.time_manager.simulated_time = self.current_window.end_time
                #self.time_manager.simulated_time = (
                #    self.time_manager.current_time + hours_simulated
                #)
                self._log_time_context("after_model_set_simulated_time")


                # Convert simulated_time to (days, seconds) for DA systems
                self.set_days_seconds_model()
                logger.info(
                    "[TIME] gregorian_conversion simulated_time=%s days=%s seconds=%s",
                    self._fmt_timestamp(self.time_manager.simulated_time),
                    self.days_model,
                    self.seconds_model,
                )

                # Optional hook: prepare model outputs for assimilation
                self.after_model()

                logger.info(
                    "[TIME] increment current_time %s -> %s",
                    self._fmt_timestamp(self.time_manager.current_time),
                    self._fmt_timestamp(self.current_window.end_time),
                )
                self.time_manager.current_time = self.current_window.end_time
                self._log_time_context("after_increment_before_assimilation")

                # Perform data assimilation if applicable
                try:
                   self.run_assimilation_if_needed()
                except pipeline_errors.SkipAssimilation as e:
                    logger.info(f"[DART] Skipped: {e}")

                # Optional hook: map analysis back to the model
                self.after_assimilation()

                # Optional hook: cleanup, logging, archiving
                self.finalize_step()
                self._log_time_context("step_end")

            except pipeline_errors.FatalPipelineError as e:
                logger.critical(f"[PIPELINE] Fatal error: {e}")
                raise

            except pipeline_errors.PipelineError as e:
                logger.error(f"[PIPELINE] Error: {e}")
                raise

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
