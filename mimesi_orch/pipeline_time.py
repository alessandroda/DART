from __future__ import annotations

from dataclasses import dataclass
from datetime import timedelta
import logging

import pandas as pd


logger = logging.getLogger(__name__)


@dataclass
class AssimWindow:
    start_time: pd.Timestamp
    end_time: pd.Timestamp
    run_hours: int
    has_assimilation: bool
    obs_time: pd.Timestamp | None = None
    orbit_filename: str | None = None


class TimeManager:
    def __init__(self, start_time: str, end_time: str, dt_seconds: int):
        """
        Initializes the TimeManager with the start, end times and the time delta.
        """
        self.start_time = pd.to_datetime(start_time)
        self.end_time = pd.to_datetime(end_time)
        self.current_time = self.start_time
        self.simulated_time = None
        self.timestamp_farm_run = None
        self.sat_obs = None
        self.dt = pd.Timedelta(dt_seconds, unit="s")
        self.last_perturbed_day = None
        self.end_file_date_control_run = self.start_time - timedelta(days=1)
        self.end_file_date = None
        self.slot_time = None
        self.run_hours = None
        self.end_file_datetime = self.current_time - timedelta(hours=1)
        self.prev_run_hours = None

        self.check_start_ahead_end()

    def check_start_ahead_end(self):
        """check if end time follows start time"""
        if self.start_time > self.end_time:
            logger.warning(
                "End time (%s) is ahead start time %s",
                self.end_time,
                self.start_time,
            )

    def increment_time(self):
        """
        Increments the current time by the delta (dt).
        """
        self.current_time += self.dt

    def shifted_time(self, step_offset: int = 0) -> pd.Timestamp:
        """
        Return current_time shifted by N assimilation timesteps.

        step_offset=0 -> current_time
        step_offset=1 -> current_time + dt
        step_offset=-1 -> current_time - dt
        """
        return self.current_time + (step_offset * self.dt)

    def formatted_time(
        self,
        step_offset: int = 0,
        time_format: str = "%Y%m%d%H",
    ) -> str:
        """Return a shifted timestamp already formatted as string."""
        return self.shifted_time(step_offset).strftime(time_format)

    @staticmethod
    def round_to_closest_hour(timestamp):
        if timestamp.minute >= 30:
            rounded_timestamp = timestamp.replace(minute=0, second=0) + timedelta(
                hours=1
            )
        else:
            rounded_timestamp = timestamp.replace(minute=0, second=0)
        return rounded_timestamp

    def is_within_bounds(self):
        """
        Checks if the current time is within the start and end bounds.
        """
        return self.current_time <= self.end_time

    def get_formatted_time(self, time_format: str = "%Y%m%d_%H%M%S"):
        """
        Returns the current time formatted as a string according to the specified format.
        Default format is YYYYMMDD_HHMMSS.
        """
        return self.current_time.strftime(time_format)

    def update_simulated_time(self, new_time):
        """
        Updates the simulated time to track the last simulated timestamp.
        """
        self.simulated_time = new_time
