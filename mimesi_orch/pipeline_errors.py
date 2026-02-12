class PipelineError(Exception):
    """Base class for all pipeline errors"""


class FatalPipelineError(PipelineError):
    """Non-recoverable error: stop everything"""


class ModelRunError(PipelineError):
    """Forward model failed or not submitted correctly"""


class AssimilationError(PipelineError):
    """Assimilation step failed"""


class SkipAssimilation(Exception):
    """Legitimate skip (no obs, disabled, etc.)"""


class SchedulerError(Exception):
    """Error in submitting the job"""
