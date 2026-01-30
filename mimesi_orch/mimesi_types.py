from enum import Enum

class ModelType(str, Enum):
    FARM = "farm"
    CHIMERE = "chimere"

class Scheduler(Enum):
    SLURM = "slurm"
    LSF = "lsf"


