# Installation

## Conda environment for the orchestrator (CINECA G100)

For your operational setup, the orchestrator is expected to run using the conda environment located at:

- `/g100_work/ARPAE_AQM/MIMESI/miniconda3/envs/mimesi_env`

In practice, you only need to activate it before running `main.py`:

```bash
conda activate mimesi_env
```

If you are inside a batch script, do this after sourcing the CINECA module environment (see the section below).

## External tools (used when available)

Some steps may call external command-line tools if they are available on your system:

- `ncks` (NCO) to subset NetCDF files (with an xarray fallback)
- `cdo` is invoked in places with warnings filtering (see `_run_cdo_suppress_5d` in `pipelines/chimere2017/pipeline.py`)

If these tools are missing, behavior may degrade or fail depending on which code paths your run exercises.

## CINECA G100: required module environment (DART/CHIMERE)

On CINECA G100, DART/CHIMERE runs typically require a specific module environment to be loaded **before job submission** (and therefore inside the batch script that runs the model).

If your site provides an environment script (for example `env_cineca`), make sure it sources a module block like the following (as provided for your operational setup):

```bash
echo "------------------------"
module purge
module load profile/base
module load /cineca/prod/opt/modulefiles/base/compilers/gcc/10.2.0
module load autoload intelmpi/oneapi-2022--binary \
                     intel-oneapi-compilers/2021.4.0 \
                     zlib/1.2.11--intel--2021.4.0 \
                     libszip/2.1.1--intel--2021.4.0 \
                     hdf5/1.10.7--intel-oneapi-mpi--2021.4.0--intel--2021.4.0 \
                     netcdf-c/4.8.1--intel-oneapi-mpi--2021.4.0--intel--2021.4.0 \
                     netcdf-fortran/4.5.3--intel-oneapi-mpi--2021.4.0--intel--2021.4.0 \
                     nco/5.0.1--intel-oneapi-mpi--2021.4.0--intel--2021.4.0
module load cdo/1.9.10--gcc--10.2.0
echo "------------------------"
```

Operational rule:

- Always source this environment (directly or via `env_cineca`) **inside your Slurm job script** before running DART/CHIMERE or the orchestrator.



## Configuration validation

Configuration is validated strictly. Unknown YAML keys cause startup to fail.

!!! note
    The file `config/config_mimesi_ITA7.yaml` is the recommended example configuration to start from in this repository. Some other example configs (for example `config_orchestrator.yaml`) may be legacy/older and may not validate against the current configuration schema.
