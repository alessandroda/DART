# Quickstart

This quickstart is written for **scientists and operational users**. It focuses on what to prepare and which commands to run.

## Step 1 — Start from the provided example config

Use the repository example `config/config_mimesi_ITA7.yaml` as your starting point.

Recommended workflow:

```bash
cp config/config_mimesi_ITA7.yaml run_config.yaml
```

Then edit `run_config.yaml` for your environment. In practice, the fields you most often need to change are:

- `paths.base_path`: the “root” used to interpret relative paths in `paths.*`
- `paths.path_data`: where run data and default logs are written
- `paths.listing_file`: observation listing file location (absolute, or relative to `paths.base_path`)
- `time.start_time`, `time.end_time`, `time.dt_seconds`: your run period and step
- `cluster.*`: queue/project/scheduler values for your system

!!! note
    In `config/config_mimesi_ITA7.yaml`, `time.start_time` and `time.end_time` are the same timestamp. For a multi-step run, set `end_time` later than `start_time`.

## Step 2 — Place required input files where the config expects them

At minimum, ensure the file referenced by `paths.listing_file` exists.

In the CHIMERE 2017 pipeline, this file is read as a `;`-separated table and is expected to contain a `start_time` column (used to decide when assimilation is performed).

### Create the ensemble run folders (under `paths.path_data`)

For the ITA7 example, `assimilation.no_mems: 12`. Operationally, you typically prepare one run folder per member under `paths.path_data`:

```text
<paths.path_data>/
  RUN_0/
  RUN_1/
  ...
  RUN_11/
```

If you also need to provide an initial-condition (IC) file to start the assimilation, one operational pattern is to link the same IC into each member folder (adjust the IC filename/path to your case).

Example commands to run **inside** `paths.path_data`:

```bash
for i in {0..11}; do mkdir -p RUN_$i; done
for i in {0..11}; do ln -s /g100_work/ARPAE_AQM/MIMESI/runs_chimere/input_testcase_2025/end.2025110100_2025110200_ITA7.nc RUN_$i/end.2025110100_2025110200_ITA7.nc; done
```

!!! important
    The IC path above is site-specific. Replace it with the correct IC file for your run, and ensure the number of members (`0..11`) matches `assimilation.no_mems - 1`.

### Satellite data (observation inputs)

In the ITA7 example config, the observation listing is:

- `paths.listing_file: SAT_OBS/C03__listing.csv`

Because this is a **relative** path, the file must exist at:

- `<paths.base_path>/SAT_OBS/C03__listing.csv`

Operational expectations inferred from the enabled workflow:

- the file is `;`-separated
- it must contain a `start_time` column (timestamps)
- it must contain a `filename` column (used to locate an orbit/product file)

The workflow treats `filename` as a path relative to `paths.base_path` (it is joined with `base_path`). So if the listing contains `SAT_OBS/some_orbit_file.nc`, that file must exist at:

- `<paths.base_path>/SAT_OBS/some_orbit_file.nc`

The time matching uses `satellite_data.search_window_seconds` (in the ITA7 example: `1800`, i.e. ±30 minutes around a window end).

### Freerun / DART working directories (what gets written during a run)

The ITA7 example config includes a `dart:` section with paths such as `RUN/data_freerun/analysis`.

Operationally, what matters for this repository is where the workflow actually writes its per-cycle products.
In the enabled CHIMERE 2017 workflow here, DART-related products are written under `paths.path_data` using these folders:

- `analysis/<timestamp>/`
- `prior/<timestamp>/`
- `preassim/<timestamp>/`
- `posteriors/<timestamp>/`

So, operationally, make sure:

- `paths.path_data` points to a location where you have write permission
- there is enough disk space for per-cycle outputs

### If you use emission perturbations (recommended for the ITA7 example)

The example config `config/config_mimesi_ITA7.yaml` uses:

- `assimilation.emi_perturbations` (for example `NO2: NO2_2000`)
- `paths.path_perturbed_emi` (where the perturbation files are expected)

This means you must generate (or provide) **perturbed emission NetCDF files** for each ensemble member and each emission variable listed in `assimilation.emi_perturbations`.

#### Expected folder layout

Based on the enabled CHIMERE 2017 workflow, the orchestrator expects perturbations under:

```text
<paths.path_perturbed_emi>/
  emi_0/
    <emi_dir>/
      AEMISSIONS.<YYYYMMDDHH>_<YYYYMMDDHH>_ITA7_0.nc
  emi_1/
    <emi_dir>/
      AEMISSIONS.<YYYYMMDDHH>_<YYYYMMDDHH>_ITA7_1.nc
  ...
```

Where:

- `<emi_dir>` is the value from `assimilation.emi_perturbations` (for example `NO2_2000`)
- the date stamp must match the day being processed (daily file, from 00:00 to 00:00 next day)
- the final `_N.nc` suffix is the ensemble member index

#### Scripts referenced in operations

This repository includes bash helpers for perturbations under:

- `perturbation/bash_emissions/`

The main scripts are:

- `perturbation/bash_emissions/perturb_fields_emission_mimesi.bsh`
  - sets `EMISSION_...` variables and runs `perturbation/perturb_emi_opt.py` for one day (and one variable)
- `perturbation/bash_emissions/submit_perturb_fields_emission_mimesi.bsh`
  - example Slurm submission wrapper to run `perturb_fields_emission_mimesi.bsh` over a date range
- `perturbation/bash_emissions/submit_replace_perturb_into_original_emission_arg.sh`
  - helper used by the workflow to apply perturbations into emission files (this path is referenced by `paths.run_submit_replace_perturbations` in the ITA7 config)

If you use these scripts, you typically:

1. Edit the script headers and paths (e.g. `BASE=...`, conda activation, account/partition) to match your environment.
2. Run `perturb_fields_emission_mimesi.bsh` for each variable you want to perturb (for the ITA7 config, this is at least `NO` and `NO2`).

Example (single day, one variable):

```bash
cd perturbation/bash_emissions
./perturb_fields_emission_mimesi.bsh 2025-12-01
```

Then check that files are created under the directory configured by `paths.path_perturbed_emi` (in the ITA7 config this is `perturbations/emi_mems`, relative to `paths.base_path`).

#### The perturbation generator: `perturbation/perturb_emi_opt.py`

This repository includes a Python program that creates emission perturbations:

- `perturbation/perturb_emi_opt.py`

Operationally, it:

1. Finds one or more “base” emission NetCDF files (glob pattern).
2. For each time step in the file, generates correlated random perturbation factors.
3. Applies those factors to the selected emission variable.
4. Writes one output NetCDF per ensemble member, under a member directory.

It is configured via environment variables with prefix `EMISSION_`.

Common environment variables (defaults are defined inside the script):

- `EMISSION_EMISSION_BASE_DIR`: folder containing the base emission NetCDF files
- `EMISSION_NAME_NETCDFS`: file pattern (for example `AEMISSIONS.*.nc`)
- `EMISSION_VAR`: emission variable name to perturb (for example `NO2` or `NO`)
- `EMISSION_MEMS`: number of ensemble members
- `EMISSION_SUB_DIR_EMI`: output “tag” used in the output directory name (for example `2000`)
- `EMISSION_PATH_EMISSIONS`: output root folder (this should contain the `emi_mems/` directory referenced by `paths.path_perturbed_emi`)

Example (adjust paths to your system):

```bash
export EMISSION_PATH_EMISSIONS="/path/to/perturbations"
export EMISSION_EMISSION_BASE_DIR="/path/to/base_emissions/"
export EMISSION_NAME_NETCDFS="AEMISSIONS.*_ITA7.nc"
export EMISSION_VAR="NO2"
export EMISSION_MEMS="12"
export EMISSION_SUB_DIR_EMI="2000"

python -u perturbation/perturb_emi_opt.py
```

Output layout produced by `perturb_emi_opt.py`:

```text
<EMISSION_PATH_EMISSIONS>/
  emi_mems/
    emi_0/
      <EMISSION_VAR>_<EMISSION_SUB_DIR_EMI>/
        <base_filename>_0.nc
    emi_1/
      <EMISSION_VAR>_<EMISSION_SUB_DIR_EMI>/
        <base_filename>_1.nc
    ...
```

To use these outputs in the orchestrator run, ensure the **final files on disk** match what the orchestrator will search for:

- location: `<paths.path_perturbed_emi>/emi_<mem>/<emi_dir>/...`
- file name: `AEMISSIONS.<YYYYMMDDHH>_<YYYYMMDDHH>_ITA7_<mem>.nc`

The provided `perturb_fields_emission_mimesi.bsh` script sets `EMISSION_NAME_NETCDFS` to an `AEMISSIONS.<stamp>_ITA7.nc` filename, so the outputs produced by `perturb_emi_opt.py` follow the same naming pattern with an added `_<mem>.nc` suffix.

## Step 3 — Run the orchestrator

The repository provides a CLI entry point in `main.py` that expects a YAML config path:

```bash
python -u main.py -c path/to/config.yaml
```

Example using the copied config:

```bash
python -u main.py -c run_config.yaml
```

### Example: run on CINECA G100 (Slurm batch script)

If you run on CINECA G100, you typically submit the orchestrator as a Slurm job and source the CINECA environment before running.

Example batch script:

```bash
#!/bin/bash
#SBATCH --partition=g100_usr_prod
#SBATCH --account=arpae_aqm
#SBATCH --time=23:59:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --job-name=orchestrator

ulimit -s unlimited
source /g100_work/ARAPE_AQM/MIMESI/env_cineca
conda activate mimesi_env
cd /g100_work/ARAPE_AQM/MIMESI/DART/mimesi_orch

python -u main.py -c config/config_mimesi_ITA7.yaml
```

!!! note
    Paths, account name, and environment name are site-specific. Adjust them to your system and user account.

### Activate the recommended environment (ITA7 example)

The ITA7 example configuration points to a specific Python installation via:

- `paths.env_python: miniconda3/bin/python3.13` (relative to `paths.base_path`)

Operationally, you should run the orchestrator using the **same conda installation** (or an equivalent environment with the required packages installed).

If your site uses a `conda init` block like the one below (from your environment setup), source it (typically via your shell startup files) and activate the conda environment used for MIMESI runs before launching the orchestrator:

```bash
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/g100_work/ARPAE_AQM/MIMESI/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/g100_work/ARPAE_AQM/MIMESI/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/g100_work/ARPAE_AQM/MIMESI/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/g100_work/ARPAE_AQM/MIMESI/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<
```

Then (example):

```bash
conda activate <your_mimesi_env>
python -u main.py -c run_config.yaml
```

If you do not use conda, ensure your Python environment can import the required packages (see [Installation](installation.md)).

## Step 4 — Follow progress in the log

When the run starts, it prints and writes status messages, including:

- which configuration file was loaded
- where the run data directory is (`paths.path_data`)
- where logs are being written
- the time being processed at each step

See the “Logs” section at the end of this page for where log files are created.

## Logs

`main.py` writes a timestamped log file under a log directory derived from your configuration.

If you don’t specify a custom log name, the default folder under `paths.path_data` is `mimesi_orchestrator_logs`.
