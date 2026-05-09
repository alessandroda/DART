# Explaining `config/config_mimesi_ITA7.yaml` (Example)

This page explains the repository example configuration file `config/config_mimesi_ITA7.yaml`.

It is intended for **scientists and operational users**: what each section controls, which values you typically change, and where files must be placed.

!!! important
    This file contains **site-specific paths** (for a particular system). You will almost certainly need to adapt paths and cluster settings for your environment.

## How to use this example

1. Copy it to a run-specific config:

   ```bash
   cp config/config_mimesi_ITA7.yaml run_config.yaml
   ```

2. Edit `run_config.yaml` to match your environment (paths, times, queue/project, etc.).
3. Run:

   ```bash
   python -u main.py -c run_config.yaml
   ```

## Section-by-section guide

### `assimilation`

This section defines the “what” of the assimilation run (variable names, ensemble size, and which workflow is selected).

From the example:

- `model_type: CHIMERE`
  - This repository only supports the CHIMERE workflow (CHIMERE 2017 pipeline).
- `ass_var: NO2`
  - Name of the state variable being assimilated (as used by the workflow and your external setup).
- `state_variable_qty: QTY_NO2`
  - Quantity label associated with the state variable.
- `obs_type: SAT_NO2_TROPOMI`
  - Observation type label used by the workflow and your external inputs.
- `case_dir: ITA7`
  - Case name/directory identifier used by the workflow.
- `no_mems: 12`
  - Ensemble size (number of members).
- `run_assimilation_flag: true`
  - Whether assimilation is enabled.

#### `emi_perturbations`

The example includes:

```yaml
emi_perturbations:
  'NO': NO_2000
  NO2: NO2_2000
```

This is a mapping from an emission variable name to a perturbation subdirectory name.

Operational note:

- YAML can accidentally interpret tokens like `NO` as boolean values unless quoted; the example correctly quotes `'NO'`.

### `time`

Controls the time period and step.

- `dt_seconds: 3600`
  - The time step (seconds). `3600` means hourly stepping.
- `start_time`, `end_time`
  - Start and end timestamps.

!!! important
    In the example, `start_time` and `end_time` are the same (`'2025-12-01 00:00:00'`). That results in a single time step. For a multi-step run, set `end_time` later than `start_time`.

The example also configures optional “backup” settings:

- `backup_perturb_days: 2`
- `backup_ic_hours: 48`
- `backup_ic_option: hourly`

These options are passed through configuration. Their operational meaning depends on how your external run environment and workflow use them.

### `paths`

This is the most operationally important section: it tells the orchestrator where inputs, templates, and run directories are.

#### Absolute vs relative paths

In this repository:

- if a `paths.*` value is absolute (starts with `/`), it is used as-is
- if it is relative, it is interpreted relative to `paths.base_path`

#### Key entries in the example

- `base_path: /g100_work/ARPAE_AQM/MIMESI/`
  - Root used to interpret relative `paths.*` entries.
- `path_data: /g100_scratch/userexternal/adausili/ITA7`
  - Run data directory. Default log directory is created under here (see “Logs” below).
- `listing_file: SAT_OBS/C03__listing.csv`
  - Observation “listing” file.
  - Because this path is relative, the file must exist at:
    - `<base_path>/SAT_OBS/C03__listing.csv`
- `env_python: miniconda3/bin/python3.13`
  - Python path used by parts of the workflow that launch Python externally.
  - Because it is relative, it refers to:
    - `<base_path>/miniconda3/bin/python3.13`

Other example entries point to scripts/templates used by the workflow:

- `run_submit_model_template`
- `path_submit_bsh`
- `path_filter`
- `chimere_par_template`
- `run_submit_replace_perturbations`
- `path_perturbed_emi`
- `chimere_input_emissions_dir`, `chimere_input_atm_dir`, `chimere_input_ibc_dir`

Operationally, you should ensure:

- each referenced file exists at the expected absolute location
- the run user has read permissions for input folders and write permissions for `path_data`

### `cluster`

Controls how jobs are submitted/managed (scheduler and queue/project identifiers).

Example values:

- `scheduler: SLURM`
- `cluster_queue: g100_usr_prod`
- `project_name: mimesi`

Adjust these to match your HPC environment.

### `logging`

Controls log verbosity and message format.

- `level: INFO`
- `format: ...`

### `satellite_data`

Observation search window setting used by the workflow:

- `search_window_seconds: 1800`

### `cleanup`

The example enables output/input cleanup to reduce disk usage:

- `enabled: true`
- `delete_window_*`: remove per-window inputs after they are no longer needed
- `retain_cycles`, `retention_interval_cycles`: keep only the most recent completed windows (“cycles”)
- `trim_end`, `trim_out` and `keep_*_vars`: write reduced NetCDF outputs with a minimal set of variables
- `keep_daily_first_hour_end`, `keep_daily_last_hour_end`: preserve some daily boundary outputs

Operational guidance:

- Start with cleanup **disabled** if you are debugging a new configuration.
- Enable cleanup once you are confident outputs are correct and disk usage needs to be controlled.

### `dart`

This section provides paths for DART-related directories used by the workflow:

- `work_dir`
- `analysis_dir`
- `preassim_dir`
- `posteriors_dir`
- `obs_converters_dir`

These are paths within your run environment and must match how your external DART installation and run directories are arranged.

## Outputs and logs

When you run `python -u main.py -c run_config.yaml`, a log file is created automatically.

Default behavior:

- logs are written under `paths.path_data` in a folder named `mimesi_orchestrator_logs`
- log filename includes the model type and a timestamp

See [Quickstart](quickstart.md) for basic log guidance and example monitoring commands.
