# Troubleshooting

## Where to look first

1. Find the log file for your run (see [Quickstart](quickstart.md)).
2. Search for the first `ERROR` / `CRITICAL` message in the log.
3. If the program stops immediately, it is usually a configuration or input-path problem.

## Startup fails with configuration errors

The orchestrator validates the YAML config strictly. Common symptoms:

- errors about missing required sections/fields (e.g. `paths`, `time`, `assimilation`, `cluster`)
- errors about unexpected keys (legacy configs)
- errors about invalid values (e.g. `dt_seconds` must be `> 0`)

If you are using `config_orchestrator.yaml` as a starting point, be aware it appears to be a legacy example and may require updates to match current field names in `config_models.py`.

## Startup fails with “only chimere2017 is enabled”

This repository only supports the CHIMERE 2017 workflow. In practice:

- if you set `pipeline.name`, it must be `chimere2017`
- `assimilation.model_type` must be `chimere`

If you see an error like "Only 'chimere2017' is enabled", adjust the config accordingly.

## Assimilation is skipped (but the run continues)

Skips can be legitimate. For the enabled workflow, assimilation depends on finding matching entries in the observation listing file near the end of a time window. If there are no matches, the program can run model windows without assimilation.

Checks:

- confirm `paths.listing_file` points to the expected file
- confirm it contains `start_time` values in the time period you are running
- confirm the file is `;`-separated and readable

## Missing external tools

If your run exercises code paths that call external tools (for example `ncks` or `cdo` in `pipelines/chimere2017/pipeline.py`) you may see failures if those commands are not available in your environment.

Operational checks:

```bash
command -v ncks || echo "ncks not found"
command -v cdo  || echo "cdo not found"
```

## Scheduler / submission failures

Scheduler-related failures are represented by `pipeline_errors.SchedulerError` / `pipeline_errors.ModelRunError`. The exact submission commands and checks are implemented in `orchestrator_utils.py` and `scheduler.py`, and depend on your runtime environment.

If you see scheduler errors:

- confirm `cluster.scheduler` matches your system (`slurm` or `lsf` are the available values in this repository)
- confirm `cluster.cluster_queue` and `cluster.project_name` are valid for your environment

## “File not found” / missing-path errors

Most file locations come from `paths.*` in your YAML config. If a file is missing:

1. Check whether the path is absolute or relative.
2. If relative, confirm `paths.base_path` is correct (relative paths are interpreted under it).
3. Confirm the referenced file/folder exists and has the expected permissions.
