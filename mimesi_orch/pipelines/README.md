# Pipeline Abstraction Migration

This folder now contains a CHIMERE 2017-only pipeline setup.

## Layout

- `contracts.py`: base interfaces (`BasePipeline`, `PipelineBuilder`)
- `registry.py`: pipeline registration and lookup
- `factory.py`: pipeline selection from config
- `chimere2017/`: active pipeline (`config.py`, `paths.py`, `pipeline.py`)

## Scope

- Enabled: CHIMERE 2017 pipeline.
- Disabled/removed from this branch:
  - FARM pipeline
  - CHIMERE 2023 pipeline
  - legacy wrapper pipeline modules

## Current status

- CHIMERE 2017:
  - Logic moved into `chimere2017/pipeline.py`.
  - Pipeline now uses `chimere2017/paths.py` directly.
  - `chimere2017/paths.py` no longer depends on global `PathManager`.

## Configuration notes

- `factory.py` now enforces CHIMERE 2017 only.
- If a config requests FARM or CHIMERE 2023-style settings, startup fails with an explicit error.
- Optional explicit selection:

```yaml
pipeline:
  name: chimere2017
```
