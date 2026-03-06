# Pipeline Abstraction Migration

This folder now contains a pipeline abstraction layer that isolates model-specific code.

## Layout

- `contracts.py`: base interfaces (`BasePipeline`, `PipelineBuilder`)
- `registry.py`: pipeline registration and lookup
- `factory.py`: pipeline selection from config
- `farm/`, `chimere2017/`, `chimere2023/`: pipeline-specific `config.py`, `paths.py`, `pipeline.py`

## Backward compatibility

- Existing execution logic remains in:
  - `farm_pipeline.py`
  - `chimere_pipeline.py`
  - `chimere_irene_pipeline.py`
- New pipeline modules wrap those classes and are instantiated through `factory.build_pipeline(...)`.
- `config_models.AppConfig` remains valid for legacy YAMLs.

## Current status

- CHIMERE 2017:
  - Logic copied from legacy `chimere_pipeline.py` into `chimere2017/pipeline.py`.
  - Pipeline now uses `chimere2017/paths.py` directly.
  - `chimere2017/paths.py` no longer depends on global `PathManager`.
- FARM and CHIMERE 2023:
  - Still use wrappers over legacy implementations and `PathManager`-backed path facades.

## Next migration steps

1. Apply the same migration pattern to `farm`:
   - copy logic into `pipelines/farm/pipeline.py`
   - remove `PathManager` coupling from `pipelines/farm/paths.py`.
2. Apply the same migration pattern to `chimere2023`:
   - copy logic into `pipelines/chimere2023/pipeline.py`
   - remove `PathManager` coupling from `pipelines/chimere2023/paths.py`.
3. After all pipelines are migrated, convert shared helpers in root `paths.py` into model-agnostic utilities only.
4. Remove legacy files (`farm_pipeline.py`, `chimere_pipeline.py`, `chimere_irene_pipeline.py`) after import checks and transition period.
5. Optionally migrate to explicit config section:

```yaml
pipeline:
  name: chimere2023
```

If `pipeline.name` is absent, `factory.py` infers the pipeline from existing fields.
