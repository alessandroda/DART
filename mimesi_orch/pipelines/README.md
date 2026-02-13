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

## Next migration steps

1. Move logic from legacy `*_pipeline.py` files into `pipelines/<name>/pipeline.py`.
2. Replace usage of legacy `PathManager` with pipeline-local path classes.
3. Remove mixed methods from `paths.py` once all pipelines are migrated.
4. Optionally migrate to explicit config section:

```yaml
pipeline:
  name: chimere2023
```

If `pipeline.name` is absent, `factory.py` infers the pipeline from existing fields.
