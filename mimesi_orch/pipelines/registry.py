"""Pipeline registry and factory helpers."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Dict, Iterable

from pipelines.contracts import BasePipeline


@dataclass(frozen=True)
class RegisteredPipeline:
    name: str
    build: Callable[[object, object], BasePipeline]


_REGISTRY: Dict[str, RegisteredPipeline] = {}


def register_pipeline(name: str, build: Callable[[object, object], BasePipeline]) -> None:
    key = name.strip().lower()
    if key in _REGISTRY:
        raise ValueError(f"Pipeline already registered: {name}")
    _REGISTRY[key] = RegisteredPipeline(name=key, build=build)


def get_registered_pipeline_names() -> Iterable[str]:
    return _REGISTRY.keys()


def create_pipeline(name: str, config: object, time_manager: object) -> BasePipeline:
    key = name.strip().lower()
    if key not in _REGISTRY:
        available = ", ".join(sorted(_REGISTRY))
        raise ValueError(f"Unknown pipeline '{name}'. Available: {available}")
    return _REGISTRY[key].build(config, time_manager)
