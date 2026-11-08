"""Backend registry"""

# rq-17848072

from __future__ import annotations

import logging
from importlib.metadata import entry_points

from ._errors import UnknownForceFieldError
from ._model import ForceFieldBackend

log = logging.getLogger(__name__)
ENTRY_POINT_GROUP = "zpoly.forcefields"

_BACKENDS: dict[str, ForceFieldBackend] = {}
_plugins_loaded = False


def _key(name: str) -> str:
    return name.replace("-", "").replace("_", "").casefold()


def _load_plugins() -> None:
    global _plugins_loaded
    if _plugins_loaded:
        return
    _plugins_loaded = True
    for ep in entry_points(group=ENTRY_POINT_GROUP):
        try:
            backend = ep.load()
            register_backend(backend() if isinstance(backend, type) else backend)
        except Exception as exc:  # a broken plugin must not break the built-ins
            log.warning("Could not load force-field plugin %s: %s", ep.name, exc)


# rq-adb49473
def register_backend(backend: ForceFieldBackend) -> None:
    """Register (or replace) a backend under its name (case-insensitive)"""
    if not isinstance(backend, ForceFieldBackend):
        raise TypeError("A backend needs name, check_available() and parametrize().")
    _BACKENDS[_key(backend.name)] = backend


def get_backend(name: str) -> ForceFieldBackend:
    _load_plugins()
    try:
        return _BACKENDS[_key(name)]
    except KeyError:
        raise UnknownForceFieldError(name, available_backends()) from None


def available_backends() -> list[str]:
    _load_plugins()
    return sorted(b.name for b in _BACKENDS.values())
