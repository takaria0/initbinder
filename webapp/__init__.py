"""Web UI backend package for InitBinder pipeline orchestration."""

from importlib.metadata import version, PackageNotFoundError

__all__ = ["__version__"]

try:
    __version__ = version("initbinder-webapp")  # optional, fallback below
except PackageNotFoundError:  # pragma: no cover - local dev
    __version__ = "0.1.0"

