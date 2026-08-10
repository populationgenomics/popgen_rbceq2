"""Helper functions for the tests."""

from typing import Any, Protocol, runtime_checkable

import toml
from cpg_utils import Path
from cpg_utils.config import set_config_paths


@runtime_checkable
class IDictRepresentable(Protocol):
    def as_dict(self) -> dict[str, Any]: ...


class TomlAnyPathEncoder(toml.TomlEncoder):
    """Support for CPG path objects in TOML.

    A CPG path is either a regular pathlib path or a cloud path.
    """

    def dump_value(self, v):
        if isinstance(v, Path):
            v = str(v)
        return super().dump_value(v)


def set_config(
    config: str | dict[str, Any] | IDictRepresentable,
    path: Path,
    merge_with: list[Path] | None = None,
) -> None:
    """Write a config to `path` and point `CPG_CONFIG_PATH` at it.

    That environment variable is how `cpg_utils` finds a config. If `merge_with` is provided,
    the config is merged with the configs at the given paths. Merging happens right to left, so
    values in the right config override values in the left one.

    Args:
        config (str | dict[str, Any] | IDictRepresentable):
            A valid TOML string, a dictionary to be converted to TOML, or an object
            which implements the `IDictRepresentable` protocol.

        path (Path):
            Path to write the config to.

        merge_with (list[Path] | None, optional):
            A list of paths to merge with the config. Merging happens right to left,
            so that values in the right config will override values in the left config.
            Defaults to `None`.
    """
    with path.open('w') as f:
        if isinstance(config, dict):
            toml.dump(config, f, encoder=TomlAnyPathEncoder())
        elif isinstance(config, IDictRepresentable):
            toml.dump(config.as_dict(), f, encoder=TomlAnyPathEncoder())
        elif isinstance(config, str):
            f.write(config)
        else:
            raise TypeError(f'Expected config to be a string, dict, or IDictRepresentable, butgot {type(config)}')

        f.flush()

    return set_config_paths([*[str(s) for s in (merge_with or [])], str(path)])
