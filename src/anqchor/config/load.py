"""Load, merge, validate and resolve configuration files.

Precedence (highest first): command-line flags, the YAML file, model defaults.
All relative paths are resolved against the directory of the YAML file, or against
the current directory for flag-only invocations.
"""

from __future__ import annotations

import copy
from pathlib import Path
from typing import Any

import yaml
from pydantic import BaseModel, ValidationError

from anqchor.config.model import FILE_PREFIX, REFERENCE_PREFIX, Config


class ConfigError(Exception):
    """Invalid or inconsistent configuration (exit code 2)."""


def deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    """Recursively merge ``override`` into a copy of ``base``; ``None`` values in
    ``override`` are skipped so unset CLI flags never clobber the file."""
    out = copy.deepcopy(base)
    for key, value in override.items():
        if value is None:
            continue
        if isinstance(value, dict):
            # Recurse even when the base has no mapping here, so that None values
            # nested in a new section are dropped as well.
            base_value = out.get(key)
            out[key] = deep_merge(base_value if isinstance(base_value, dict) else {}, value)
        else:
            out[key] = copy.deepcopy(value)
    return out


def read_yaml(path: Path) -> dict[str, Any]:
    """Parse a YAML mapping; raises :class:`ConfigError` for a missing file or bad YAML."""
    try:
        with open(path, encoding="utf-8") as fh:
            data = yaml.safe_load(fh)
    except FileNotFoundError as e:
        raise ConfigError(f"configuration file not found: {path}") from e
    except OSError as e:
        raise ConfigError(f"cannot read configuration file {path}: {e}") from e
    except yaml.YAMLError as e:
        raise ConfigError(f"cannot parse YAML in {path}: {e}") from e
    if data is None:
        return {}
    if not isinstance(data, dict):
        raise ConfigError(f"{path}: the top level must be a mapping")
    return data


def format_validation_error(err: ValidationError) -> str:
    """One line per pydantic error: ``  dotted.location: message``."""
    lines = []
    for e in err.errors():
        loc = ".".join(str(p) for p in e["loc"]) or "(root)"
        lines.append(f"  {loc}: {e['msg']}")
    return "invalid configuration:\n" + "\n".join(lines)


def load_config(
    path: Path | None = None,
    overrides: dict[str, Any] | None = None,
    base_dir: Path | None = None,
) -> Config:
    """Build a validated :class:`Config` with absolute paths.

    Args:
        path: YAML file (optional when ``overrides`` describe everything).
        overrides: nested dict of values from CLI flags, merged over the file.
        base_dir: directory for resolving relative paths in the file; defaults to the
            file's directory, or the current directory without a file. Paths given in
            ``overrides`` are resolved against the current directory before merging.
    """
    raw: dict[str, Any] = {}
    if path is not None:
        path = Path(path).expanduser().resolve()
        raw = read_yaml(path)
        if base_dir is None:
            base_dir = path.parent
    if base_dir is None:
        base_dir = Path.cwd()
    # The file's relative paths resolve against base_dir; flag values resolve
    # against the current directory. Both are made absolute before merging.
    try:
        merged = _walk_raw(copy.deepcopy(raw), base_dir, ())
        if overrides:
            ov = _walk_raw(copy.deepcopy(overrides), Path.cwd(), ())
            # A genotype flag replaces the file's genotype section: 'pgen_list' and
            # 'pfile' are alternatives, so merging them would give both.
            if isinstance(ov.get("genotypes"), dict) and ov["genotypes"]:
                merged.pop("genotypes", None)
            merged = deep_merge(merged, ov)
        cfg = Config.model_validate(merged)
    except ValidationError as e:
        raise ConfigError(format_validation_error(e)) from e
    return resolve_paths(cfg, base_dir)


_PATH_KEYS = {
    "binary",
    "pgen_list",
    "pfile",
    "pairs",
    "variants",
    "pheno",
    "trait_annot",
    "cov",
    "sample_map",
    "trait_map",
    "prs_samples",
    "file",
    "out_dir",
    "traits",
}
_NON_PATH_LITERALS = {"default", "reference", "all"}


def _abs(s: str, base: Path) -> str:
    p = Path(s).expanduser()
    return str(p if p.is_absolute() else (base / p).resolve())


def _walk_raw(node: Any, base: Path, trail: tuple[str, ...]) -> Any:
    """Make path-valued entries of a raw (YAML or flag) mapping absolute against ``base``.

    Entries are recognised by key (:data:`_PATH_KEYS`). Literals that are not paths
    (``default``, ``reference``, ``all``, ``reference:<name>``) are kept; a ``file:``
    prefix is removed and the rest resolved. ``trait_id_rule: file:<path>`` keeps its
    prefix (the model turns it into ``trait_map``) but its path is resolved here.
    """
    if isinstance(node, dict):
        return {k: _walk_raw(v, base, (*trail, str(k))) for k, v in node.items()}
    if isinstance(node, list):
        return [_walk_raw(v, base, trail) for v in node]
    if not trail or not isinstance(node, (str, Path)):
        return node
    key = trail[-1]
    s = str(node)
    if key == "trait_id_rule" and s.startswith(FILE_PREFIX):
        return FILE_PREFIX + _abs(s[len(FILE_PREFIX) :], base)
    if key in _PATH_KEYS:
        if s in _NON_PATH_LITERALS or s.startswith(REFERENCE_PREFIX):
            return s
        if s.startswith(FILE_PREFIX):
            s = s[len(FILE_PREFIX) :]
        return _abs(s, base)
    return node


def resolve_paths(cfg: Config, base_dir: Path) -> Config:
    """Return a copy of ``cfg`` with every :class:`~pathlib.Path` field made absolute."""
    base_dir = Path(base_dir).resolve()

    def walk(obj: Any) -> Any:
        if isinstance(obj, BaseModel):
            updates = {}
            for name in type(obj).model_fields:
                value = getattr(obj, name)
                new = walk(value)
                if new is not value:
                    updates[name] = new
            return obj.model_copy(update=updates) if updates else obj
        if isinstance(obj, Path):
            p = obj.expanduser()
            return p if p.is_absolute() else (base_dir / p).resolve()
        if isinstance(obj, dict):
            return {k: walk(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [walk(v) for v in obj]
        return obj

    out = walk(cfg)
    assert isinstance(out, Config)
    return out


def config_to_dict(cfg: Config) -> dict[str, Any]:
    """JSON/YAML-friendly dict (paths as strings, aliases such as ``lambda`` used)."""
    data: dict[str, Any] = cfg.model_dump(mode="json", by_alias=True)
    return data


def write_resolved_config(cfg: Config, path: Path) -> None:
    """Write the fully resolved configuration (defaults filled, paths absolute)."""
    from anqchor.util.files import write_text_if_changed

    text = "# Resolved by anqchor: defaults filled in, paths absolute.\n" + yaml.safe_dump(
        config_to_dict(cfg), sort_keys=False, default_flow_style=False
    )
    write_text_if_changed(path, text)


def load_resolved_config(run_dir: Path) -> Config:
    """Read ``<run_dir>/anqchor.resolved.yaml`` (already absolute)."""
    raw = read_yaml(Path(run_dir) / "anqchor.resolved.yaml")
    try:
        return Config.model_validate(raw)
    except ValidationError as e:  # pragma: no cover - written by us
        raise ConfigError(format_validation_error(e)) from e
