"""Verified execution and publication of atlas extraction artifacts."""

import hashlib
import json
import os
import platform
import subprocess
import tempfile
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path

import h5py

from hypomap_exploration.extraction import extract_atlas_counts


def sha256_file(path):
    """Hash a file sequentially without loading it into memory."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _git_state(root):
    """Record Git state, explicitly allowing non-repository test fixtures."""
    try:
        commit = subprocess.run(
            ["git", "-C", str(root), "rev-parse", "HEAD"],
            check=True, capture_output=True, text=True,
        ).stdout.strip()
        status = subprocess.run(
            ["git", "-C", str(root), "status", "--porcelain"],
            check=True, capture_output=True, text=True,
        ).stdout
        return {
            "commit": commit,
            "dirty": bool(status.strip()),
            "status": status.splitlines(),
        }
    except (FileNotFoundError, subprocess.CalledProcessError):
        return {"commit": None, "dirty": None, "status": None}


def _publish_without_overwrite(temporary, destination):
    """Publish a complete file without replacing an existing destination."""
    os.link(temporary, destination)


def run_extraction(config_path, project_root, *, max_nnz_per_chunk):
    """Verify inputs, extract counts and publish a table with provenance."""
    root = Path(project_root).resolve()
    config_path = Path(config_path)
    if not config_path.is_absolute():
        config_path = root / config_path

    output = root / "data/derived/all_nuclei_gene_counts.tsv.gz"
    manifest_path = root / "provenance/extraction_manifest.json"

    for destination in (output, manifest_path):
        if destination.exists() or destination.is_symlink():
            raise FileExistsError(
                f"Refusing to overwrite existing artifact: {destination}"
            )

    config_bytes = config_path.read_bytes()
    config = json.loads(config_bytes)
    config_checksum = hashlib.sha256(config_bytes).hexdigest()

    source = Path(config["source_file"])
    if not source.is_absolute():
        source = root / source

    expected = config.get("source_sha256")
    if (
        not isinstance(expected, str)
        or len(expected) != 64
        or any(c not in "0123456789abcdefABCDEF" for c in expected)
    ):
        raise ValueError("Source checksum must be a 64-character SHA-256")

    source_checksum = sha256_file(source)
    if source_checksum != expected.lower():
        raise ValueError(
            f"Source checksum mismatch: expected {expected}, "
            f"found {source_checksum}"
        )

    git_state = _git_state(root)
    with h5py.File(source, "r") as handle:
        table, mapping = extract_atlas_counts(
            handle, config, max_nnz_per_chunk=max_nnz_per_chunk
        )

    if len(table) != config["extraction"]["expected_nuclei"]:
        raise ValueError("Extracted row count differs from configuration")

    missing_marker = "__MISSING__"
    metadata_fields = config["extraction"]["preserve_obs_fields"]
    text_fields = ["cell_id", *metadata_fields]
    for field in text_fields:
        if table[field].eq(missing_marker).any():
            raise ValueError(
                f"Reserved missing-value marker occurs in column {field}"
            )

    missing_counts = {
        field: int(table[field].isna().sum())
        for field in metadata_fields
    }

    temporary_paths = []
    published_paths = []

    try:
        output.parent.mkdir(parents=True, exist_ok=True)
        manifest_path.parent.mkdir(parents=True, exist_ok=True)

        def temporary_file(directory, suffix):
            descriptor, name = tempfile.mkstemp(
                prefix=".extraction-", suffix=suffix, dir=directory
            )
            os.close(descriptor)
            path = Path(name)
            temporary_paths.append(path)
            return path

        temporary_table = temporary_file(output.parent, ".tsv.gz")
        table.to_csv(
            temporary_table,
            sep="\t",
            index=False,
            na_rep=missing_marker,
            compression={"method": "gzip", "mtime": 0},
        )

        manifest = {
            "status": "complete",
            "source_file": config["source_file"],
            "source_sha256": source_checksum,
            "config_sha256": config_checksum,
            "output_file": str(output.relative_to(root)),
            "output_sha256": sha256_file(temporary_table),
            "n_nuclei": len(table),
            "n_columns": len(table.columns),
            "gene_columns": mapping,
            "missing_metadata_counts": missing_counts,
            "table_read_options": {
                "sep": "\t",
                "keep_default_na": False,
                "na_values": [missing_marker],
            },
            "git": git_state,
            "software_versions": {
                "python": platform.python_version(),
                **{
                    package: version(package)
                    for package in ("numpy", "pandas", "h5py")
                },
            },
            "max_nnz_per_chunk": max_nnz_per_chunk,
            "validation": {
                "source_checksum_matches": True,
                "wrapper_and_csr_validation": "passed",
                "expected_row_count_matches": True,
                "reserved_missing_marker_absent": True,
            },
            "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        }

        temporary_manifest = temporary_file(
            manifest_path.parent, ".json"
        )
        temporary_manifest.write_text(
            json.dumps(manifest, indent=2, allow_nan=False) + "\n"
        )

        _publish_without_overwrite(temporary_table, output)
        published_paths.append(output)

        _publish_without_overwrite(temporary_manifest, manifest_path)
        published_paths.append(manifest_path)

        return manifest

    except Exception:
        for path in reversed(published_paths):
            path.unlink(missing_ok=True)
        raise

    finally:
        for path in temporary_paths:
            path.unlink(missing_ok=True)
