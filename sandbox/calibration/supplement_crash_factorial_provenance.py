"""Record post-run source provenance for a 1991 crash-mechanism screen.

This supplement records revisions, dirty flags and source hashes without
modifying immutable simulation outputs. Earlier Julia runs sometimes recorded
Git revisions as `unknown` in this Windows environment.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import tomllib
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CORE = ROOT.parent / "COTSMod.jl"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git(repo: Path, *args: str) -> str:
    return subprocess.check_output(
        ["git", "-C", str(repo), *args], text=True, stderr=subprocess.PIPE
    ).strip()


def main(run: Path) -> None:
    output = run / "provenance_supplement.json"
    if output.exists():
        raise FileExistsError(output)
    with (run / "metadata.toml").open("rb") as stream:
        metadata = tomllib.load(stream)
    if metadata.get("experiment") not in ("crash_factorial", "lagged_hazard"):
        raise ValueError("Not a supported 1991 crash-mechanism run")
    for name, expected in metadata["outputs"].items():
        if sha256(run / name) != expected:
            raise ValueError(f"Run output changed after simulation: {name}")
    for name in ("script", "adapter", "scenario", "core"):
        path = {
            "script": ROOT / "sandbox/calibration/run_1991_initialization_test.jl",
            "adapter": ROOT / "ADRIA/src/ecosystem/cots.jl",
            "scenario": ROOT / "ADRIA/src/scenario.jl",
            "core": CORE / "src/COTSMod.jl",
        }[name]
        if name in metadata["inputs"] and sha256(path) != metadata["inputs"][name]:
            raise ValueError(f"Runtime source changed since the simulation: {name}")
    files = {
        "cotsmod_core": CORE / "src/COTSMod.jl",
        "adria_adapter": ROOT / "ADRIA/src/ecosystem/cots.jl",
        "adria_scenario": ROOT / "ADRIA/src/scenario.jl",
        "run_script": ROOT / "sandbox/calibration/run_1991_initialization_test.jl",
        "score": ROOT / "sandbox/calibration/cots_cycle_metrics.jl",
        "reef_area_mapping": ROOT / "sandbox/calibration/reef_observation_mapping.jl",
        "audit": ROOT / "sandbox/calibration" / (
            "audit_1991_lagged_hazard.py" if metadata["experiment"] == "lagged_hazard"
            else "audit_1991_crash_factorial.py"),
        "project": ROOT / "sandbox/Project.toml",
        "manifest": ROOT / "sandbox/Manifest.toml",
        "protocol": ROOT / "sandbox/calibration/crash_mechanism_protocol.md",
    }
    report = {
        "status": "post_run_provenance_supplement_not_clean_revision_replay",
        "captured_utc": datetime.now(timezone.utc).isoformat(),
        "note": "Captured after simulation and audit. Available runtime script/adapter/scenario/core hashes were checked against simulation metadata. Both repositories may contain uncommitted changes, so revisions alone do not reproduce this run.",
        "run_metadata_sha256": sha256(run / "metadata.toml"),
        "audited_results_sha256": sha256(run / (
            "lagged_hazard_audit.json" if metadata["experiment"] == "lagged_hazard"
            else "crash_factorial_audit.json")),
        "adria_revision": git(ROOT, "rev-parse", "HEAD"),
        "cotsmod_revision": git(CORE, "rev-parse", "HEAD"),
        "adria_dirty": bool(git(ROOT, "status", "--porcelain")),
        "cotsmod_dirty": bool(git(CORE, "status", "--porcelain")),
        "source_sha256": {name: sha256(path) for name, path in files.items()},
        "python_version": sys.version.split()[0],
        "script_sha256": sha256(Path(__file__)),
    }
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    main(parser.parse_args().run)
