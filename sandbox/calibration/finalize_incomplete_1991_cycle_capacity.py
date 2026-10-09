"""Write immutable provenance for an interrupted 1991 capacity screen."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
CAL = ROOT / "sandbox/calibration"
CORE = ROOT.parent / "COTSMod.jl"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def revision(path: Path) -> str:
    return subprocess.check_output(["git", "-C", str(path), "rev-parse", "HEAD"],
                                   text=True).strip()


def main(run: Path, reason: str) -> None:
    output = run / "incomplete_provenance.json"
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite {output}")
    design = pd.read_csv(run / "design.csv")
    scores = pd.read_csv(run / "scores.csv")
    failures = pd.read_csv(run / "failures.csv")
    expected = {f"pair_{int(pair):03d}_theta_{theta}"
                for pair in design.pair for theta in (3, 1)}
    completed = set(scores.treatment)
    bad_groups = scores.groupby("treatment").size()
    if (bad_groups != 4).any():
        raise ValueError("A treatment has incomplete per-reef scores")
    missing = sorted(expected - completed)
    if not missing:
        raise ValueError("Run is complete; use normal Julia metadata")
    inputs = {
        "1991_surface": CAL / "runs/20261007T199104_1991_initial_surfaces/surface.csv",
        "1991_surface_provenance": CAL / "runs/20261007T199104_1991_initial_surfaces/metadata.json",
        "1991_reference_trajectories": CAL / "runs/20261007T199111_owen_1991_cohort_history_final/trajectories.csv",
        "reference_parameters": CAL / "runs/expanded_cotsconn_pilot_seed20260930/best_summary.csv",
        "screen_metadata": CAL / "runs/20261001T162000_owen_embedded_mortality_screen/metadata.toml",
        "observations": ROOT / "sandbox/data/reef_cots.csv",
        "manifest": ROOT / "sandbox/Manifest.toml",
        "owen_forcing": ROOT / "sandbox/data/Lizard_Historical_v0.2/cots_connectivity/datasets/owen_global_2018_2023/provenance.json",
        "runner": CAL / "screen_1991_cycle_capacity.jl",
        "audit": CAL / "audit_1991_cycle_capacity.py",
        "adapter": ROOT / "ADRIA/src/ecosystem/cots.jl",
        "core": CORE / "src/COTSMod.jl",
    }
    outputs = [run / name for name in
               ("design.csv", "scores.csv", "trajectories.csv", "flows.csv",
                "failures.csv", "capacity_audit.csv",
                "capacity_treatment_summary.csv", "capacity_audit.json")]
    outputs.extend(sorted(run.glob("initial_state_pair_*.csv")))
    info = {
        "status": "interrupted_incomplete_no_promotion",
        "reason": reason,
        "completed_treatments": len(completed),
        "expected_treatments": len(expected),
        "missing_treatments": missing,
        "failed_evaluations_recorded": len(failures),
        "design_pairs_including_control": len(design),
        "design_seed": 20261008,
        "hydrodynamic_seed": 20260930,
        "julia_version": "1.11.4",
        "calibration_allee_adults_ha": 3.0,
        "counterfactual_allee_adults_ha": 1.0,
        "orientation": "Owen source rows, sink columns",
        "larval_survival_assumption": "embedded in Owen connectivity; no extra modifier",
        "q_units": "COTS/tow per adult COTS/ha, same factor for initial state and observation",
        "year_window": [1991, 2024],
        "score_start_year": 1992,
        "adria_revision": revision(ROOT),
        "cotsmod_revision": revision(CORE),
        "inputs_sha256": {name: sha256(path) for name, path in inputs.items()},
        "outputs_sha256": {path.name: sha256(path) for path in outputs},
        "finalizer_sha256": sha256(Path(__file__)),
    }
    output.write_text(json.dumps(info, indent=2, sort_keys=True) + "\n")
    print(output)
    print(f"Completed {len(completed)}/{len(expected)}; missing: {missing}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("--reason", required=True)
    args = parser.parse_args()
    main(args.run, args.reason)
