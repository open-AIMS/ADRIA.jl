"""Recover and verify metadata for a completed seed-boundary screen.

The first 2026-10-07 run completed all six scenarios but the Julia process
could not serialize a Matrix of treatment dictionaries as TOML. This script
does not rerun or alter trajectories; it verifies the saved control against
the archived run and writes a complete, non-overwriting JSON manifest.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import tomllib
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / "sandbox/calibration/runs"
FILES = ("scores.csv", "trajectories.csv", "flows.csv", "initialization.csv", "failures.csv")
TREATMENTS = {f"{seed}_{boundary}" for seed in ("reference", "quarter")
              for boundary in ("full", "10pct", "zero")}


def digest(path: Path) -> str:
    sha = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            sha.update(block)
    return sha.hexdigest()


def revision(path: Path) -> str:
    return subprocess.check_output(["git", "-C", str(path), "rev-parse", "HEAD"],
                                   text=True).strip()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-id", required=True)
    args = parser.parse_args()
    if not args.run_id.replace("_", "").replace("-", "").isalnum():
        raise ValueError("Invalid run ID")
    run = RUNS / args.run_id
    output = run / "metadata.json"
    if output.exists():
        raise FileExistsError(output)
    scores = pd.read_csv(run / "scores.csv")
    traj = pd.read_csv(run / "trajectories.csv")
    flows = pd.read_csv(run / "flows.csv")
    init = pd.read_csv(run / "initialization.csv")
    if (run / "failures.csv").read_text(encoding="utf-8").strip() != "treatment,error":
        raise ValueError("Run contains failed evaluations")
    if (set(scores.treatment) != TREATMENTS or set(traj.treatment) != TREATMENTS or
            set(flows.treatment) != TREATMENTS or set(init.treatment) != TREATMENTS or
            len(scores) != 48 or len(traj) != 960 or len(flows) != 960 or len(init) != 24):
        raise ValueError("Six-cell screen incomplete")
    archived_path = RUNS / "20261001T163000_owen_embedded_replicates/trajectories.csv"
    archived = pd.read_csv(archived_path)
    archived = archived[(archived.seed == 20260930) &
                        (archived.treatment == "embedded_joint_4x")]
    reference = traj[traj.treatment == "reference_full"]
    paired = reference.merge(archived, on=["reef_name", "year"], suffixes=("_new", "_old"),
                             validate="one_to_one")
    if len(paired) != 160 or not np.array_equal(
            paired.forcing_year_new.to_numpy(), paired.forcing_year_old.to_numpy()):
        raise ValueError("Reference forcing years differ")
    adult_error = float(np.max(np.abs(paired.adults_ha_new - paired.adults_ha_old)))
    coral_error = float(np.max(np.abs(paired.coral_cover_new - paired.coral_cover_old)))
    if adult_error > 1e-12 or coral_error > 1e-12:
        raise ValueError("Reference replay changed")

    screen_path = RUNS / "20261001T162000_owen_embedded_mortality_screen/metadata.toml"
    screen = tomllib.loads(screen_path.read_text(encoding="utf-8"))
    candidate_path = RUNS / "expanded_cotsconn_pilot_seed20260930/best_summary.csv"
    candidate = pd.read_csv(candidate_path).sort_values("loss").iloc[0]
    full_boundary = 4 * float(screen["matched_boundary_pre_pelagic_ha_year"])
    fecundity = 4 * float(screen["matched_fecundity_pre_pelagic_ha_year"])
    provenance = ROOT / ("sandbox/data/Lizard_Historical_v0.2/cots_connectivity/datasets/"
                         "owen_global_2018_2023/provenance.json")
    sources = {
        "archived_trajectories": archived_path,
        "candidate": candidate_path,
        "gate_scores": RUNS / "20261001T113401_lizard_domain_gate/scores.csv",
        "screen_metadata": screen_path,
        "observations": ROOT / "sandbox/data/reef_cots.csv",
        "owen_provenance": provenance,
        "protocol": ROOT / "sandbox/calibration/CALIBRATION_PROTOCOL.md",
        "screen_script_corrected_serialization": ROOT / "sandbox/calibration/screen_owen_seed_boundary.jl",
        "adapter": ROOT / "ADRIA/src/ecosystem/cots.jl",
        "scenario": ROOT / "ADRIA/src/scenario.jl",
        "core": ROOT / "../COTSMod.jl/src/COTSMod.jl",
        "finalizer": Path(__file__),
    }
    metadata = {
        "status": "completed_bounded_diagnostic_not_promoted",
        "seed": 20260930,
        "model_years": [1985, 2024],
        "hydrodynamic_years": [2018, 2019, 2020, 2021, 2022],
        "hydrodynamic_mode": "sample",
        "hydrodynamic_seed": 20260930,
        "connectivity_dataset": "owen_global_2018_2023",
        "connectivity_orientation": "source rows, sink columns",
        "source_survival_applied_separately": False,
        "fixed_allee_adults_ha": 3.0,
        "fixed_scoring": "three-year-smoothed observed CPUE; V1 frozen reef-specific scales",
        "seed_multiplier_reference": float(candidate.seed_mult),
        "seed_multiplier_quarter": float(candidate.seed_mult) * 0.25,
        "seed_stage_ratio_juvenile_subadult_adult": [6, 3, 1],
        "fecundity_pre_pelagic_recruits_per_adult_year": fecundity,
        "boundary_full_pre_pelagic_recruits_per_ha_year": full_boundary,
        "boundary_fractions": [1.0, 0.1, 0.0],
        "first_cots_log_row_caveat": "1985 is unpopulated; initial state reconstructed separately",
        "reference_max_adult_difference_ha": adult_error,
        "reference_max_coral_difference_fraction": coral_error,
        "source_script_note": "After simulations, only treatment-list vectorization for TOML was patched; numeric logic did not change",
        "adria_revision": revision(ROOT),
        "cotsmod_revision": revision(ROOT / "../COTSMod.jl"),
        "source_sha256": {name: digest(path) for name, path in sources.items()},
        "output_sha256": {name: digest(run / name) for name in FILES},
    }
    output.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print("Reference max differences: adults", adult_error, "coral", coral_error)
    print("Wrote", output)


if __name__ == "__main__":
    main()
