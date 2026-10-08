"""Audit Owen input provenance and reconstruct the spatial COTS seed.

The first COTS log row is unpopulated (the scenario loop starts at timestep 2),
so the actual initial state is reconstructed from the domain field and the
archived spatial-initialization formula. No source data or model state is edited.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / "sandbox/calibration/runs"
DOMAIN = ROOT / "sandbox/data/Lizard_Historical_v0.2"
DATASET = DOMAIN / "cots_connectivity/datasets/owen_global_2018_2023"
TARGETS = ("Lizard Island Reef", "MacGillivray Reef", "North Direction Reef", "Eyrie Reef")
ALLEE = 3.0  # COTS adults/ha; frozen for this diagnostic


def digest(path: Path) -> str:
    sha = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            sha.update(block)
    return sha.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-id", required=True)
    args = parser.parse_args()
    if not args.run_id.replace("_", "").replace("-", "").isalnum():
        raise ValueError("Invalid run ID")
    out = RUNS / args.run_id
    out.mkdir(exist_ok=False)

    provenance_path = DATASET / "provenance.json"
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    expected_dir = "sandbox/data/OwenHydro/connectivityMatsYearlyGlobal/"
    matrix_sources = [entry for entry in provenance["sources"]
                      if entry["path"].startswith(expected_dir + "conMatCots")]
    if len(matrix_sources) != 5 or provenance["orientation"] != (
        "rows are larval sources; columns are settlement sinks"
    ):
        raise ValueError("Unexpected Owen global input contract")
    input_rows = []
    for entry in provenance["sources"] + provenance["outputs"]:
        path = ROOT / entry["path"] if entry in provenance["sources"] else DATASET / entry["path"]
        if not path.is_file():
            raise FileNotFoundError(path)
        actual = digest(path)
        input_rows.append({"path": str(path.relative_to(ROOT)), "bytes": path.stat().st_size,
                           "sha256": actual, "provenance_sha256": entry["sha256"],
                           "matches_provenance": actual == entry["sha256"]})
    inputs = pd.DataFrame(input_rows)
    if not inputs.matches_provenance.all():
        raise ValueError("Owen source or derived output changed since build")
    inputs.to_csv(out / "input_hashes.csv", index=False)

    candidate_path = RUNS / "expanded_cotsconn_pilot_seed20260930/best_summary.csv"
    gate_path = RUNS / "20261001T113401_lizard_domain_gate/scores.csv"
    trajectory_path = RUNS / "20261002T210000_owen_low_cover_mortality_screen/trajectories.csv"
    candidate = pd.read_csv(candidate_path).sort_values("loss").iloc[0]
    multiplier = float(candidate.seed_mult)
    scales = pd.read_csv(gate_path)
    scales = scales[(scales.seed == 20260930) & (scales.treatment == "v1_legacy_closed") &
                    (scales.observation_treatment == "raw")].set_index("reef_name")
    trajectories = pd.read_csv(trajectory_path)
    trajectories = trajectories[trajectories.treatment == "joint4_reference"]
    sites = gpd.read_file(DOMAIN / "spatial/lizard_cluster.gpkg",
                          columns=["reef_name", "area", "cots_density"])
    if len(sites) != 2895 or len(provenance["hydrodynamic_years"]) != 5:
        raise ValueError("Unexpected site or hydrodynamic-year count")
    sites["model_reef"] = sites.reef_name.str.split(" (", regex=False).str[0]
    probability = sites.cots_density.to_numpy(dtype=float)
    if not np.all(np.isfinite(probability)):
        raise ValueError("Non-finite spatial COTS seed field")
    scale = np.where(probability > 0.1, np.minimum(probability / 0.5, 1.0), 0.0)
    sites["initial_juveniles_ha"] = 0.6 * scale * multiplier
    sites["initial_subadults_ha"] = 0.3 * scale * multiplier
    sites["initial_adults_ha"] = 0.1 * scale * multiplier
    sites["initial_allee_factor"] = sites.initial_adults_ha.pow(2) / (
        ALLEE ** 2 + sites.initial_adults_ha.pow(2))

    rows = []
    for reef in TARGETS:
        subset = sites[sites.model_reef == reef]
        if subset.empty or reef not in scales.index:
            raise ValueError(f"Missing mapped target or frozen scale: {reef}")
        weights = subset["area"].to_numpy(dtype=float)
        weights /= weights.sum()
        avg = lambda column: float(np.dot(weights, subset[column].to_numpy(dtype=float)))
        scale_cpue = float(scales.loc[reef, "observation_scale"])
        logged = trajectories[trajectories.reef_name == reef].set_index("year")
        if sorted(logged.index) != list(range(1985, 2025)):
            raise ValueError(f"Incomplete archived trajectory: {reef}")
        rows.append({
            "reef_name": reef, "site_count": len(subset), "seed_multiplier": multiplier,
            "initial_juveniles_ha": avg("initial_juveniles_ha"),
            "initial_subadults_ha": avg("initial_subadults_ha"),
            "initial_adults_ha": avg("initial_adults_ha"),
            "initial_total_ha": sum(avg(c) for c in (
                "initial_juveniles_ha", "initial_subadults_ha", "initial_adults_ha")),
            "initial_adult_above_allee_site_fraction": float(np.mean(
                subset.initial_adults_ha >= ALLEE)),
            "initial_mean_allee_factor": avg("initial_allee_factor"),
            "fixed_cpue_per_adult_ha": scale_cpue,
            "allee_equivalent_cpue": ALLEE * scale_cpue,
            "initial_adult_equivalent_cpue": avg("initial_adults_ha") * scale_cpue,
            "logged_1985_adults_ha": float(logged.loc[1985, "adults_ha"]),
            "modeled_1986_adults_ha": float(logged.loc[1986, "adults_ha"]),
            "modeled_1989_adults_ha": float(logged.loc[1989, "adults_ha"]),
            "modeled_1989_cpue": float(logged.loc[1989, "simulated_cpue"]),
            "adults_ha_per_0p1_cpue": 0.1 / scale_cpue,
            "adults_ha_per_0p4_cpue": 0.4 / scale_cpue,
        })
    summary = pd.DataFrame(rows)
    summary.to_csv(out / "initialization_by_reef.csv", index=False)

    metadata = {
        "status": "read_only_initialization_and_owen_source_audit",
        "initialization_equations": {
            "scale": "p>0.1 ? min(p/0.5,1) : 0",
            "juveniles_ha": "0.6 * scale * seed_multiplier",
            "subadults_ha": "0.3 * scale * seed_multiplier",
            "adults_ha": "0.1 * scale * seed_multiplier",
        },
        "log_first_year_caveat": "1985 cots_log row is zero-initialized; scenario loop begins at timestep 2",
        "allee_adults_ha": ALLEE,
        "hydrodynamic_years": provenance["hydrodynamic_years"],
        "matrix_paths": [entry["path"] for entry in matrix_sources],
        "source_to_sink_orientation": provenance["orientation"],
        "source_policy": provenance["policies"],
        "external_reef_count": 3692,
        "sources": {str(path.relative_to(ROOT)): digest(path) for path in
                    (provenance_path, candidate_path, gate_path, trajectory_path,
                     DOMAIN / "spatial/lizard_cluster.gpkg", Path(__file__))},
        "outputs": {name: digest(out / name) for name in
                    ("input_hashes.csv", "initialization_by_reef.csv")},
    }
    (out / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(summary[["reef_name", "initial_adults_ha", "initial_subadults_ha",
                   "initial_juveniles_ha", "initial_mean_allee_factor",
                   "allee_equivalent_cpue", "modeled_1989_cpue"]].round(4).to_string(index=False))
    print("Verified", len(inputs), "Owen source/output hashes; saved", out)


if __name__ == "__main__":
    main()
