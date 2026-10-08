"""Build an opt-in, source-specific outbreak-history condition proxy boundary.

The proxy is NOT measured coral or maternal health. The linear Owen boundary
must replay before writing a candidate. No domain forcing is overwritten.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "sandbox/domain_building"))
import build_lizard_cots_connectivity_owen as builder  # noqa: E402

DOMAIN = ROOT / "sandbox/data/Lizard_Historical_v0.2"
FORCING = DOMAIN / "cots_connectivity/datasets/owen_global_2018_2023/annual_forcing"
REFERENCE = ROOT / "sandbox/calibration/runs/20261002T120000_owen_recruitment_blackout"
BEST = ROOT / "sandbox/calibration/runs/expanded_cotsconn_pilot_seed20260930/best_summary.csv"
TARGETS = ("Lizard Island Reef", "MacGillivray Reef", "North Direction Reef", "Eyrie Reef")
WINDOWS = {"first_wave": range(1994, 1998), "trough": range(2004, 2010),
           "late_wave": range(2012, 2016)}


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
    out = ROOT / "sandbox/calibration/runs" / args.run_id
    out.mkdir(exist_ok=False)
    sources = {
        "prediction": builder.PREDICTIONS,
        "shape_shp": builder.REEF_SHAPE,
        "shape_dbf": builder.REEF_SHAPE.with_suffix(".dbf"),
        "rme_id_list": builder.RME_ID_LIST,
        "rme_survival": builder.DEFAULT_SURVIVAL,
        "site_map": DOMAIN / "site_to_reef.csv",
        "archived_boundary": FORCING / "external_outbreak_boundary.csv",
        "reference_trajectories": REFERENCE / "trajectories.csv",
        "candidate": BEST,
        "builder": Path(builder.__file__),
    }
    reef_shape = gpd.read_file(sources["shape_shp"])
    global_ids = [builder.reef_id_from_name(name) for name in reef_shape.reefName]
    if len(global_ids) != 3861 or len(set(global_ids)) != len(global_ids):
        raise ValueError("Unexpected global reef IDs")
    global_lookup = {reef_id: index for index, reef_id in enumerate(global_ids)}
    p, prediction_years, _ = builder.load_prediction_matrix(global_ids)
    if prediction_years != list(range(1991, 2026)):
        raise ValueError("Unexpected prediction years")
    tau = float(pd.read_csv(BEST).sort_values("loss").iloc[0].tau_condition)
    alpha = 1.0 / tau
    q = np.zeros_like(p)
    q[:, 0] = 0.8
    for t in range(1, len(prediction_years)):
        q[:, t] = (1.0 - alpha) * q[:, t-1] + alpha * (1.0 - p[:, t-1])
    if not np.all((q >= 0.0) & (q <= 1.0)):
        raise ValueError("Condition proxy out of bounds")
    proxy_weight = p * q**2

    sites = pd.read_csv(sources["site_map"])
    local_ids = set(sites.reef_id.astype(str))
    id_list = pd.read_csv(sources["rme_id_list"], comment="#", header=None,
                          dtype={0: str})
    old_area = {str(row.iloc[0]): float(row.iloc[1]) * 100.0 *
                (1.0 - float(row.iloc[2])) for _, row in id_list.iterrows()}
    survival_ids = set(pd.read_csv(sources["rme_survival"], usecols=["id"],
                                   dtype=str).id.astype(str))
    mapped_ids = [reef_id if reef_id in old_area else
                  builder.NEW_TO_OLD_ALIASES.get(reef_id) for reef_id in global_ids]
    mapped = np.array([old_id in old_area and old_id in survival_ids
                       for old_id in mapped_ids])
    external = mapped & ~np.isin(global_ids, list(local_ids))
    if int(external.sum()) != 3692:
        raise ValueError("External source count changed")
    source_area = np.array([old_area[old_id] if okay else 0.0
                            for old_id, okay in zip(mapped_ids, mapped)])
    local_order = sorted(local_ids, key=global_lookup.__getitem__)
    local_indices = np.array([global_lookup[reef_id] for reef_id in local_order])
    sink_area = np.array([old_area[mapped_ids[index]] for index in local_indices])
    rows = []
    matrix_paths = sorted(builder.SOURCE_ROOT.glob("conMatCots*.csv"),
                          key=builder.matrix_year)
    for path in matrix_paths:
        sources[path.name] = path
        hydro = builder.matrix_year(path)
        matrix = np.loadtxt(path, delimiter=",", dtype=np.float64)
        if matrix.shape != (len(global_ids), len(global_ids)):
            raise ValueError(f"Bad matrix shape: {path}")
        nonfinite = ~np.isfinite(matrix)
        full_bad_rows = nonfinite.all(axis=1)
        if np.any(nonfinite & ~full_bad_rows[:, None]):
            raise ValueError(f"Partly non-finite source row: {path}")
        matrix[full_bad_rows, :] = 0.0
        pressure = matrix[np.ix_(external, local_indices)] * (
            source_area[external, None] / sink_area[None, :])
        linear = pressure.T @ p[external, :]
        proxy = pressure.T @ proxy_weight[external, :]
        for sink, reef_id in enumerate(local_order):
            for year_idx, year in enumerate(prediction_years):
                rows.append((reef_id, year, hydro, linear[sink, year_idx],
                             proxy[sink, year_idx]))
        del matrix
    result = pd.DataFrame(rows, columns=["reef_id", "calendar_year",
                                         "hydrodynamic_year", "linear", "proxy"])
    archived = pd.read_csv(sources["archived_boundary"])
    merged = result.merge(archived, on=["reef_id", "calendar_year",
                                        "hydrodynamic_year"], validate="one_to_one")
    if len(merged) != len(result) or len(archived) != len(result):
        raise ValueError("Boundary rows do not match")
    replay_error = float(np.max(np.abs(
        merged.linear - merged.probability_area_connectivity)))
    if replay_error > 1e-10:
        raise ValueError(f"Archived boundary replay failed: {replay_error}")

    reference = pd.read_csv(sources["reference_trajectories"])
    reference = reference[reference.treatment == "joint4_reference"]
    reference = reference[reference.year.isin(prediction_years)]
    if len(reference) != len(TARGETS) * 34:  # Model ends in 2024; predictions include 2025.
        raise ValueError("Incomplete reference reef-years")
    sites["target_reef"] = sites.reef_name.str.split(" (", regex=False).str[0]
    target_weights = {}
    for target in TARGETS:
        subset = sites[sites.target_reef == target]
        shares = subset.groupby("reef_id").area_m2.sum()
        target_weights[target] = (shares / shares.sum()).to_dict()
    weighted_rows = []
    lookup = result.set_index(["reef_id", "calendar_year", "hydrodynamic_year"])
    for row in reference.itertuples():
        values = {name: sum(weight * float(lookup.loc[
            reef_id, int(row.year), int(row.forcing_year)][name])
            for reef_id, weight in target_weights[row.reef_name].items())
            for name in ("linear", "proxy")}
        weighted_rows.append({"reef_name": row.reef_name, "year": int(row.year),
                              **values})
    weighted = pd.DataFrame(weighted_rows)
    late = weighted[weighted.year.isin(WINDOWS["late_wave"])]
    normalizer = float(late.linear.mean() / late.proxy.mean())
    result["condition_proxy_coefficient"] = result.proxy * normalizer
    weighted["condition_proxy"] = weighted.proxy * normalizer
    result = result[["reef_id", "calendar_year", "hydrodynamic_year",
                     "linear", "condition_proxy_coefficient"]]
    result.to_csv(out / "external_source_condition_boundary.csv", index=False)
    summary = []
    for target in TARGETS:
        subset = weighted[weighted.reef_name == target]
        for window, years in WINDOWS.items():
            for response in ("linear", "condition_proxy"):
                summary.append({"reef_name": target, "window": window,
                                "response": response,
                                "mean_coefficient": float(subset[
                                    subset.year.isin(years)][response].mean())})
    pd.DataFrame(summary).to_csv(out / "response_window_summary.csv", index=False)
    metadata = {
        "status": "external_source_condition_proxy_counterfactual_not_promoted",
        "equation": "q[1991]=0.8; q[t]=(1-1/tau)*q[t-1]+(1/tau)*(1-p[t-1]); source weight=p[t]*q[t]^2",
        "tau_years": tau,
        "probability_definition": "P(upstream reef COTS > 0.22/tow)",
        "proxy_warning": "1-p is not measured coral or maternal condition; adult-density fertilisation at external sources is untested",
        "assumption": "Owen connectivity embeds pelagic mortality; no extra pelagic factor",
        "source_to_sink_orientation": "global matrix rows are sources, columns are sinks",
        "coefficient_units": "dimensionless; times boundary production gives pre-settlement recruits/ha/year",
        "global_late_wave_normalizer": normalizer,
        "max_linear_boundary_replay_error": replay_error,
        "external_source_count": int(external.sum()),
        "inputs": {name: {"path": str(path), "sha256": digest(path)}
                   for name, path in sources.items()},
        "script_sha256": digest(Path(__file__)),
        "outputs": {path.name: digest(path) for path in
                    (out / "external_source_condition_boundary.csv",
                     out / "response_window_summary.csv")},
    }
    (out / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n",
                                       encoding="utf-8")
    print("Archived boundary max replay error:", replay_error)
    print("Global late-wave normalizer:", normalizer)
    print(pd.DataFrame(summary).pivot_table(index=["reef_name", "window"],
        columns="response", values="mean_coefficient").round(4).to_string())
    print("Saved", out)


if __name__ == "__main__":
    main()
