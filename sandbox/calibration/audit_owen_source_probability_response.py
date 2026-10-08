"""Reconstruct Owen upstream sources and screen probability-response leverage.

This is an offline boundary audit, not a model run or proposed biological
estimate. It leaves the frozen COTS score, input data, and core equations
unchanged. The p^2 and p^4 responses are labelled counterfactuals and use one
global normalizer to match the linear response's four-reef 2012-15 mean
boundary coefficient. No per-reef or trough-window fitting is performed.
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
TARGETS = ("Lizard Island Reef", "MacGillivray Reef", "North Direction Reef", "Eyrie Reef")
WINDOWS = {"first_wave": range(1994, 1998), "trough": range(2004, 2010),
           "late_wave_supply": range(2012, 2016)}
RESPONSES = {"linear": 1, "p_squared": 2, "p_fourth": 4}
PRODUCTION = 142.66240759512007  # recruits/ha/year, archived joint4 value


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
        "site_map": DOMAIN / "site_to_reef.csv",
        "boundary": FORCING / "external_outbreak_boundary.csv",
        "reference_trajectories": REFERENCE / "trajectories.csv",
        "reference_flows": REFERENCE / "flows.csv",
        "reference_metadata": REFERENCE / "metadata.toml",
    }
    matrix_paths = sorted(builder.SOURCE_ROOT.glob("conMatCots*.csv"),
                          key=builder.matrix_year)
    for path in matrix_paths:
        sources[path.name] = path

    reef_shape = gpd.read_file(sources["shape_shp"])
    global_ids = [builder.reef_id_from_name(name) for name in reef_shape.reefName]
    if len(global_ids) != 3861 or len(set(global_ids)) != len(global_ids):
        raise ValueError("Unexpected Owen global reef IDs")
    global_lookup = {reef_id: index for index, reef_id in enumerate(global_ids)}
    probabilities, years, _ = builder.load_prediction_matrix(global_ids)
    year_lookup = {year: index for index, year in enumerate(years)}

    id_list = pd.read_csv(sources["rme_id_list"], comment="#", header=None,
                          dtype={0: str})
    old_area_ha = {str(row.iloc[0]): float(row.iloc[1]) * 100.0 *
                   (1.0 - float(row.iloc[2])) for _, row in id_list.iterrows()}
    rme_survival = pd.read_csv(builder.DEFAULT_SURVIVAL, usecols=["id"], dtype=str)
    sources["rme_survival"] = builder.DEFAULT_SURVIVAL
    survival_ids = set(rme_survival.id.astype(str))
    mapped_old_ids = [reef_id if reef_id in old_area_ha else
                      builder.NEW_TO_OLD_ALIASES.get(reef_id) for reef_id in global_ids]
    mapped = np.array([old_id in old_area_ha and old_id in survival_ids
                       for old_id in mapped_old_ids])
    source_area = np.array([old_area_ha[old_id] if okay else 0.0
                            for old_id, okay in zip(mapped_old_ids, mapped)])

    sites = pd.read_csv(sources["site_map"])
    sites["target_reef"] = sites.reef_name.str.split(" (", regex=False).str[0]
    local_ids = set(sites.reef_id.astype(str))
    external = mapped & ~np.isin(global_ids, list(local_ids))
    if int(external.sum()) != 3692:
        raise ValueError("External source count changed")
    sink_specs: dict[str, list[tuple[int, float, float]]] = {}
    for target in TARGETS:
        rows = sites[sites.target_reef == target]
        if rows.empty:
            raise ValueError(f"Target reef missing: {target}")
        shares = rows.groupby("reef_id", sort=True).area_m2.sum()
        shares = shares / shares.sum()
        sink_specs[target] = [(global_lookup[reef_id], float(weight),
                               old_area_ha[mapped_old_ids[global_lookup[reef_id]]])
                              for reef_id, weight in shares.items()]

    # Source coefficients before multiplying annual outbreak probability.
    factors: dict[tuple[str, int], np.ndarray] = {}
    for path in matrix_paths:
        hydro = builder.matrix_year(path)
        matrix = np.loadtxt(path, delimiter=",", dtype=np.float64)
        if matrix.shape != (len(global_ids), len(global_ids)):
            raise ValueError(f"Bad matrix shape: {path}")
        nonfinite = ~np.isfinite(matrix)
        full_bad_rows = nonfinite.all(axis=1)
        if np.any(nonfinite & ~full_bad_rows[:, None]):
            raise ValueError(f"Partly non-finite source row: {path}")
        matrix[full_bad_rows, :] = 0.0
        for target, sinks in sink_specs.items():
            factor = np.zeros(len(global_ids), dtype=np.float64)
            for sink_index, weight, sink_area in sinks:
                factor += weight * matrix[:, sink_index] * source_area / sink_area
            factor[~external] = 0.0
            factors[target, hydro] = factor
        del matrix

    boundary = pd.read_csv(sources["boundary"])
    max_reconstruction_error = 0.0
    for target, sinks in sink_specs.items():
        for hydro in sorted({key[1] for key in factors}):
            recorded = np.zeros(len(years), dtype=np.float64)
            for sink_index, weight, _ in sinks:
                sink_id = global_ids[sink_index]
                subset = boundary[(boundary.reef_id == sink_id) &
                                  (boundary.hydrodynamic_year == hydro)].sort_values(
                                      "calendar_year")
                if subset.calendar_year.tolist() != years:
                    raise ValueError("Boundary years do not align")
                recorded += weight * subset.probability_area_connectivity.to_numpy()
            reconstructed = factors[target, hydro] @ probabilities
            max_reconstruction_error = max(max_reconstruction_error,
                                           float(np.max(np.abs(recorded-reconstructed))))
    if max_reconstruction_error > 1e-10:
        raise ValueError(f"Boundary reconstruction failed: {max_reconstruction_error}")

    trajectories = pd.read_csv(sources["reference_trajectories"])
    trajectories = trajectories[trajectories.treatment == "joint4_reference"]
    flows = pd.read_csv(sources["reference_flows"])
    flows = flows[flows.treatment == "joint4_reference"]
    selected = trajectories[trajectories.year.isin(years)]
    if len(selected) != len(TARGETS) * 34:  # 1991-2024
        raise ValueError("Incomplete selected-hydrodynamic-year trajectory")

    raw_rows = []
    contributions: dict[str, list[np.ndarray]] = {target: [] for target in TARGETS}
    for row in selected.itertuples():
        target, year, hydro = row.reef_name, int(row.year), int(row.forcing_year)
        if target not in TARGETS:
            raise ValueError("Unexpected target reef")
        p = probabilities[:, year_lookup[year]]
        factor = factors[target, hydro]
        coefficients = {name: float(factor @ (p ** power))
                        for name, power in RESPONSES.items()}
        if year in WINDOWS["trough"]:
            contributions[target].append(factor * p)
        logged = flows[(flows.reef_name == target) & (flows.year == year)]
        if len(logged) != 1:
            raise ValueError("Missing logged reef-year flux")
        modeled = coefficients["linear"] * PRODUCTION
        if abs(modeled - float(logged.external_pelagic.iloc[0])) > 1e-9:
            raise ValueError(f"Prediction boundary and logged flux differ: {target}, {year}")
        raw_rows.append({"reef_name": target, "calendar_year": year,
                         "hydrodynamic_year": hydro,
                         **{name: value for name, value in coefficients.items()}})
    raw = pd.DataFrame(raw_rows)
    late = raw[raw.calendar_year.isin(WINDOWS["late_wave_supply"])]
    normalizers = {"linear": 1.0}
    for name in RESPONSES:
        if name != "linear":
            normalizers[name] = float(late.linear.mean() / late[name].mean())
    for name in RESPONSES:
        raw[name + "_normalized_flux_recruits_ha_year"] = (
            raw[name] * normalizers[name] * PRODUCTION
        )
    raw.to_csv(out / "boundary_response_by_reef_year.csv", index=False)

    summary_rows = []
    for target in TARGETS:
        reef_rows = raw[raw.reef_name == target]
        for name in RESPONSES:
            for window, window_years in WINDOWS.items():
                values = reef_rows[reef_rows.calendar_year.isin(window_years)][
                    name + "_normalized_flux_recruits_ha_year"]
                summary_rows.append({"reef_name": target, "response": name,
                                     "window": window, "mean_flux_recruits_ha_year":
                                     float(values.mean()), "minimum_flux_recruits_ha_year":
                                     float(values.min()), "global_late_normalizer":
                                     normalizers[name]})
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(out / "response_window_summary.csv", index=False)

    source_rows = []
    for target in TARGETS:
        mean_contribution = np.mean(contributions[target], axis=0)
        total = float(mean_contribution.sum())
        order = np.argsort(mean_contribution)[::-1]
        for rank, index in enumerate(order[:20], start=1):
            source_rows.append({
                "sink_reef": target, "rank": rank,
                "source_reef_id": global_ids[index],
                "source_reef_name": str(reef_shape.reefName.iloc[index]),
                "mean_trough_probability": float(np.mean(probabilities[
                    index, [year_lookup[year] for year in WINDOWS["trough"]]])),
                "mean_late_probability": float(np.mean(probabilities[
                    index, [year_lookup[year] for year in WINDOWS["late_wave_supply"]]])),
                "mean_trough_boundary_coefficient": float(mean_contribution[index]),
                "share_of_trough_boundary": float(mean_contribution[index] / total),
                "cumulative_share": float(mean_contribution[order[:rank]].sum() / total),
            })
    top = pd.DataFrame(source_rows)
    top.to_csv(out / "top_trough_sources.csv", index=False)

    metadata = {
        "status": "offline_counterfactual_not_promoted",
        "assumption": "Owen connectivity includes larval mortality; use potential coefficient",
        "source_count": int(external.sum()),
        "source_to_sink_orientation": "global matrix rows are sources, columns are sinks",
        "coefficient_units": "dimensionless; times production gives pre-settlement recruits/ha/year",
        "production_recruits_ha_year": PRODUCTION,
        "window_years": {name: list(years) for name, years in WINDOWS.items()},
        "response_functions": {name: f"p^{power}" for name, power in RESPONSES.items()},
        "normalization": "one global scalar per alternative preserves four-reef 2012-15 mean boundary coefficient",
        "normalizers": normalizers,
        "max_boundary_reconstruction_error": max_reconstruction_error,
        "inputs": {name: {"path": str(path), "sha256": digest(path)}
                   for name, path in sources.items()},
        "script_sha256": digest(Path(__file__)),
        "builder_sha256": digest(Path(builder.__file__)),
        "outputs": {path.name: digest(path) for path in
                    (out / "boundary_response_by_reef_year.csv",
                     out / "response_window_summary.csv",
                     out / "top_trough_sources.csv")},
    }
    (out / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n",
                                       encoding="utf-8")
    print("Reconstruction max error:", max_reconstruction_error)
    print("Global late-wave normalizers:", normalizers)
    print(summary[summary.window == "trough"].pivot(
        index="reef_name", columns="response", values="mean_flux_recruits_ha_year"
    ).round(3).to_string())
    print("Top-five source shares:")
    print(top[top["rank"] == 5][["sink_reef", "cumulative_share"]].to_string(index=False))
    print("Saved", out)


if __name__ == "__main__":
    main()
