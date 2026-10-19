"""Which reefs supply the external COTS larvae reaching the calibration reefs? (lead's question,
2026-10-15, fifteenth round).

Re-computes, source reef by source reef, the external outbreak boundary that
build_lizard_cots_connectivity_owen.py sums for each domain reef:

    B[i, y, h] = sum_j  C_h[j, i] * A_j / A_i * s_j * P_j(y),   j outside the 113 domain reefs

(C_h: Owen connectivity for hydrodynamic year h, rows sources; A: RME habitable area; s: RME
larval survival; P: predicted probability that manta CPUE > 0.22). ADRIA multiplies B by the
outbreak production constant x Ext. Contributions are averaged over the five hydrodynamic years
(the expectation over the seeded random year draw) and checked against the dataset's
external_outbreak_boundary.csv. Calibration reefs are aggregated as sink-area-weighted means.

Outputs (licensed / internal inputs, so git-ignored): sandbox/data/derived/external_supply_sources/

    python analyse_external_sources.py [--output <dir>]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import geopandas as gpd
from shapely.geometry import MultiPoint, Point
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "sandbox" / "data"
OWEN = DATA / "OwenHydro" / "connectivityMatsYearlyGlobal"
SHAPE = OWEN / "reefShapefiles" / "gbrShapeLL.shp"
DOMAIN = DATA / "Lizard_Historical_v0.4"
DATASETS = DOMAIN / "cots_connectivity" / "datasets"
LAYERS = {"original": "owen_global_2018_2023", "corrected": "owen_global_2018_2023_corrected"}
RME = DATA / "rme_ml_2025_06_05" / "data_files"
CAL = {"Lizard Island Reef": ["14-116b", "14-116c", "14-116d"], "MacGillivray Reef": ["14-114"],
       "North Direction Reef": ["14-143"], "Eyrie Reef": ["14-118"]}
YEARS = (1995, 2005, 2012, 2018, 2020, 2022, 2024, 2025, 2026)
ALIASES = {"11-325": "10-441", "11-244e": "11-288", "11-244f": "11-303", "11-244g": "11-310",
           "11-244h": "11-311", "20198": "20-198"}
DEFAULT_OUT = DATA / "derived" / "external_supply_sources"
TOP = 10


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def reef_id(name: str) -> str:
    return re.search(r"\(([^()]*)\)\s*$", str(name)).group(1).strip()


def haversine_km(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    a = np.sin((p2 - p1) / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(np.radians(lon2 - lon1) / 2) ** 2
    return 6371.0 * 2 * np.arcsin(np.sqrt(a))


def main(out: Path) -> None:
    if out.exists():
        raise FileExistsError(f"Refusing to overwrite {out}")
    shape = gpd.read_file(SHAPE)
    ids = [reef_id(n) for n in shape.reefName]
    names = {i: re.sub(r"\s*\([^()]*\)\s*$", "", str(n)) for i, n in zip(ids, shape.reefName)}
    cent = shape.geometry.set_crs(4326).to_crs(3577).centroid.to_crs(4326)  # gbrShapeLL: lon/lat, no .prj CRS
    lat, lon = dict(zip(ids, cent.y)), dict(zip(ids, cent.x))
    index = {r: k for k, r in enumerate(ids)}
    sites = pd.read_csv(DOMAIN / "site_to_reef.csv", dtype=str)
    local = sorted({index[r] for r in sites.reef_id})
    if len(local) != 113:
        raise ValueError("Expected 113 domain reefs")
    local_ids = [ids[k] for k in local]
    lpos = {r: k for k, r in enumerate(local_ids)}

    idl = pd.read_csv(RME / "id" / "id_list_2024_12_01.csv", comment="#", header=None, dtype={0: str})
    area_old = {str(r.iloc[0]): float(r.iloc[1]) * 100.0 * (1.0 - float(r.iloc[2])) for _, r in idl.iterrows()}
    surv = pd.read_csv(RME / "water_csv" / "COTS_LARVAL_REDUCTION_q3baseline.csv", dtype={"id": str}).set_index("id")
    old = [r if r in area_old else ALIASES.get(r) for r in ids]
    mapped = np.array([o is not None and o in area_old and o in surv.index for o in old])
    external = mapped.copy()
    external[local] = False
    area = np.array([area_old[o] if m else 0.0 for o, m in zip(old, mapped)])

    prov = json.loads((DATASETS / LAYERS["original"] / "provenance.json").read_text())
    surv_year = {int(k): int(v) for k, v in prov["survival_year_by_hydrodynamic_year"].items()}
    mats = sorted(OWEN.glob("conMatCots*.csv"))
    hydro = [int(re.fullmatch(r"conMatCots(\d{4})-\d{2}\.csv", p.name).group(1)) for p in mats]

    preds = {}
    for layer, ds in LAYERS.items():
        p = json.loads((DATASETS / ds / "provenance.json").read_text())
        files = [ROOT / s["path"] for s in p["sources"] if s["path"].endswith(".csv") and ("Pred" in s["path"])]
        d = pd.concat([pd.read_csv(f)[["year", "reefName", "outbrProb"]] for f in files], ignore_index=True)
        d["reef_id"] = d.reefName.map(reef_id)
        preds[layer] = d.pivot(index="reef_id", columns="year", values="outbrProb").loc[ids, list(YEARS)].to_numpy()

    sinks = {cal: [lpos[r] for r in rs] for cal, rs in CAL.items()}
    sink_cols = {cal: [local[k] for k in ks] for cal, ks in sinks.items()}
    # Per-source weight for each calibration reef: mean over hydro years of
    # sum_i C[j, i] A_j s_j / sum_i A_i (sink-area-weighted mean of the per-ha boundary).
    weight = {cal: np.zeros(len(ids)) for cal in CAL}
    boundary_check = []
    for path, h in zip(mats, hydro):
        m = np.loadtxt(path, delimiter=",", dtype=np.float64)
        m[~np.isfinite(m).all(axis=1)] = 0.0
        s = np.where(mapped, [float(surv.loc[o, str(surv_year[h])]) if mm else 0.0 for o, mm in zip(old, mapped)], 0.0)
        for cal, cols in sink_cols.items():
            w = (m[:, cols].sum(axis=1) * area * s) / area[cols].sum()
            w[~external] = 0.0
            weight[cal] += w / len(mats)
            for layer in LAYERS:
                for yk, y in enumerate(YEARS):
                    boundary_check.append({"layer": layer, "hydro": h, "cal": cal, "year": y,
                                           "recomputed": float((w * preds[layer][:, yk]).sum()),
                                           "area_weighted_by": cols})
    # Check against the files ADRIA reads.
    worst = 0.0
    for layer, ds in LAYERS.items():
        b = pd.read_csv(DATASETS / ds / "annual_forcing" / "external_outbreak_boundary.csv", dtype={"reef_id": str})
        for row in boundary_check:
            if row["layer"] != layer:
                continue
            cols = row["area_weighted_by"]
            rs = [ids[c] for c in cols]
            sub = b[(b.hydrodynamic_year == row["hydro"]) & (b.calendar_year == row["year"]) & b.reef_id.isin(rs)].set_index("reef_id")
            ref = float(sum(sub.loc[r, "survival_weighted_probability_area_connectivity"] * area[index[r]] for r in rs) / area[cols].sum())
            worst = max(worst, abs(ref - row["recomputed"]) / max(abs(ref), 1e-12))
    if worst > 1e-9:
        raise ValueError(f"Recomputed boundary differs from the dataset file (max rel {worst:.2e})")

    lizard = (np.mean([lat[r] for r in CAL["Lizard Island Reef"]]), np.mean([lon[r] for r in CAL["Lizard Island Reef"]]))
    dist_liz = np.array([haversine_km(lizard[0], lizard[1], lat[r], lon[r]) for r in ids])
    domain_radius = float(dist_liz[local].max())
    # Domain footprint: convex hull of the 113 domain-reef centroids. External reefs inside it are
    # unmodelled reefs interspersed among the modelled ones; outside it, truly beyond the domain.
    hull = MultiPoint([(lon[ids[k]], lat[ids[k]]) for k in local]).convex_hull
    in_hull = np.array([hull.contains(Point(lon[r], lat[r])) for r in ids]) & external
    rows, totals = [], []
    for cal in CAL:
        c_lat, c_lon = np.mean([lat[r] for r in CAL[cal]]), np.mean([lon[r] for r in CAL[cal]])
        dist = np.array([haversine_km(c_lat, c_lon, lat[r], lon[r]) for r in ids])
        for yk, y in enumerate(YEARS):
            for layer in LAYERS:
                contrib = weight[cal] * preds[layer][:, yk]
                total = contrib.sum()
                order = np.argsort(-contrib)
                share = contrib / total if total > 0 else contrib * 0
                within = lambda km: float(share[dist <= km].sum())  # noqa: E731
                totals.append({"calibration_reef": cal, "year": y, "layer": layer, "boundary_per_ha": float(total),
                               "top1_share": float(share[order[0]]), "top5_share": float(share[order[:5]].sum()),
                               "top10_share": float(share[order[:10]].sum()), "n_sources_for_90pct": int(np.searchsorted(np.cumsum(share[order]), 0.9) + 1),
                               "share_within_50km": within(50), "share_within_100km": within(100),
                               "share_within_200km": within(200), "share_beyond_500km": float(share[dist > 500].sum()),
                               "share_external_inside_domain_footprint": float(share[in_hull].sum())})
                for rank, j in enumerate(order[:TOP], start=1):
                    rows.append({"calibration_reef": cal, "year": y, "layer": layer, "rank": rank, "source_reef_id": ids[j],
                                 "source_name": names[ids[j]], "distance_km": float(dist[j]),
                                 "source_lat": float(lat[ids[j]]), "source_lon": float(lon[ids[j]]),
                                 "p_outbreak_original": float(preds["original"][j, yk]),
                                 "p_outbreak_corrected": float(preds["corrected"][j, yk]),
                                 "connection_weight": float(weight[cal][j]), "share": float(share[j])})
    out.mkdir(parents=True)
    top, tot = pd.DataFrame(rows), pd.DataFrame(totals)
    top.to_csv(out / "top_sources.csv", index=False)
    tot.to_csv(out / "supply_concentration.csv", index=False)
    # Source weights (connection strength regardless of outbreak state), all calibration reefs.
    wt = pd.DataFrame({"source_reef_id": ids, "source_name": [names[r] for r in ids], "lat": [lat[r] for r in ids],
                       "lon": [lon[r] for r in ids], "distance_from_lizard_km": dist_liz, "in_domain": ~external & mapped,
                       "external_inside_domain_footprint": in_hull,
                       **{f"weight_{CAL_SHORT}": weight[cal] for cal, CAL_SHORT in zip(CAL, ("liz", "macg", "nd", "eyr"))}})
    wt.to_csv(out / "source_weights.csv", index=False)

    # Map: connection weight to Lizard (all years) and 2025 contribution under each layer.
    fig, axes = plt.subplots(1, 3, figsize=(17, 6.4), sharex=True, sharey=True)
    yk = YEARS.index(2025)
    panels = (("A  Connection to Lizard (outbreak-free weight)", weight["Lizard Island Reef"]),
              ("B  Lizard supply 2025, original layer", weight["Lizard Island Reef"] * preds["original"][:, yk]),
              ("C  Lizard supply 2025, corrected layer", weight["Lizard Island Reef"] * preds["corrected"][:, yk]))
    lats, lons = np.array([lat[r] for r in ids]), np.array([lon[r] for r in ids])
    vmax = max(p[1].max() for p in panels[1:])
    for ax, (title, v) in zip(axes, panels):
        ax.scatter(lons, lats, s=3, color="#d6d5d0", linewidth=0)
        ax.scatter(lons[local], lats[local], s=6, color="#86b6ef", linewidth=0, label="113 domain reefs")
        sel = v > 0
        sc = ax.scatter(lons[sel], lats[sel], s=8 + 140 * v[sel] / v.max(), c=v[sel], cmap="Oranges",
                        vmin=0, vmax=v.max() if title.startswith("A") else vmax, linewidth=0.3, edgecolor="#52514e")
        ax.scatter([lizard[1]], [lizard[0]], marker="*", s=120, color="#0d366b", label="Lizard")
        ax.set_title(title, loc="left", fontsize=10)
        ax.set_xlim(143.0, 147.5)
        ax.set_ylim(-17.5, -12.0)
        ax.grid(color="#e4e3df", linewidth=0.6)
        ax.set_xlabel("Longitude")
        plt.colorbar(sc, ax=ax, fraction=0.04, pad=0.02)
    axes[0].set_ylabel("Latitude")
    axes[0].legend(frameon=False, fontsize=8, loc="lower left")
    fig.tight_layout()
    fig.savefig(out / "external_sources_lizard.png", dpi=150, facecolor="#fcfcfb")
    plt.close(fig)

    meta = {"status": "diagnostic", "boundary_check_max_rel_delta": worst, "domain_radius_km_from_lizard": domain_radius,
            "external_reefs_inside_domain_footprint": int(in_hull.sum()),
            "years": list(YEARS), "hydrodynamic_years": hydro, "calibration_reef_ids": CAL,
            "inputs": {str(p.relative_to(ROOT)): sha256(p) for p in [SHAPE, DOMAIN / "site_to_reef.csv", *mats]
                       + [DATASETS / ds / "annual_forcing" / "external_outbreak_boundary.csv" for ds in LAYERS.values()]},
            "script_sha256": sha256(Path(__file__)),
            "outputs": {n: sha256(out / n) for n in ("top_sources.csv", "supply_concentration.csv", "source_weights.csv",
                                                     "external_sources_lizard.png")}}
    (out / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
    pd.set_option("display.width", 250)
    print(f"boundary check max rel delta {worst:.2e}; domain radius from Lizard {domain_radius:.0f} km")
    print(tot.round(3).to_string(index=False))
    show = top[(top["rank"] <= 6) & top.year.isin([2012, 2025])]
    print(show[["calibration_reef", "year", "layer", "rank", "source_reef_id", "source_name", "distance_km",
                "p_outbreak_original", "p_outbreak_corrected", "share"]].round(3).to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    main(parser.parse_args().output.resolve())
