"""Are the never-COTS ribbon reefs (Yonge, Carter, No Name) connected to the reefs that actually had
COTS? (lead's question, 2026-10-19, twentieth round.)

The model's internal COTS dispersal uses the Owen hydrodynamic matrices (the dataset's
reef_connectivity_<year>.csv equal the Owen global sub-matrix to 1e-16), not the distance-based
coral connectivity. Weights per sink i and source j (per ha of sink, per unit per-ha production at j),
averaged over the five hydrodynamic years:
    w[j, i] = C[j, i] * A_j / A_i * s_j        (C rows sources; A habitable area; s RME larval survival)

"Known COTS" source sets per outbreak wave (W3 1994-2001, W4 2009-2016, W5 2021-2027):
- predicted: outbreak probability P_j (original manta-based layer, gbrPredsAdj + 2026), wave mean;
- LTMP: AIMS reefs with manta CPUE >= 0.22 in the wave (matched to the nearest Owen reef <= 3 km);
- culled: reefs culled at >= 1 COTS/ha/yr in any year of the wave (cull data from 2002; W4, W5).
Pressure on sink i = sum_j w[j, i] x indicator_j (or P_j), split into domain / non-domain sources.

Outputs (licensed / internal inputs; git-ignored): sandbox/data/derived/ribbon_connectivity/

    python analyse_ribbon_connectivity.py [--output <dir>]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import geopandas as gpd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "sandbox" / "data"
OWEN = DATA / "OwenHydro" / "connectivityMatsYearlyGlobal"
DOMAIN = DATA / "Lizard_Historical_v0.4"
RME = DATA / "rme_ml_2025_06_05" / "data_files"
PREDS = [DATA / "gbrPredsAdj_20262408.csv", DATA / "gbrPreds2026Adj_20262408.csv"]
LTMP = DATA / "reef_cots.csv"
AIMS_REEFS = Path(r"C:\Users\smatthew\Documents\GitKraken\COTS-Coral-OutcomesV2\data\raw\aims-ltmp\2026-10-16\reefs.parquet")
CULL = DATA / "derived" / "cots_cull_removals_lizard_v0_3_v2" / "cull_reef_year.csv"
PROV = DOMAIN / "cots_connectivity" / "datasets" / "owen_global_2018_2023" / "provenance.json"
ALIASES = {"11-325": "10-441", "11-244e": "11-288", "11-244f": "11-303", "11-244g": "11-310",
           "11-244h": "11-311", "20198": "20-198"}
SINKS = {"Yonge": ["14-138"], "Carter": ["14-137"], "No Name": ["14-139"], "Martin": ["14-123"], "Linnet": ["14-126"],
         "Lizard": ["14-116b", "14-116c", "14-116d"], "MacGillivray": ["14-114"], "North Direction": ["14-143"],
         "Eyrie": ["14-118"]}
RIBBONS = ("Yonge", "Carter", "No Name")
WAVES = {"W3": (1994, 2001), "W4": (2009, 2016), "W5": (2021, 2027)}
DEFAULT_OUT = DATA / "derived" / "ribbon_connectivity"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def haversine_km(lat1, lon1, lat2, lon2):
    p1, p2 = np.radians(lat1), np.radians(lat2)
    a = np.sin((p2 - p1) / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(np.radians(lon2 - lon1) / 2) ** 2
    return 6371.0 * 2 * np.arcsin(np.sqrt(a))


def main(out: Path) -> None:
    if out.exists():
        raise FileExistsError(f"Refusing to overwrite {out}")
    shape = gpd.read_file(OWEN / "reefShapefiles" / "gbrShapeLL.shp")
    ids = [re.search(r"\(([^()]*)\)\s*$", str(n)).group(1).strip() for n in shape.reefName]
    names = [re.sub(r"\s*\([^()]*\)\s*$", "", str(n)) for n in shape.reefName]
    cent = shape.geometry.set_crs(4326).to_crs(3577).centroid.to_crs(4326)
    lat, lon = cent.y.to_numpy(), cent.x.to_numpy()
    index = {r: k for k, r in enumerate(ids)}
    domain_ids = set(pd.read_csv(DOMAIN / "site_to_reef.csv", dtype=str).reef_id)
    in_domain = np.array([r in domain_ids for r in ids])

    idl = pd.read_csv(RME / "id" / "id_list_2024_12_01.csv", comment="#", header=None, dtype={0: str})
    area_old = {str(r.iloc[0]): float(r.iloc[1]) * 100.0 * (1.0 - float(r.iloc[2])) for _, r in idl.iterrows()}
    surv = pd.read_csv(RME / "water_csv" / "COTS_LARVAL_REDUCTION_q3baseline.csv", dtype={"id": str}).set_index("id")
    old = [r if r in area_old else ALIASES.get(r) for r in ids]
    mapped = np.array([o is not None and o in area_old and o in surv.index for o in old])
    area = np.array([area_old[o] if m else 0.0 for o, m in zip(old, mapped)])
    surv_year = {int(k): int(v) for k, v in json.loads(PROV.read_text())["survival_year_by_hydrodynamic_year"].items()}

    sink_idx = {s: [index[r] for r in rs] for s, rs in SINKS.items()}
    mats = sorted(OWEN.glob("conMatCots*.csv"))
    W = {s: np.zeros(len(ids)) for s in SINKS}
    for path in mats:
        h = int(re.fullmatch(r"conMatCots(\d{4})-\d{2}\.csv", path.name).group(1))
        c = np.loadtxt(path, delimiter=",", dtype=np.float64)
        c[~np.isfinite(c).all(axis=1)] = 0.0
        s = np.array([float(surv.loc[o, str(surv_year[h])]) if m else 0.0 for o, m in zip(old, mapped)])
        for sink, cols in sink_idx.items():
            # Arrivals per ha of the sink (group) per unit per-ha production at each source.
            W[sink] += c[:, cols].sum(axis=1) * area * s / area[cols].sum() / len(mats)

    # Source sets.
    pr = pd.concat([pd.read_csv(p)[["year", "reefName", "outbrProb"]] for p in PREDS], ignore_index=True)
    pr["reef_id"] = pr.reefName.map(lambda n: re.search(r"\(([^()]*)\)\s*$", str(n)).group(1).strip())
    P = pr.pivot(index="reef_id", columns="year", values="outbrProb").reindex(ids)
    ltmp = pd.read_csv(LTMP)
    aims = pd.read_parquet(AIMS_REEFS)
    match = {}
    for r in aims.itertuples(index=False):
        d = haversine_km(r.latitude, r.longitude, lat, lon)
        k = int(np.argmin(d))
        if d[k] <= 3.0:
            match[r.aims_reef_name] = ids[k]
    cull = pd.read_csv(CULL, dtype={"reef_id": str})
    sets = {}
    for wave, (y0, y1) in WAVES.items():
        p = P.loc[:, [y for y in P.columns if y0 <= y <= y1]].mean(axis=1).fillna(0.0).to_numpy()
        lt = ltmp[ltmp.year.between(y0, y1)].groupby("reef_name").cotsptow.max()
        lt_ids = {match[n] for n, v in lt.items() if v >= 0.22 and n in match}
        cu = cull[cull.calendar_year.between(y0, y1) & (cull.culled_ha >= 1.0)].reef_id
        sets[wave] = {"predicted": p, "ltmp": np.array([r in lt_ids for r in ids], float),
                      "culled": np.array([r in set(cu) for r in ids], float),
                      "ltmp_reefs": sorted(lt_ids), "culled_reefs": sorted(set(cu))}

    rows = []
    for sink in SINKS:
        w = W[sink]
        row = {"sink": sink, "total_weight": float(w.sum()), "domain_share_of_weight": float(w[in_domain].sum() / w.sum())}
        for wave in WAVES:
            for key in ("predicted", "ltmp", "culled"):
                v = w * sets[wave][key]
                row[f"{wave}_{key}"] = float(v.sum())
                row[f"{wave}_{key}_domain_share"] = float(v[in_domain].sum() / v.sum()) if v.sum() > 0 else np.nan
        rows.append(row)
    table = pd.DataFrame(rows).set_index("sink")
    rel = table[[c for c in table if c.endswith(("_predicted", "_ltmp", "_culled"))]].div(table.loc["Martin", [c for c in table if c.endswith(("_predicted", "_ltmp", "_culled"))]])
    rel.columns = [f"{c}_vs_Martin" for c in rel.columns]
    table = table.join(rel)

    top = []
    for sink in RIBBONS + ("Martin", "Lizard"):
        w = W[sink]
        for rank, j in enumerate(np.argsort(-w)[:10], start=1):
            top.append({"sink": sink, "rank": rank, "source_id": ids[j], "source_name": names[j], "in_domain": bool(in_domain[j]),
                        "share_of_weight": float(w[j] / w.sum()), "lon": float(lon[j]), "lat": float(lat[j]),
                        **{f"P_{wv}": float(sets[wv]["predicted"][j]) for wv in WAVES},
                        **{f"ltmp_outbreak_{wv}": bool(sets[wv]["ltmp"][j]) for wv in WAVES},
                        **{f"culled_{wv}": bool(sets[wv]["culled"][j]) for wv in ("W4", "W5")}})
    top = pd.DataFrame(top)

    out.mkdir(parents=True)
    table.to_csv(out / "ribbon_pressure.csv")
    top.to_csv(out / "ribbon_top_sources.csv", index=False)
    fig, axes = plt.subplots(1, 3, figsize=(17, 6.2), sharex=True, sharey=True)
    p3 = sets["W3"]["predicted"]
    for ax, sink in zip(axes, ("Yonge", "No Name", "Martin")):
        w = W[sink] * p3
        sel = w > 0
        ax.scatter(lon, lat, s=3, color="#d6d5d0", linewidth=0)
        ax.scatter(lon[in_domain], lat[in_domain], s=6, color="#86b6ef", linewidth=0, label="domain reefs")
        sc = ax.scatter(lon[sel], lat[sel], s=8 + 160 * w[sel] / w.max(), c=w[sel], cmap="Oranges", linewidth=0.3, edgecolor="#52514e")
        for r in RIBBONS:
            j = sink_idx[r][0]
            ax.scatter(lon[j], lat[j], marker="^", s=70, color="#0b0b0b")
        j = sink_idx[sink][0]
        ax.scatter(lon[j], lat[j], marker="*", s=160, color="#0d366b", label=f"sink: {sink}")
        ax.set_xlim(144.2, 146.6)
        ax.set_ylim(-16.4, -13.4)
        ax.set_title(f"{sink}: outbreak-weighted larval sources, 1994-2001 (predicted P)", loc="left", fontsize=9)
        ax.grid(color="#e4e3df", linewidth=0.6)
        plt.colorbar(sc, ax=ax, fraction=0.04, pad=0.02)
    axes[0].legend(frameon=False, fontsize=8, loc="lower left")
    fig.tight_layout()
    fig.savefig(out / "ribbon_sources_map.png", dpi=150, facecolor="#fcfcfb")
    plt.close(fig)
    meta = {"status": "diagnostic", "weights": "w[j,i] = C[j,i] A_j / A_i s_j, mean over hydrodynamic years 2018-22",
            "waves": WAVES, "ltmp_matched_reefs": len(match),
            "source_sets": {wv: {"ltmp_reefs": sets[wv]["ltmp_reefs"], "culled_reefs": sets[wv]["culled_reefs"]} for wv in WAVES},
            "inputs": {p.name: sha256(p) for p in [*mats, *PREDS, LTMP, AIMS_REEFS, CULL, DOMAIN / "site_to_reef.csv"]},
            "script_sha256": sha256(Path(__file__)),
            "outputs": {n: sha256(out / n) for n in ("ribbon_pressure.csv", "ribbon_top_sources.csv", "ribbon_sources_map.png")}}
    (out / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
    pd.set_option("display.width", 250)
    cols = ["total_weight", "domain_share_of_weight"] + [f"{w}_{k}_vs_Martin" for w in WAVES for k in ("predicted", "ltmp", "culled")]
    print(table[cols].round(3).to_string())
    print(top[top.sink.isin(RIBBONS)][["sink", "rank", "source_name", "source_id", "in_domain", "share_of_weight", "lon",
                                       "P_W3", "P_W4", "P_W5", "ltmp_outbreak_W3", "ltmp_outbreak_W4"]].round(3).to_string(index=False))
    print({wv: (sets[wv]["ltmp_reefs"], len(sets[wv]["culled_reefs"])) for wv in WAVES})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    main(parser.parse_args().output.resolve())
