"""What keeps Yonge, Carter and No Name reefs free of COTS? (lead's decisions, 2026-10-19,
twentieth round: most likely connectivity plus strong wave exposure on the eastern flank; check
in-water chl because the satellite series may obscure lower outer-shelf levels.)

Data checks only (no model run), comparing the never-COTS reefs with outbreak reefs:
1. Wave exposure per site (GBR ReefGuidance processed rasters, ~10 m): mean significant wave
   height Hs, peak period Tp and a wave-power index Hs^2 Tp, over each site polygon; per reef the
   area-weighted mean and the most exposed site; domain percentile ranks.
2. Chlorophyll: AIMS LTCMP in-water (1992-2009) median at the nearest station (>= 100 samples)
   against the satellite-derived in-water-equivalent p90 (chla_index_lizard_1993_2026_inwater, RMA).
3. Connectivity (Owen 2018-22, mean over hydrodynamic years): larval arrivals per ha at each reef if
   every reef produced equally per ha (sum_j C[j, i] A_j / A_i), split into own reef, other domain
   reefs and non-domain reefs.

Outputs (licensed / internal inputs; git-ignored): sandbox/data/derived/never_cots_diagnostic/

    python analyse_never_cots_reefs.py [--output <dir>]
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
import rasterio
from rasterio.mask import mask

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "sandbox" / "data"
SITES = DATA / "Lizard_Historical_v0.4" / "spatial" / "lizard_cluster.gpkg"
WAVES = Path(r"C:\Users\smatthew\Documents\GitKraken\COTS Fine Scale SDM\DataProcessing\data\GBR-ReefGuidance-processed_2025-11-19")
REGIONS = ("Cairns-Cooktown", "FarNorthern")
INWATER = DATA / "GBR_chlorophyll_monitoring_1992-2009.csv"
CHL_INDEX = DATA / "derived" / "chla_index_lizard_1993_2026_inwater" / "chla_index_reef_year.csv"
OWEN = DATA / "OwenHydro" / "connectivityMatsYearlyGlobal"
RME = DATA / "rme_ml_2025_06_05" / "data_files"
ALIASES = {"11-325": "10-441", "11-244e": "11-288", "11-244f": "11-303", "11-244g": "11-310",
           "11-244h": "11-311", "20198": "20-198"}
GROUPS = {"absence": ["14-138", "14-137", "14-139"],
          "outbreak": ["14-123", "14-126", "14-116b", "14-116c", "14-116d", "14-114", "14-143", "14-118"]}
LABEL = {"14-138": "Yonge", "14-137": "Carter", "14-139": "No Name", "14-123": "Martin", "14-126": "Linnet",
         "14-116b": "Lizard NE", "14-116c": "Lizard Coconut Bay", "14-116d": "Lizard Lagoon",
         "14-114": "MacGillivray", "14-143": "North Direction", "14-118": "Eyrie"}
DEFAULT_OUT = DATA / "derived" / "never_cots_diagnostic"


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


def site_waves(sites: gpd.GeoDataFrame) -> pd.DataFrame:
    rows = []
    sources = {(r, v): rasterio.open(WAVES / f"{r}_waves_{v}.tif") for r in REGIONS for v in ("Hs", "Tp")}
    try:
        geoms = sites.to_crs(sources[(REGIONS[0], "Hs")].crs)
        for site, geom in zip(geoms.site_id, geoms.geometry):
            out = {"site_id": site}
            for region in REGIONS:
                vals = {}
                for v in ("Hs", "Tp"):
                    try:
                        arr, _ = mask(sources[(region, v)], [geom], crop=True, filled=False)
                    except ValueError:  # polygon outside this raster
                        arr = np.ma.masked_all((1, 1))
                    a = arr.compressed()
                    vals[v] = a[np.isfinite(a)]
                if len(vals["Hs"]) and len(vals["Tp"]):
                    out.update({"region": region, "n_pixels": len(vals["Hs"]), "hs_mean": float(vals["Hs"].mean()),
                                "hs_p90": float(np.quantile(vals["Hs"], 0.9)), "tp_mean": float(vals["Tp"].mean()),
                                "power_mean": float(np.mean(vals["Hs"] ** 2) * vals["Tp"].mean())})
                    break
            rows.append(out)
    finally:
        for s in sources.values():
            s.close()
    return pd.DataFrame(rows)


def main(out: Path) -> None:
    if out.exists():
        raise FileExistsError(f"Refusing to overwrite {out}")
    sites = gpd.read_file(SITES)
    sites["GBRMPA_ID"] = sites.GBRMPA_ID.astype(str)

    # 1. Waves.
    w = site_waves(sites).merge(sites[["site_id", "GBRMPA_ID", "reef_name", "area", "x_coord", "y_coord"]], on="site_id")
    valid = w.dropna(subset=["hs_mean"])
    reef_w = valid.groupby("GBRMPA_ID").apply(lambda g: pd.Series({
        "hs_mean": np.average(g.hs_mean, weights=g.area), "tp_mean": np.average(g.tp_mean, weights=g.area),
        "power_mean": np.average(g.power_mean, weights=g.area), "hs_max_site": g.hs_mean.max(),
        "power_max_site": g.power_mean.max(), "sites_with_waves": len(g)}), include_groups=False)
    for c in ("hs_mean", "power_mean", "hs_max_site", "power_max_site"):
        reef_w[f"{c}_pct_rank"] = reef_w[c].rank(pct=True)

    # 2. Chlorophyll.
    d = pd.read_csv(INWATER, skiprows=26, low_memory=False)
    d["CHL_A"] = pd.to_numeric(d.CHL_A, errors="coerce")
    st = d[d.CHL_A > 0].groupby(["STATION_NAME", "TRANSECT"]).agg(lat=("LATITUDE", "mean"), lon=("LONGITUDE", "mean"),
                                                                 n=("CHL_A", "size"), chl_median=("CHL_A", "median")).reset_index()
    st = st[st.n >= 100]
    cent = sites.groupby("GBRMPA_ID")[["x_coord", "y_coord"]].mean()
    idx = pd.read_csv(CHL_INDEX, dtype={"RME_UNIQUE_ID": str})
    idx = idx[(idx.conversion == "rma") & (idx.year >= 2003)]
    uid = sites.drop_duplicates("GBRMPA_ID").set_index("GBRMPA_ID").UNIQUE_ID.astype(str)
    sat = idx.groupby("RME_UNIQUE_ID").p90_iw.median()
    chl_rows = []
    for reef in [r for g in GROUPS.values() for r in g]:
        lat, lon = cent.loc[reef, "y_coord"], cent.loc[reef, "x_coord"]
        dist = haversine_km(lat, lon, st.lat.to_numpy(), st.lon.to_numpy())
        k = int(np.argmin(dist))
        chl_rows.append({"GBRMPA_ID": reef, "nearest_station": st.STATION_NAME.iloc[k], "station_transect": st.TRANSECT.iloc[k],
                         "station_distance_km": float(dist[k]), "inwater_chl_median": float(st.chl_median.iloc[k]),
                         "satellite_p90_iw_median": float(sat.get(uid.loc[reef], np.nan))})
    chl = pd.DataFrame(chl_rows).set_index("GBRMPA_ID")

    # 3. Connectivity.
    shape = gpd.read_file(OWEN / "reefShapefiles" / "gbrShapeLL.shp")
    ids = [re.search(r"\(([^()]*)\)\s*$", str(n)).group(1).strip() for n in shape.reefName]
    index = {r: k for k, r in enumerate(ids)}
    in_domain = np.zeros(len(ids), bool)
    in_domain[[index[r] for r in sites.GBRMPA_ID.unique()]] = True
    idl = pd.read_csv(RME / "id" / "id_list_2024_12_01.csv", comment="#", header=None, dtype={0: str})
    area_old = {str(r.iloc[0]): float(r.iloc[1]) * 100.0 * (1.0 - float(r.iloc[2])) for _, r in idl.iterrows()}
    area = np.array([area_old.get(r, area_old.get(ALIASES.get(r, ""), 0.0)) for r in ids])
    mats = sorted(OWEN.glob("conMatCots*.csv"))
    reefs_all = [r for g in GROUPS.values() for r in g]
    conn = {r: np.zeros(4) for r in reefs_all}
    for path in mats:
        c = np.loadtxt(path, delimiter=",", dtype=np.float64)
        c[~np.isfinite(c).all(axis=1)] = 0.0
        for r in reefs_all:
            j = index[r]
            arrivals = area * c[:, j] / area[j]
            own = arrivals[j]
            dom = arrivals[in_domain].sum() - own
            ext = arrivals[~in_domain].sum()
            conn[r] += np.array([own + dom + ext, own, dom, ext]) / len(mats)
    con = pd.DataFrame(conn, index=["arrivals_equal_production", "own", "other_domain", "external"]).T
    for c in ("own", "other_domain", "external"):
        con[f"{c}_share"] = con[c] / con.arrivals_equal_production
    # Domain-wide percentile of total arrivals.
    allc = np.zeros(len(ids))
    for path in mats:
        c = np.loadtxt(path, delimiter=",", dtype=np.float64)
        c[~np.isfinite(c).all(axis=1)] = 0.0
        with np.errstate(divide="ignore", invalid="ignore"):
            allc += np.where(area > 0, (area[:, None] * c).sum(axis=0) / area, np.nan) / len(mats)
    dom_vals = pd.Series(allc[in_domain], index=[ids[k] for k in np.where(in_domain)[0]])
    con["arrivals_pct_rank_domain"] = [float((dom_vals <= con.loc[r, "arrivals_equal_production"]).mean()) for r in con.index]

    table = con.join(reef_w, how="left").join(chl, how="left")
    table.insert(0, "group", ["absence" if r in GROUPS["absence"] else "outbreak" for r in table.index])
    table.insert(0, "reef", [LABEL[r] for r in table.index])
    out.mkdir(parents=True)
    table.to_csv(out / "never_cots_reef_comparison.csv")
    w.to_csv(out / "site_waves.csv", index=False)

    fig, axes = plt.subplots(1, 4, figsize=(17, 4.6), sharey=True)
    order = table.sort_values(["group", "reef"]).index
    y = np.arange(len(order))
    colors = ["#eb6834" if table.loc[r, "group"] == "absence" else "#0d366b" for r in order]
    panels = (("arrivals_equal_production", "Larval arrivals per ha\n(equal production everywhere)"),
              ("hs_mean", "Mean significant wave height (m)"), ("power_max_site", "Wave power index, most exposed site\n(Hs² Tp)"),
              ("inwater_chl_median", "In-water chl-a, nearest LTCMP station (µg/L)"))
    for ax, (col, title) in zip(axes, panels):
        ax.barh(y, table.loc[order, col], color=colors, height=0.6)
        ax.set_title(title, loc="left", fontsize=9)
        ax.grid(axis="x", color="#e4e3df", linewidth=.8)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_yticks(y)
    axes[0].set_yticklabels([table.loc[r, "reef"] for r in order])
    axes[0].text(0.98, 0.02, "orange: never-COTS reefs", transform=axes[0].transAxes, ha="right", fontsize=8, color="#eb6834")
    fig.tight_layout()
    fig.savefig(out / "never_cots_comparison.png", dpi=150, facecolor="#fcfcfb")
    plt.close(fig)

    meta = {"status": "diagnostic", "groups": GROUPS,
            "inputs": {"sites": sha256(SITES), "inwater": sha256(INWATER), "chl_index": sha256(CHL_INDEX),
                       **{f"{r}_waves_{v}": sha256(WAVES / f"{r}_waves_{v}.tif") for r in REGIONS for v in ("Hs", "Tp")},
                       **{p.name: sha256(p) for p in mats}},
            "script_sha256": sha256(Path(__file__)),
            "outputs": {n: sha256(out / n) for n in ("never_cots_reef_comparison.csv", "site_waves.csv", "never_cots_comparison.png")}}
    (out / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
    pd.set_option("display.width", 250)
    print(table.drop(columns=["own", "other_domain", "external"]).round(3).to_string())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    main(parser.parse_args().output.resolve())
