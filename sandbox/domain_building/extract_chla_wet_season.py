"""Reef-based wet-season chlorophyll-a for the Lizard v0.3 domain reefs (lead's decision,
2026-10-15, ninth round: reef-based wet-season chl-a anomaly plus the wet-season 90th
percentile as an indicator of the maximum; index for the initiation-region larval-survival
hypothesis).

Method (as BleachingMortality/src/data/extract_sst_chla_features.py): IMOS SRS MODIS-Aqua
daily chlorophyll-a, OC3 algorithm, 4 km (public GBR DMS zarr, anonymous). Each reef
centroid (mean of its model-site centroids) uses the nearest cell with any valid value in
that season within 5 km. Wet season = 1 December (year - 1) to 30 April (year), labelled by
`year` (the COTS spawning and larval season). Per reef-season: n valid days, mean, median,
90th percentile (positive values only). Anomalies on the log scale against each reef's
2003-2026 mean of the seasonal values:
    log_anom_mean = log(mean / reef climatology of means), likewise for p90.

    python extract_chla_wet_season.py [--output <dir>] [--first-year 2003] [--last-year 2026]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import warnings
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr

ROOT = Path(__file__).resolve().parents[2]
GPKG = ROOT / "sandbox" / "data" / "Lizard_Historical_v0.3" / "spatial" / "lizard_cluster.gpkg"
STORE = "s3://gbr-dms-data-public/imos-srs-aqua-oc-chla-4km/data.zarr"
VARIABLE = "chl_oc3"
DEFAULT_OUT = ROOT / "sandbox" / "data" / "derived" / "chla_wet_season_lizard"
RADIUS_KM = 5.0


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def haversine_km(lat0, lon0, lats, lons):
    lat0_r, lats_r = np.radians(lat0), np.radians(lats)
    dlat, dlon = np.radians(np.asarray(lats) - lat0), np.radians(np.asarray(lons) - lon0)
    a = np.sin(dlat / 2) ** 2 + np.cos(lat0_r) * np.cos(lats_r) * np.sin(dlon / 2) ** 2
    return 6371.0 * 2 * np.arcsin(np.sqrt(a))


def nearest_valid(valid, glat, glon, lat, lon):
    lat2, lon2 = np.meshgrid(glat, glon, indexing="ij")
    d = haversine_km(lat, lon, lat2, lon2)
    ok = valid & (d <= RADIUS_KM)
    if not ok.any():
        return None
    i, j = np.unravel_index(np.where(ok, d, np.inf).argmin(), ok.shape)
    return int(i), int(j), float(d[i, j])


def main(out: Path, first: int, last: int) -> None:
    if out.exists():
        raise FileExistsError(f"Refusing to overwrite {out}")
    sites = gpd.read_file(GPKG)[["UNIQUE_ID", "GBRMPA_ID", "x_coord", "y_coord"]]
    reefs = sites.groupby("UNIQUE_ID").agg(reef_id=("GBRMPA_ID", "first"), lon=("x_coord", "mean"),
                                           lat=("y_coord", "mean")).reset_index()
    ds = xr.open_zarr(STORE, storage_options={"anon": True})
    pad = RADIUS_KM / 111.0 + 0.05
    lat_vals = ds.latitude.values
    lat_slice = (slice(reefs.lat.max() + pad, reefs.lat.min() - pad) if lat_vals[0] > lat_vals[-1]
                 else slice(reefs.lat.min() - pad, reefs.lat.max() + pad))
    box = ds[VARIABLE].sel(latitude=lat_slice, longitude=slice(reefs.lon.min() - pad, reefs.lon.max() + pad))
    coverage = (str(ds.time.values[0])[:10], str(ds.time.values[-1])[:10])
    rows = []
    for year in range(first, last + 1):
        arr = box.sel(time=slice(f"{year - 1}-12-01", f"{year}-04-30")).compute()
        values = arr.values.astype(float)
        values[values <= 0] = np.nan
        n_days = values.shape[0]
        valid = np.isfinite(values).any(axis=0)
        for r in reefs.itertuples(index=False):
            pick = nearest_valid(valid, arr.latitude.values, arr.longitude.values, r.lat, r.lon)
            if pick is None:
                rows.append({"RME_UNIQUE_ID": str(r.UNIQUE_ID), "reef_id": r.reef_id, "year": year,
                             "season_days": n_days, "n": 0, "mean": np.nan, "median": np.nan, "p90": np.nan,
                             "pixel_distance_km": np.nan})
                continue
            i, j, dist = pick
            series = values[:, i, j]
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                rows.append({"RME_UNIQUE_ID": str(r.UNIQUE_ID), "reef_id": r.reef_id, "year": year,
                             "season_days": n_days, "n": int(np.isfinite(series).sum()),
                             "mean": float(np.nanmean(series)), "median": float(np.nanmedian(series)),
                             "p90": float(np.nanquantile(series, 0.9)), "pixel_distance_km": dist})
        print(f"{year}: {n_days} days; domain median of reef means "
              f"{np.nanmedian([x['mean'] for x in rows if x['year'] == year]):.3f} mg m-3", flush=True)
    table = pd.DataFrame(rows)
    usable = table.n >= 10
    table.loc[~usable, ["mean", "median", "p90"]] = np.nan
    clim = table.groupby("RME_UNIQUE_ID")[["mean", "p90"]].agg(lambda v: float(np.exp(np.nanmean(np.log(v)))))
    table = table.join(clim.rename(columns={"mean": "clim_mean", "p90": "clim_p90"}), on="RME_UNIQUE_ID")
    table["log_anom_mean"] = np.log(table["mean"] / table.clim_mean)
    table["log_anom_p90"] = np.log(table["p90"] / table.clim_p90)

    out.mkdir(parents=True)
    table.to_csv(out / "chla_wet_season_reef_year.csv", index=False)
    meta = {"status": "derived_input_not_promoted", "store": STORE, "variable": VARIABLE,
            "store_time_coverage": coverage, "season": "1 Dec (year-1) to 30 Apr (year)",
            "pixel": f"nearest cell with any valid value in the season within {RADIUS_KM} km of the reef centroid",
            "usable_rule": "n >= 10 valid days, otherwise NaN", "climatology": f"geometric mean over {first}-{last}",
            "years": [first, last], "reefs": int(len(reefs)),
            "missing_reef_seasons": int(table["mean"].isna().sum()),
            "inputs": {"gpkg_sha256": sha256(GPKG)}, "script_sha256": sha256(Path(__file__)),
            "outputs": {"chla_wet_season_reef_year.csv": sha256(out / "chla_wet_season_reef_year.csv")}}
    (out / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
    summary = table.groupby("year")[["mean", "p90", "log_anom_mean", "log_anom_p90", "n"]].median()
    print(summary.round(3).to_string())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--first-year", type=int, default=2003)
    parser.add_argument("--last-year", type=int, default=2026)
    a = parser.parse_args()
    main(a.output.resolve(), a.first_year, a.last_year)
