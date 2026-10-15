"""VIIRS wet-season chlorophyll-a for the Lizard v0.3 domain reefs, calibrated to MODIS-Aqua,
to extend the reef-based chl-a index to 2025-2026 (lead's decision, 2026-10-15, tenth round).

Source: NOAA CoastWatch ERDDAP `noaacwNPPVIIRSSQchlaDaily` (S-NPP VIIRS, Science Quality,
global 4 km, daily, 2012-present; public). Same reef pixel rule and season as
extract_chla_wet_season.py (nearest cell with any valid value in the season within 5 km of
the reef centroid; 1 Dec (year-1) - 30 Apr (year); n >= 10 days). VIIRS log anomalies are
computed against each reef's VIIRS geometric-mean climatology over the MODIS overlap
(2013-2024) and compared with the MODIS log anomalies over those seasons. The combined index
uses MODIS for 2003-2024 and, for 2025-2026, the VIIRS anomaly mapped to the MODIS scale by
the pooled least-squares slope (no intercept) fitted on the overlap.

    python extract_chla_viirs_wet_season.py [--output <dir>]
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import time
import urllib.request
import warnings
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr

ROOT = Path(__file__).resolve().parents[2]
GPKG = ROOT / "sandbox" / "data" / "Lizard_Historical_v0.3" / "spatial" / "lizard_cluster.gpkg"
MODIS = ROOT / "sandbox" / "data" / "derived" / "chla_wet_season_lizard" / "chla_wet_season_reef_year.csv"
ERDDAP = "https://coastwatch.noaa.gov/erddap/griddap/noaacwNPPVIIRSSQchlaDaily"
DEFAULT_OUT = ROOT / "sandbox" / "data" / "derived" / "chla_index_lizard_2003_2026"
FIRST, LAST, OVERLAP = 2013, 2026, (2013, 2024)
RADIUS_KM = 5.0


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def haversine_km(lat0, lon0, lats, lons):
    lat0_r, lats_r = np.radians(lat0), np.radians(lats)
    dlat, dlon = np.radians(np.asarray(lats) - lat0), np.radians(np.asarray(lons) - lon0)
    a = np.sin(dlat / 2) ** 2 + np.cos(lat0_r) * np.cos(lats_r) * np.sin(dlon / 2) ** 2
    return 6371.0 * 2 * np.arcsin(np.sqrt(a))


def fetch_season(year: int, box: tuple[float, float, float, float]) -> xr.DataArray:
    lat_min, lat_max, lon_min, lon_max = box
    query = (f"{ERDDAP}.nc?chlor_a[({year - 1}-12-01T00:00:00Z):1:({year}-04-30T23:59:59Z)][(0.0)]"
             f"[({lat_max}):1:({lat_min})][({lon_min}):1:({lon_max})]")
    for attempt in range(4):
        try:
            with urllib.request.urlopen(query, timeout=600) as response:
                payload = response.read()
            break
        except Exception:
            if attempt == 3:
                raise
            time.sleep(10 * (attempt + 1))
    ds = xr.open_dataset(io.BytesIO(payload), engine="h5netcdf" if payload[:4] == b"\x89HDF" else "scipy")
    return ds["chlor_a"].isel(altitude=0).load()


def main(out: Path) -> None:
    if out.exists():
        raise FileExistsError(f"Refusing to overwrite {out}")
    sites = gpd.read_file(GPKG)[["UNIQUE_ID", "GBRMPA_ID", "x_coord", "y_coord"]]
    reefs = sites.groupby("UNIQUE_ID").agg(reef_id=("GBRMPA_ID", "first"), lon=("x_coord", "mean"),
                                           lat=("y_coord", "mean")).reset_index()
    pad = RADIUS_KM / 111.0 + 0.05
    box = (float(reefs.lat.min() - pad), float(reefs.lat.max() + pad),
           float(reefs.lon.min() - pad), float(reefs.lon.max() + pad))
    rows = []
    for year in range(FIRST, LAST + 1):
        arr = fetch_season(year, box)
        values = arr.values.astype(float)
        values[~np.isfinite(values) | (values <= 0)] = np.nan
        valid = np.isfinite(values).any(axis=0)
        lat2, lon2 = np.meshgrid(arr.latitude.values, arr.longitude.values, indexing="ij")
        for r in reefs.itertuples(index=False):
            d = haversine_km(r.lat, r.lon, lat2, lon2)
            ok = valid & (d <= RADIUS_KM)
            if not ok.any():
                rows.append({"RME_UNIQUE_ID": str(r.UNIQUE_ID), "year": year, "n": 0})
                continue
            i, j = np.unravel_index(np.where(ok, d, np.inf).argmin(), ok.shape)
            series = values[:, i, j]
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                rows.append({"RME_UNIQUE_ID": str(r.UNIQUE_ID), "year": year, "n": int(np.isfinite(series).sum()),
                             "mean": float(np.nanmean(series)), "p90": float(np.nanquantile(series, 0.9)),
                             "pixel_distance_km": float(d[i, j])})
        print(f"{year}: {values.shape[0]} days", flush=True)
    viirs = pd.DataFrame(rows)
    viirs.loc[viirs.n < 10, ["mean", "p90"]] = np.nan
    base = viirs[viirs.year.between(*OVERLAP)]
    clim = base.groupby("RME_UNIQUE_ID")[["mean", "p90"]].agg(lambda v: float(np.exp(np.nanmean(np.log(v)))))
    viirs = viirs.join(clim.rename(columns={"mean": "clim_mean", "p90": "clim_p90"}), on="RME_UNIQUE_ID")
    viirs["log_anom_mean"] = np.log(viirs["mean"] / viirs.clim_mean)
    viirs["log_anom_p90"] = np.log(viirs["p90"] / viirs.clim_p90)

    modis = pd.read_csv(MODIS, dtype={"RME_UNIQUE_ID": str})
    both = modis.merge(viirs, on=["RME_UNIQUE_ID", "year"], suffixes=("_modis", "_viirs"))
    both = both[both.year.between(*OVERLAP)]
    calibration = {}
    for k in ("mean", "p90"):
        x, y = both[f"log_anom_{k}_viirs"], both[f"log_anom_{k}_modis"]
        ok = np.isfinite(x) & np.isfinite(y)
        slope = float((x[ok] * y[ok]).sum() / (x[ok] ** 2).sum())
        dom = both[ok].groupby("year")[[f"log_anom_{k}_viirs", f"log_anom_{k}_modis"]].median()
        calibration[k] = {"slope_no_intercept": slope, "reef_season_r": float(np.corrcoef(x[ok], y[ok])[0, 1]),
                          "domain_median_r": float(np.corrcoef(dom.iloc[:, 0], dom.iloc[:, 1])[0, 1]),
                          "n": int(ok.sum())}

    index = modis[["RME_UNIQUE_ID", "reef_id", "year", "log_anom_mean", "log_anom_p90", "n"]].copy()
    index["source"] = np.where(np.isfinite(index.log_anom_mean), "MODIS-Aqua OC3", "missing")
    extra = viirs[viirs.year > OVERLAP[1]][["RME_UNIQUE_ID", "year", "log_anom_mean", "log_anom_p90", "n"]].copy()
    extra["log_anom_mean"] *= calibration["mean"]["slope_no_intercept"]
    extra["log_anom_p90"] *= calibration["p90"]["slope_no_intercept"]
    extra = extra.merge(reefs.rename(columns={"UNIQUE_ID": "RME_UNIQUE_ID"})[["RME_UNIQUE_ID", "reef_id"]], on="RME_UNIQUE_ID")
    extra["source"] = "VIIRS S-NPP SQ (MODIS-calibrated)"
    index = pd.concat([index[index.year <= OVERLAP[1]], extra], ignore_index=True).sort_values(["RME_UNIQUE_ID", "year"])

    out.mkdir(parents=True)
    viirs.to_csv(out / "viirs_wet_season_reef_year.csv", index=False)
    index.to_csv(out / "chla_index_reef_year.csv", index=False)
    meta = {"status": "derived_input_not_promoted", "viirs_source": ERDDAP, "modis_input": str(MODIS.relative_to(ROOT)),
            "modis_input_sha256": sha256(MODIS), "season": "1 Dec (year-1) - 30 Apr (year)", "overlap": OVERLAP,
            "calibration": calibration, "index": "MODIS log anomalies 2003-2024; VIIRS x slope for 2025-2026",
            "inputs": {"gpkg_sha256": sha256(GPKG)}, "script_sha256": sha256(Path(__file__)),
            "outputs": {f.name: sha256(f) for f in sorted(out.iterdir()) if f.suffix == ".csv"}}
    (out / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps(calibration, indent=1))
    print(index.groupby(["year", "source"])[["log_anom_mean", "log_anom_p90"]].median().round(3).to_string())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    main(parser.parse_args().output.resolve())
