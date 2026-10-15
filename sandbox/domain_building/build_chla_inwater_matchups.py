"""Same-day matchups between AIMS in-water chlorophyll-a (Long Term Chlorophyll Monitoring
Program, 1992-2009) and MODIS-Aqua OC3, to put the satellite record on the in-water scale
(lead's decision, 2026-10-15, thirteenth round: in-water data are the reference because surface
satellite chl-a under-represents the chl-a larvae are exposed to).

In-water: samples north of 17.5 S from July 2002 (MODIS start) to January 2009; per station
visit (station x date) the mean CHL_A over all sampled depths (larval-exposure proxy) and, as a
sensitivity, over samples <= 10 m. Satellite: IMOS MODIS-Aqua OC3 (DMS zarr) on the same UTC
day at the nearest valid pixel within 5 km of the station. Fit (OLS):
    log10(in-water) = a + b * log10(OC3)

    python build_chla_inwater_matchups.py [--output <dir>]
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

ROOT = Path(__file__).resolve().parents[2]
INWATER = ROOT / "sandbox" / "data" / "GBR_chlorophyll_monitoring_1992-2009.csv"
STORE = "s3://gbr-dms-data-public/imos-srs-aqua-oc-chla-4km/data.zarr"
DEFAULT_OUT = ROOT / "sandbox" / "data" / "derived" / "chla_inwater_matchups"
RADIUS_KM = 5.0
BOX = (-17.6, -10.4, 142.0, 147.6)


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


def fit(x: np.ndarray, y: np.ndarray) -> dict:
    lx, ly = np.log10(x), np.log10(y)
    b, a = np.polyfit(lx, ly, 1)
    resid = ly - (a + b * lx)
    return {"a": float(a), "b": float(b), "n": int(len(x)), "r": float(np.corrcoef(lx, ly)[0, 1]),
            "resid_sd_log10": float(resid.std(ddof=2)),
            "median_ratio_inwater_to_oc3": float(np.median(y / x))}


def main(out: Path) -> None:
    if out.exists():
        raise FileExistsError(f"Refusing to overwrite {out}")
    d = pd.read_csv(INWATER, skiprows=26, low_memory=False)
    d["date"] = pd.to_datetime(d.SAMPLE_TIME, dayfirst=True, errors="coerce")
    d["CHL_A"] = pd.to_numeric(d.CHL_A, errors="coerce")
    d["SAMPLE_DEPTH"] = pd.to_numeric(d.SAMPLE_DEPTH, errors="coerce")
    d = d[(d.LATITUDE > BOX[0]) & (d.date >= "2002-07-01") & (d.CHL_A > 0)].copy()
    d["day"] = d.date.dt.floor("D")
    visits = (d.groupby(["STATION_NAME", "TRANSECT", "day"])
              .agg(lat=("LATITUDE", "mean"), lon=("LONGITUDE", "mean"), n=("CHL_A", "size"),
                   inwater_all=("CHL_A", "mean"),
                   inwater_shallow=("CHL_A", lambda v: v[d.loc[v.index, "SAMPLE_DEPTH"] <= 10].mean()))
              .reset_index())
    ds = xr.open_zarr(STORE, storage_options={"anon": True})
    lat_vals = ds.latitude.values
    lat_slice = slice(BOX[1], BOX[0]) if lat_vals[0] > lat_vals[-1] else slice(BOX[0], BOX[1])
    box = ds["chl_oc3"].sel(latitude=lat_slice, longitude=slice(BOX[2], BOX[3]))
    times = pd.to_datetime(box.time.values).floor("D")
    rows = []
    for day, g in visits.groupby("day"):
        idx = np.where(times == day)[0]
        if len(idx) == 0:
            continue
        grid = box.isel(time=int(idx[0])).compute()
        values = grid.values.astype(float)
        values[values <= 0] = np.nan
        valid = np.isfinite(values)
        if not valid.any():
            continue
        lat2, lon2 = np.meshgrid(grid.latitude.values, grid.longitude.values, indexing="ij")
        for v in g.itertuples(index=False):
            dist = haversine_km(v.lat, v.lon, lat2, lon2)
            ok = valid & (dist <= RADIUS_KM)
            if not ok.any():
                continue
            i, j = np.unravel_index(np.where(ok, dist, np.inf).argmin(), ok.shape)
            rows.append(v._asdict() | {"oc3": float(values[i, j]), "pixel_distance_km": float(dist[i, j])})
    m = pd.DataFrame(rows)
    fits = {"all_depths": fit(m.oc3.to_numpy(), m.inwater_all.to_numpy())}
    s = m[np.isfinite(m.inwater_shallow)]
    fits["shallow_le_10m"] = fit(s.oc3.to_numpy(), s.inwater_shallow.to_numpy())
    fits["by_transect_all_depths"] = {t: fit(g.oc3.to_numpy(), g.inwater_all.to_numpy())
                                      for t, g in m.groupby("TRANSECT") if len(g) >= 15}
    wet = m[pd.to_datetime(m.day).dt.month.isin([12, 1, 2, 3, 4])]
    fits["wet_season_all_depths"] = fit(wet.oc3.to_numpy(), wet.inwater_all.to_numpy())
    out.mkdir(parents=True)
    m.to_csv(out / "matchups.csv", index=False)
    meta = {"status": "derived_input_not_promoted", "inwater": str(INWATER.relative_to(ROOT)), "inwater_sha256": sha256(INWATER),
            "satellite": STORE, "pixel": f"nearest valid same-day OC3 pixel within {RADIUS_KM} km",
            "visits": int(len(visits)), "matchups": int(len(m)), "fits": fits,
            "model": "log10(inwater) = a + b log10(OC3), OLS", "script_sha256": sha256(Path(__file__)),
            "outputs": {"matchups.csv": sha256(out / "matchups.csv")}}
    (out / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps({k: v for k, v in meta.items() if k in ("visits", "matchups", "fits")}, indent=1))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    main(parser.parse_args().output.resolve())
