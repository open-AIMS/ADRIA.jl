"""Annual maximum NOAA CRW DHW for the Lizard v0.3 domain reefs, 2016-2026, to extend the
historical forcing beyond 2024 (lead's decision, 2026-10-14, seventh round).

Method (as BleachingMortality/src/data/add_repeated_dhw_exposure.py): NOAA Coral Reef Watch
daily 5 km DHW from the public GBR DMS (s3://gbr-dms-data-public/noaa-crw-chs-dhw/data.zarr,
anonymous), nearest cell to each reef centroid (mean of the reef's model-site centroids),
calendar-year maximum. Years 2016-2024 are extracted only to validate the method against the
domain's existing `sandbox/data/dhw_historical.csv`; 2025-2026 are the extension.

Output (derived data; not yet a domain): sandbox/data/derived/dhw_2016_2026_lizard/

    python extract_dhw_2025_2026.py [--output <dir>]
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "sandbox" / "data"
GPKG = DATA / "Lizard_Historical_v0.3" / "spatial" / "lizard_cluster.gpkg"
HISTORICAL = DATA / "dhw_historical.csv"
STORE = "s3://gbr-dms-data-public/noaa-crw-chs-dhw/data.zarr"
DEFAULT_OUT = DATA / "derived" / "dhw_2016_2026_lizard"
YEARS = list(range(2016, 2027))
VALIDATION_YEARS = list(range(2016, 2025))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main(out: Path) -> None:
    if out.exists():
        raise FileExistsError(f"Refusing to overwrite {out}")
    sites = gpd.read_file(GPKG)[["UNIQUE_ID", "GBRMPA_ID", "x_coord", "y_coord"]]
    reefs = sites.groupby("UNIQUE_ID").agg(reef_id=("GBRMPA_ID", "first"), lon=("x_coord", "mean"),
                                           lat=("y_coord", "mean")).reset_index()
    dataset = xr.open_zarr(STORE, storage_options={"anon": True})
    lats = xr.DataArray(reefs.lat.to_numpy(), dims="reef")
    lons = xr.DataArray(reefs.lon.to_numpy(), dims="reef")
    rows = []
    for year in YEARS:
        values = (dataset["degree_heating_week"].sel(time=slice(f"{year}-01-01", f"{year}-12-31"))
                  .sel(lat=lats, lon=lons, method="nearest"))
        n_days = int(values.sizes["time"])
        annual = values.max(dim="time").compute().values
        for uid, rid, v in zip(reefs.UNIQUE_ID, reefs.reef_id, annual):
            rows.append({"RME_UNIQUE_ID": str(uid), "reef_id": rid, "timestep": year, "dhw": float(v), "days": n_days})
        print(f"{year}: {n_days} days, domain mean max DHW {np.nanmean(annual):.2f}", flush=True)
    extracted = pd.DataFrame(rows)
    if not np.all(np.isfinite(extracted.dhw)):
        raise ValueError("Non-finite DHW (masked cell at a reef centroid)")

    hist = pd.read_csv(HISTORICAL, dtype={"RME_UNIQUE_ID": str}).drop_duplicates(["RME_UNIQUE_ID", "timestep"])
    check = extracted[extracted.timestep.isin(VALIDATION_YEARS)].merge(
        hist, on=["RME_UNIQUE_ID", "timestep"], suffixes=("", "_historical"))
    diff = check.dhw - check.dhw_historical
    validation = {"reef_years": int(len(check)), "max_abs_diff": float(diff.abs().max()),
                  "mean_abs_diff": float(diff.abs().mean()), "mean_diff": float(diff.mean()),
                  "share_within_0_1": float((diff.abs() <= 0.1).mean()),
                  "correlation": float(np.corrcoef(check.dhw, check.dhw_historical)[0, 1]),
                  "by_year_mean_diff": check.assign(d=diff).groupby("timestep").d.mean().round(3).to_dict()}

    out.mkdir(parents=True)
    extracted.to_csv(out / "dhw_reef_year_2016_2026.csv", index=False)
    check.assign(diff=diff).to_csv(out / "validation_2016_2024.csv", index=False)
    meta = {"status": "derived_input_not_promoted", "store": STORE, "variable": "degree_heating_week",
            "method": "calendar-year maximum of daily DHW at the nearest 5 km cell to each reef centroid",
            "reef_centroid": "mean of model-site centroids (lizard_cluster.gpkg x_coord, y_coord) per UNIQUE_ID",
            "years": YEARS, "validation_against_dhw_historical_csv": validation,
            "inputs": {"gpkg_sha256": sha256(GPKG), "dhw_historical_sha256": sha256(HISTORICAL)},
            "script_sha256": sha256(Path(__file__)),
            "outputs": {f.name: sha256(f) for f in sorted(out.iterdir())}}
    (out / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps(validation, indent=1))
    print(extracted[extracted.timestep >= 2023].groupby("timestep").dhw.describe().round(2).to_string())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    main(parser.parse_args().output.resolve())
