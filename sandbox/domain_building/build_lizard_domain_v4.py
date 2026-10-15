"""Build the additive Lizard_Historical_v0.4 domain package (1985-2026) from v0.3.

Lead's decisions (2026-10-14/15): extend to 2026; DHW 2025-2026 from the DMS (NOAA CRW,
validated against the domain's 2016-2024 DHW, extract_dhw_2025_2026.py); cyclones ignored
for 2025-2026 (zero mortality); 2026 outbreak predictions from gbrPreds2026Adj_20262408.csv.

Changes relative to v0.3 (everything else is copied byte-for-byte):
* DHWs/dhw_RCPhistorical.nc and dhw_RCP45.nc: 1985-2026; 1985-2024 copied from v0.3,
  2025-2026 the reef annual-maximum DHW expanded to sites by UNIQUE_ID;
* cyclones/cyclone_mortality.csv: columns 2025 and 2026 set to 0 (labelled stand-in);
* datapackage.json: domain_data_version 0.4.0.
The COTS connectivity datasets for 2026 are written afterwards by
build_lizard_cots_connectivity_owen.py --domain <v0.4> (two datasets: the lead's predictions
extended with 2026, and the cull/SALAD-corrected layer).

    python build_lizard_domain_v4.py [--output <dir>]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from scipy.io import netcdf_file

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA = REPO_ROOT / "sandbox" / "data"
SOURCE = DATA / "Lizard_Historical_v0.3"
DEFAULT_OUTPUT = DATA / "Lizard_Historical_v0.4"
DHW_EXTENSION = DATA / "derived" / "dhw_2016_2026_lizard" / "dhw_reef_year_2016_2026.csv"
YEARS = list(range(1985, 2027))
NEW_YEARS = (2025, 2026)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main(output: Path) -> None:
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite {output}")
    shutil.copytree(SOURCE, output, ignore=shutil.ignore_patterns("datasets"))

    sites = gpd.read_file(SOURCE / "spatial" / "lizard_cluster.gpkg")[["site_id", "UNIQUE_ID"]]
    with netcdf_file(SOURCE / "DHWs" / "dhw_RCPhistorical.nc", "r", mmap=False) as ds:
        old_years = ds.variables["time"][:].astype(int).tolist()
        old = np.array(ds.variables["dhw"][:], dtype=np.float32)  # scenario x location x time
    if old_years != list(range(1985, 2025)) or old.shape[1] != len(sites):
        raise ValueError("Unexpected v0.3 DHW layout")
    ext = pd.read_csv(DHW_EXTENSION, dtype={"RME_UNIQUE_ID": str})
    ext = ext[ext.timestep.isin(NEW_YEARS)].set_index(["RME_UNIQUE_ID", "timestep"]).dhw
    new = np.empty((old.shape[0], old.shape[1], len(NEW_YEARS)), dtype=np.float32)
    for site_index, uid in enumerate(sites.UNIQUE_ID.astype(str)):
        for k, year in enumerate(NEW_YEARS):
            new[:, site_index, k] = ext.loc[(uid, year)]
    dhw = np.concatenate([old, new], axis=2)
    if not np.all(np.isfinite(dhw)) or np.any(dhw < 0.0):
        raise ValueError("Extended DHW cube is invalid")
    for name in ("dhw_RCPhistorical.nc", "dhw_RCP45.nc"):
        path = output / "DHWs" / name
        path.unlink()
        with netcdf_file(path, "w") as dataset:
            dataset.createDimension("time", dhw.shape[2])
            dataset.createDimension("location", dhw.shape[1])
            dataset.createDimension("scenario", dhw.shape[0])
            dataset.createVariable("time", "i4", ("time",))[:] = YEARS
            dataset.createVariable("location", "i4", ("location",))[:] = np.arange(1, dhw.shape[1] + 1)
            dataset.createVariable("scenario", "i4", ("scenario",))[:] = np.arange(1, dhw.shape[0] + 1)
            variable = dataset.createVariable("dhw", "f4", ("scenario", "location", "time"))
            variable[:] = dhw
            variable.units = "degree heating weeks"
            dataset.history = "Historical 1985-2024 (v0.3) plus NOAA CRW 2025-2026 reef annual maxima expanded to V2 sites"

    cyclone_path = output / "cyclones" / "cyclone_mortality.csv"
    cyclones = pd.read_csv(cyclone_path)
    for year in NEW_YEARS:
        cyclones[str(year)] = 0.0
    cyclones.to_csv(cyclone_path, index=False)

    package = json.loads((SOURCE / "datapackage.json").read_text())
    package["domain_data_version"] = "0.4.0"
    package["description"] = package["description"] + "; extended to 2026 (DHW 2025-2026 NOAA CRW; cyclones 2025-2026 zero)"
    (output / "datapackage.json").write_text(json.dumps(package, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    provenance = {
        "domain_version": "Lizard_Historical_v0.4", "built": datetime.now(timezone.utc).isoformat(),
        "builder": "sandbox/domain_building/build_lizard_domain_v4.py", "builder_sha256": sha256_file(Path(__file__)),
        "source": "Lizard_Historical_v0.3 (copied except the changes below; cots_connectivity datasets rebuilt separately)",
        "changes": {
            "dhw": {"years": [YEARS[0], YEARS[-1]], "extension_source": str(DHW_EXTENSION.relative_to(REPO_ROOT)),
                    "extension_sha256": sha256_file(DHW_EXTENSION),
                    "validation": "extract_dhw_2025_2026.py reproduces v0.3 2016-2024 (r 0.9996, mean |diff| 0.018)"},
            "cyclones": "2025 and 2026 columns zero (lead: cyclones not influential in 2025-26; labelled stand-in)",
        },
        "outputs": {str(p.relative_to(output)): sha256_file(p) for p in (
            output / "DHWs" / "dhw_RCPhistorical.nc", output / "DHWs" / "dhw_RCP45.nc", cyclone_path)},
    }
    (output / "provenance_v0.4.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(json.dumps(provenance, indent=1))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    main(parser.parse_args().output.resolve())
