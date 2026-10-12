"""Build the additive Lizard_Historical_v0.3 domain package from v0.2.

Changes relative to v0.2 (everything else is copied byte-for-byte):

* spatial/lizard_cluster.gpkg: ``k`` = site habitable fraction from the GBR10
  2025 benthic map, infilled with the geomorphic map
  (``k_geomorphic_infill`` from build_site_habitable_fraction.py), floored at
  ``K_FLOOR`` because ADRIA cannot run k = 0 sites. ``k_raw`` and ``k_source``
  record the unfloored value and its origin.
* initial_cover.nc: cover set to zero at sites whose raw habitable fraction is 0.
* cyclones/cyclone_mortality.csv: site x calendar-year coral mortality fraction
  from historical cyclone wave exposure (GBRMPA RRN ``cyc_maxHrs4mw``,
  build_cyclone_hours_rrn.py), ``1 - exp(B_PER_HOUR * hours)``, applied
  uniformly to all coral groups. ``datapackage.json`` declares this file so the
  ADRIA Lizard loader uses it; v0.1/v0.2 have no such key and keep zero cyclones.

Sources are internal or licensed; the package lives in the git-ignored
sandbox/data tree. Requires scipy + geopandas (base miniconda Python):

    python sandbox/domain_building/build_lizard_domain_v3.py
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
SOURCE = DATA / "Lizard_Historical_v0.2"
DEFAULT_OUTPUT = DATA / "Lizard_Historical_v0.3"
DERIVED = SOURCE / "derived"
HABITAT_CSV = DERIVED / "site_habitable_fraction_benthic_geomorphic2025.csv"
CYCLONE_CSV = DERIVED / "cyclone_hours4m_rrn2025.csv"
K_FLOOR = 0.01
B_PER_HOUR = -0.00444  # COTS-Coral-OutcomesV2 LTMP dose-response, -0.0444 per 10 h
YEARS = list(range(1985, 2025))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main(output: Path) -> None:
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite {output}")
    shutil.copytree(SOURCE, output, ignore=shutil.ignore_patterns("derived"))

    # Habitable fraction ---------------------------------------------------------
    gpkg = output / "spatial" / "lizard_cluster.gpkg"
    sites = gpd.read_file(gpkg)
    habitat = pd.read_csv(HABITAT_CSV, dtype={"site_id": str}).set_index("site_id")
    missing = sorted(set(sites.site_id) - set(habitat.index))
    if missing:
        raise ValueError(f"Habitat fraction missing for sites: {missing[:5]}")
    raw = habitat.loc[sites.site_id, "k_geomorphic_infill"].to_numpy(dtype=float)
    if not np.all((raw >= 0.0) & (raw <= 1.0)):
        raise ValueError("Habitable fraction outside [0, 1]")
    sites["k_raw"] = raw
    sites["k"] = np.maximum(raw, K_FLOOR)
    sites["k_source"] = np.where(raw == 0.0, "floored_zero", "benthic_geomorphic_2025")
    gpkg.unlink()
    sites.to_file(gpkg, layer="lizard_cluster", driver="GPKG")

    # Initial cover ----------------------------------------------------------------
    cover_path = output / "initial_cover.nc"
    with netcdf_file(cover_path, "r", mmap=False) as dataset:
        covers = dataset.variables["covers"][:].copy()  # location x species
        species = dataset.variables["species"][:].copy()
        location = dataset.variables["location"][:].copy()
        units = dataset.variables["covers"].units
    if covers.shape != (len(sites), 35):
        raise ValueError("Unexpected initial cover dimensions")
    zero_sites = raw == 0.0
    covers[zero_sites, :] = 0.0
    cover_path.unlink()
    with netcdf_file(cover_path, "w") as dataset:
        dataset.createDimension("species", covers.shape[1])
        dataset.createDimension("location", covers.shape[0])
        dataset.createVariable("species", "i4", ("species",))[:] = species
        dataset.createVariable("location", "i4", ("location",))[:] = location
        variable = dataset.createVariable("covers", "f8", ("location", "species"))
        variable[:] = covers
        variable.units = units
        dataset.history = ("v0.2 initial cover with zero cover at sites whose benthic/geomorphic "
                           "habitable fraction is 0 (build_lizard_domain_v3.py)")

    # Cyclone mortality --------------------------------------------------------------
    hours = pd.read_csv(CYCLONE_CSV, dtype={"reef_id": str}).set_index("reef_id")
    hours.columns = hours.columns.astype(int)
    reef_ids = sites.GBRMPA_ID.astype(str)
    missing_reefs = sorted(set(reef_ids) - set(hours.index))
    if missing_reefs:
        raise ValueError(f"Cyclone hours missing for reefs: {missing_reefs}")
    site_hours = hours.reindex(columns=YEARS).loc[reef_ids].fillna(0.0).to_numpy(dtype=float)
    mortality = 1.0 - np.exp(B_PER_HOUR * site_hours)
    if not np.all((mortality >= 0.0) & (mortality <= 1.0)):
        raise ValueError("Cyclone mortality outside [0, 1]")
    (output / "cyclones").mkdir(exist_ok=True)
    cyclone_path = output / "cyclones" / "cyclone_mortality.csv"
    pd.DataFrame(mortality, columns=[str(y) for y in YEARS]).assign(
        site_id=sites.site_id.to_numpy()).set_index("site_id").to_csv(cyclone_path)

    # Package metadata --------------------------------------------------------------
    package_path = output / "datapackage.json"
    package = json.loads(package_path.read_text())
    package.update({
        "domain_data_version": "0.3.0",
        "description": "Lizard V2 sites with benthic/geomorphic site habitable fraction and "
                       "historical RRN cyclone mortality",
        "cyclone_mortality": "cyclones/cyclone_mortality.csv",
    })
    package_path.write_text(json.dumps(package, indent=2) + "\n")
    provenance_path = output / "provenance.json"
    provenance = json.loads(provenance_path.read_text())
    provenance["v0_3"] = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "builder": str(Path(__file__).relative_to(REPO_ROOT)),
        "builder_sha256": sha256_file(Path(__file__)),
        "source_package": str(SOURCE),
        "habitat_csv": str(HABITAT_CSV), "habitat_csv_sha256": sha256_file(HABITAT_CSV),
        "habitat_column": "k_geomorphic_infill", "k_floor": K_FLOOR,
        "k_zero_sites": int(zero_sites.sum()),
        "k_mean": float(sites.k.mean()),
        "cyclone_hours_csv": str(CYCLONE_CSV), "cyclone_hours_sha256": sha256_file(CYCLONE_CSV),
        "cyclone_mortality": f"1 - exp({B_PER_HOUR} * hours of >=4 m waves), all groups, "
                             "summer assigned to its ending calendar year",
        "cyclone_years": [YEARS[0], YEARS[-1]],
        "access": "derived from INTERNAL ONLY / licensed sources; do not publish or commit",
        "outputs": {name: sha256_file(output / name) for name in (
            "spatial/lizard_cluster.gpkg", "initial_cover.nc", "cyclones/cyclone_mortality.csv",
            "datapackage.json")},
    }
    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n")
    print(f"v0.3 written: mean k {sites.k.mean():.3f}; zero-k sites {int(zero_sites.sum())}; "
          f"max cyclone mortality {mortality.max():.3f}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    main(parser.parse_args().output)
