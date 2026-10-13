"""Build a site x year table of observed adult COTS density for forced-predation
coral calibration (crash_mechanism_protocol.md, 2026-10-12).

For each calendar year 1991-2024, LTMP manta COTS counts per tow at each surveyed
reef (sum(COTS)/sum(tows) within the year) are interpolated to Lizard domain sites
with the same effort-aware IDW, reef crosswalk, radius and donor rule as
build_1991_initial_surfaces.py. Site-years with too few donors are filled by linear
interpolation in time at that site (ends held constant). Adults are CPUE / q with the
approved observation unit q = 0.015 COTS/tow per adult COTS/ha. Raw files are never
modified; the output directory must not exist.

    python sandbox/calibration/build_observed_cots_forcing.py --output <new dir>
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_1991_initial_surfaces import (  # noqa: E402
    COTS, GLOBAL_REEFS, ROOT, SITE_MAP, SITES, crosswalk, haversine_km, idw, sha256,
    MIN_DONORS, RADIUS_KM, DISTANCE_FLOOR_KM, EXPLICIT_REEF_IDS)

YEARS = list(range(1991, 2025))
Q = 0.015  # COTS/tow per adult COTS/ha (approved convention)


def build(output: Path) -> None:
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite {output}")
    reef_shape = gpd.read_file(GLOBAL_REEFS)[["LABEL_ID", "reefName", "Shape_Area",
                                              "xCentroid", "yCentroid"]]
    cots = pd.read_csv(COTS)
    cots = cots[cots.year.isin(YEARS)].copy()
    if (cots.tows <= 0).any() or (cots.cots < 0).any():
        raise ValueError("Invalid COTS count or tow effort")
    yearly = cots.groupby(["reef_name", "year"], as_index=False)[["cots", "tows"]].sum()
    yearly["cots_per_tow"] = yearly.cots / yearly.tows
    xref = crosswalk(reef_shape, set(yearly.reef_name))
    donors = yearly.merge(xref[xref.status == "matched"], left_on="reef_name",
                          right_on="survey_reef", validate="many_to_one")
    sites = gpd.read_file(SITES)[["site_id", "xCentroid", "yCentroid"]]
    mapping = pd.read_csv(SITE_MAP, dtype={"UNIQUE_ID": str, "reef_id": str})
    sites = sites.merge(mapping[["site_id", "reef_id"]], on="site_id", validate="one_to_one")
    if len(sites) != 2895 or sites.site_id.duplicated().any():
        raise ValueError("Unexpected Lizard site set")
    lon = sites.xCentroid.to_numpy(dtype=float)
    lat = sites.yCentroid.to_numpy(dtype=float)

    cpue = np.full((len(sites), len(YEARS)), np.nan)
    coverage = []
    for j, year in enumerate(YEARS):
        d = donors[donors.year == year]
        if len(d) == 0:
            coverage.append({"year": year, "donor_reefs": 0, "sites_estimated": 0})
            continue
        dist = haversine_km(lon, lat, d.longitude.to_numpy(), d.latitude.to_numpy())
        estimate, count, _ = idw(dist, d.cots_per_tow.to_numpy(dtype=float),
                                 d.tows.to_numpy(dtype=float))
        cpue[:, j] = estimate
        coverage.append({"year": year, "donor_reefs": int(len(d)),
                         "sites_estimated": int(np.isfinite(estimate).sum())})
    interpolated = pd.DataFrame(cpue.T, index=YEARS).interpolate(
        method="index", limit_direction="both").to_numpy().T
    filled = np.isnan(cpue) & np.isfinite(interpolated)
    if np.isnan(interpolated).any():
        raise ValueError("Some sites have no COTS estimate in any year")
    adults = interpolated / Q

    output.mkdir(parents=True)
    table = pd.DataFrame(adults, columns=[str(y) for y in YEARS])
    table.insert(0, "site_id", sites.site_id.astype(str).to_numpy())
    table.to_csv(output / "observed_adults_ha.csv", index=False)
    cpue_table = pd.DataFrame(interpolated, columns=[str(y) for y in YEARS])
    cpue_table.insert(0, "site_id", sites.site_id.astype(str).to_numpy())
    cpue_table.to_csv(output / "observed_cpue_idw.csv", index=False)
    pd.DataFrame(coverage).to_csv(output / "coverage_by_year.csv", index=False)
    metadata = {
        "status": "experimental_forcing_not_promoted",
        "purpose": "COTS predation forcing for coral-side calibration (forced grazing)",
        "years": [YEARS[0], YEARS[-1]],
        "units": "adult COTS per ha of site area = IDW COTS/tow / q",
        "q_cots_per_tow_per_cots_ha": Q,
        "interpolation": {"method": "IDW per calendar year", "distance": "haversine km",
                          "radius_km": RADIUS_KM, "min_distinct_donor_reefs": MIN_DONORS,
                          "distance_floor_km": DISTANCE_FLOOR_KM,
                          "weight": "tows / max(distance_km, 1)^2",
                          "temporal_fill": "linear in time per site; ends held constant",
                          "explicit_reef_id_crosswalk": EXPLICIT_REEF_IDS},
        "site_years_temporally_filled": int(filled.sum()),
        "site_years_total": int(filled.size),
        "unmatched_survey_names": xref[xref.status != "matched"].survey_reef.tolist(),
        "inputs": {str(path.relative_to(ROOT)).replace("\\", "/"): sha256(path)
                   for path in (COTS, SITES, SITE_MAP, Path(__file__),
                                Path(__file__).parent / "build_1991_initial_surfaces.py",
                                *sorted(GLOBAL_REEFS.parent.glob("gbrShapeLL.*")))
                   if path.is_file()},
        "outputs": {name: sha256(output / name) for name in (
            "observed_adults_ha.csv", "observed_cpue_idw.csv", "coverage_by_year.csv")},
    }
    (output / "metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: metadata[k] for k in ("site_years_temporally_filled", "site_years_total")}))
    print(pd.DataFrame(coverage).to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    build(parser.parse_args().output)
