"""Observed COTS cull removals on the Lizard v0.3 domain (close-out C1, 2026-10-14).

Source: CSIRO workbook `260529-COTS-Manta-Cull-RHIS-Lawrence-CSIRO.xlsx`, sheet `Cull`
(one row per cull dive; Cohort1 <15 cm, Cohort2 15-<25 cm, Cohort3 25-40 cm, Cohort4 >40 cm;
sources AMPTO Legacy 2002-06, AMPTO 2012-17, EOTR 2018-). Data licence: outputs stay in the
git-ignored sandbox/data/derived; never publish or commit them.

Transformations (no modelling decisions; the cohort-to-stage mapping and the year basis are
applied by the model adapter):
1. keep dives whose ReefLabel is a domain reef (site_to_reef.csv reef_id); dives without a
   ReefLabel (inter-reefal) or on other reefs are counted and dropped;
2. map each cull site (CullSiteName centroid) to the nearest model site on the same reef
   (haversine on site centroids from spatial/lizard_cluster.gpkg); the distance is recorded;
3. aggregate counts by model site x calendar year x LTMP reporting year (Sep-Dec dives belong
   to the following reporting year, as for LTMP surveys);
4. densities: COTS culled per ha of model site area (site_to_reef.csv area_m2 / 1e4), the
   model's COTS state unit; reef tables divide by the summed model area of the reef.
5. v2 (lead's decisions, 2026-10-14, fifth round): adapter tables `removals_sigma<s>.csv`
   (site_id, report_year, phase, juveniles_ha, adults_ha) for the opt-in ADRIA switch
   COTS_CULL_REMOVALS_CSV. Each cull site's counts are fuzzy-matched to model sites on the
   same reef with weights area x exp(-d^2 / 2 sigma^2) over sites within 3 sigma (nearest
   site if none; sigma = 0 is nearest only). sigma = 210 m (primary; half the median
   cull-site spacing, 422 m) with sensitivities 0 and 420 m. Phase: Sep-Dec dives are
   `pre` (removed before spawning in the following reporting year), Jan-Aug dives `post`.
   Cohort 1 (<15 cm) -> juveniles; cohorts 2-4 (>=15 cm) -> adults.

    python build_cull_removals.py [--output <dir>]
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "sandbox" / "data"
WORKBOOK = DATA / "260529-COTS-Manta-Cull-RHIS-Lawrence-CSIRO.xlsx"
DOMAIN = DATA / "Lizard_Historical_v0.3"
SITES = DOMAIN / "site_to_reef.csv"
GPKG = DOMAIN / "spatial" / "lizard_cluster.gpkg"
DEFAULT_OUT = DATA / "derived" / "cots_cull_removals_lizard_v0_3_v2"
SIGMAS_M = (0.0, 210.0, 420.0)
COHORTS = ("Cohort1", "Cohort2", "Cohort3", "Cohort4")
COHORT_SIZES_CM = {"Cohort1": "<15", "Cohort2": "15-<25", "Cohort3": "25-40", "Cohort4": ">40"}
FAR_M = 1000.0


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def haversine_m(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = map(np.radians, (lat1, lon1, lat2, lon2))
    a = np.sin((lat2 - lat1) / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2
    return 2 * 6_371_000.0 * np.arcsin(np.sqrt(a))


def fuzzy_weights(cull_sites: pd.DataFrame, sites: pd.DataFrame, sigma: float) -> pd.DataFrame:
    """Rows (CullSiteName, site_id, weight); weights sum to 1 per cull site."""
    rows = []
    for name, s in cull_sites.iterrows():
        candidates = sites[sites.reef_id == s.reef_id]
        d = haversine_m(s.lat, s.lon, candidates.y_coord.to_numpy(), candidates.x_coord.to_numpy())
        near = d <= 3 * sigma if sigma > 0 else np.zeros(len(d), dtype=bool)
        if near.any():
            w = candidates.area_ha.to_numpy()[near] * np.exp(-d[near] ** 2 / (2 * sigma ** 2))
            ids = candidates.site_id.to_numpy()[near]
        else:
            w, ids = np.array([1.0]), np.array([candidates.site_id.iloc[int(np.argmin(d))]])
        for site_id, weight in zip(ids, w / w.sum()):
            rows.append({"CullSiteName": name, "site_id": site_id, "weight": float(weight)})
    return pd.DataFrame(rows)


def main(out: Path) -> None:
    if out.exists():
        raise FileExistsError(f"Refusing to overwrite {out}")
    sites = pd.read_csv(SITES)
    centroids = gpd.read_file(GPKG)[["site_id", "x_coord", "y_coord"]]
    sites = sites.merge(centroids, on="site_id", validate="one_to_one")
    if len(sites) != len(centroids) or sites.area_m2.le(0).any():
        raise ValueError("Site table and spatial layer disagree, or non-positive site area")
    sites["area_ha"] = sites.area_m2 / 1e4

    cull = pd.read_excel(WORKBOOK, sheet_name="Cull")
    n_all = len(cull)
    if cull[list(COHORTS)].lt(0).any().any() or cull.Bottomtime.lt(0).any():
        raise ValueError("Negative cull counts or bottom time")
    date = pd.to_datetime(cull.SurveyDate.astype(str).str[:10], format="%Y-%m-%d")
    cull = cull.assign(date=date, calendar_year=date.dt.year,
                       report_year=date.dt.year + (date.dt.month >= 9).astype(int))
    no_label = int(cull.ReefLabel.isna().sum())
    domain = cull[cull.ReefLabel.isin(set(sites.reef_id))].copy()
    if domain[["Latitude", "Longitude", "CullSiteName"]].isna().any().any():
        raise ValueError("Domain cull dives without a site or centroid")

    # Cull site -> nearest model site on the same reef.
    cull_sites = domain.groupby("CullSiteName").agg(reef_id=("ReefLabel", "first"), lat=("Latitude", "mean"),
                                                    lon=("Longitude", "mean"), n_reefs=("ReefLabel", "nunique"))
    if cull_sites.n_reefs.gt(1).any():
        raise ValueError("A cull site spans more than one reef label")
    mapped = []
    for name, s in cull_sites.iterrows():
        candidates = sites[sites.reef_id == s.reef_id]
        d = haversine_m(s.lat, s.lon, candidates.y_coord.to_numpy(), candidates.x_coord.to_numpy())
        i = int(np.argmin(d))
        mapped.append({"CullSiteName": name, "reef_id": s.reef_id, "cull_lat": s.lat, "cull_lon": s.lon,
                       "site_id": candidates.site_id.iloc[i], "distance_m": float(d[i])})
    site_map = pd.DataFrame(mapped)
    domain = domain.merge(site_map[["CullSiteName", "site_id"]], on="CullSiteName", validate="many_to_one")
    domain["culled"] = domain[list(COHORTS)].sum(axis=1)

    keys = ["site_id", "calendar_year", "report_year"]
    site_year = (domain.groupby(keys)
                 .agg(**{c: (c, "sum") for c in COHORTS}, culled=("culled", "sum"), dives=("culled", "size"),
                      bottomtime_min=("Bottomtime", "sum"),
                      sources=("Source", lambda v: ";".join(sorted(set(v)))))
                 .reset_index()
                 .merge(sites[["site_id", "reef_id", "reef_name", "area_ha"]], on="site_id", validate="many_to_one"))
    for c in (*COHORTS, "culled"):
        site_year[f"{c.lower()}_ha"] = site_year[c] / site_year.area_ha
    reef_area = sites.groupby("reef_id").area_ha.sum()
    reef_year = (domain.groupby(["ReefLabel", "calendar_year"])
                 .agg(**{c: (c, "sum") for c in COHORTS}, culled=("culled", "sum"), dives=("culled", "size"),
                      bottomtime_min=("Bottomtime", "sum"), cull_sites=("CullSiteName", "nunique"),
                      model_sites=("site_id", "nunique"), sources=("Source", lambda v: ";".join(sorted(set(v)))))
                 .reset_index().rename(columns={"ReefLabel": "reef_id"}))
    reef_year["reef_area_ha"] = reef_year.reef_id.map(reef_area)
    reef_year["culled_ha"] = reef_year.culled / reef_year.reef_area_ha

    # Mass balance: every domain dive's COTS reach exactly one site-year row.
    if int(site_year.culled.sum()) != int(domain.culled.sum()) or int(reef_year.culled.sum()) != int(domain.culled.sum()):
        raise ValueError("Cull totals not preserved")

    # Adapter tables per kernel width (lead's decisions 3a-3d).
    domain["phase"] = np.where(domain.date.dt.month >= 9, "pre", "post")
    domain["juveniles"] = domain.Cohort1
    domain["adults"] = domain[["Cohort2", "Cohort3", "Cohort4"]].sum(axis=1)
    area = sites.set_index("site_id").area_ha
    adapter, kernel_stats = {}, {}
    for sigma in SIGMAS_M:
        weights = fuzzy_weights(cull_sites, sites, sigma)
        spread = domain[["CullSiteName", "report_year", "phase", "juveniles", "adults"]].merge(
            weights, on="CullSiteName", validate="many_to_many")
        spread[["juveniles", "adults"]] = spread[["juveniles", "adults"]].mul(spread.weight, axis=0)
        table = spread.groupby(["site_id", "report_year", "phase"])[["juveniles", "adults"]].sum().reset_index()
        for col in ("juveniles", "adults"):
            if not np.isclose(table[col].sum(), domain[col].sum(), rtol=1e-12):
                raise ValueError(f"Fuzzy matching does not preserve {col} at sigma {sigma}")
            table[f"{col}_ha"] = table[col] / table.site_id.map(area)
        name = f"removals_sigma{int(sigma):03d}.csv"
        adapter[name] = table[["site_id", "report_year", "phase", "juveniles_ha", "adults_ha", "juveniles", "adults"]]
        kernel_stats[name] = {"sigma_m": sigma, "model_sites": int(table.site_id.nunique()),
                              "sites_per_cull_site_median": float(weights.groupby("CullSiteName").size().median()),
                              "max_site_year_density_ha": float((table.juveniles_ha + table.adults_ha).max())}

    out.mkdir(parents=True)
    files = {"cull_site_year.csv": site_year, "cull_reef_year.csv": reef_year, "cull_site_map.csv": site_map} | adapter
    for name, frame in files.items():
        frame.to_csv(out / name, index=False)
    meta = {
        "status": "derived_input_not_promoted",
        "licence": "CSIRO workbook under data licence: derived outputs are internal (git-ignored), never published",
        "source": {"workbook": str(WORKBOOK.relative_to(ROOT)), "sheet": "Cull", "sha256": sha256(WORKBOOK)},
        "domain": {"site_to_reef": str(SITES.relative_to(ROOT)), "site_to_reef_sha256": sha256(SITES),
                   "gpkg": str(GPKG.relative_to(ROOT)), "gpkg_sha256": sha256(GPKG)},
        "cohorts_cm": COHORT_SIZES_CM,
        "units": {"cohort*_ha, culled_ha": "COTS culled per ha of model site (or summed reef) area",
                  "bottomtime_min": "minutes"},
        "year_basis": "calendar_year = dive date year; report_year = LTMP reporting year (Sep-Dec -> year + 1)",
        "site_mapping": "cull_site_year / cull_site_map: nearest model site on the same reef (haversine); "
                        "removals_sigma*: Gaussian fuzzy match, weights area x exp(-d^2/2 sigma^2) within 3 sigma",
        "adapter_tables": kernel_stats,
        "stage_mapping": "juveniles = Cohort1 (<15 cm); adults = Cohort2-4 (>=15 cm)",
        "phase": "pre = Sep-Dec dives (before spawning, following reporting year); post = Jan-Aug dives (after spawning)",
        "counts": {"dives_all": n_all, "dives_without_reef_label": no_label, "dives_domain": int(len(domain)),
                   "domain_reefs_culled": int(domain.ReefLabel.nunique()), "cull_sites": int(len(site_map)),
                   "model_sites_culled": int(site_year.site_id.nunique()),
                   "cull_sites_over_1km_from_model_site": int((site_map.distance_m > FAR_M).sum()),
                   "culled_total": int(domain.culled.sum()),
                   "culled_by_source": domain.groupby("Source").culled.sum().astype(int).to_dict(),
                   "calendar_years": [int(domain.calendar_year.min()), int(domain.calendar_year.max())]},
        "dimensions": {name: list(frame.shape) for name, frame in files.items()},
        "script_sha256": sha256(Path(__file__)),
        "outputs": {name: sha256(out / name) for name in files},
    }
    (out / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps(meta["counts"], indent=1))
    print(json.dumps(kernel_stats, indent=1))
    print(reef_year.groupby("calendar_year")[["culled", "dives"]].sum().T.to_string())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    main(parser.parse_args().output.resolve())
