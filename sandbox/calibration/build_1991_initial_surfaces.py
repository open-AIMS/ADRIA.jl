"""Build an exploratory, observation-based 1991 Lizard initial surface.

Raw LTMP files are never modified. COTS counts/tows are pooled over 1989-1991;
9 m manta hard-coral medians are averaged over the same observed years. The two
quantities are interpolated independently using effort-aware and ordinary IDW,
respectively. The 1991 outbreak prediction is recorded, not treated as CPUE.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "sandbox" / "data"
DOMAIN = DATA / "Lizard_Historical_v0.2"
GLOBAL_REEFS = DATA / "OwenHydro" / "connectivityMatsYearlyGlobal" / "reefShapefiles" / "gbrShapeLL.shp"
COTS = DATA / "reef_cots.csv"
CORAL = DATA / "reef_manta.csv"
PREDICTIONS = DATA / "gbrPredsAdj_20262408.csv"
SITES = DOMAIN / "spatial" / "lizard_cluster.gpkg"
SITE_MAP = DOMAIN / "site_to_reef.csv"
YEARS = (1989, 1990, 1991)
RADIUS_KM = 111.0
MIN_DONORS = 3
DISTANCE_FLOOR_KM = 1.0
EXPLICIT_REEF_IDS = {
    # Audited LTMP name variants in the 1989-1991 Cape York records.
    "creech reef north": "13-118a",
    "log reef south": "12-107",
    "martin reef": "14-123",  # global shapefile also has a 19-075 Martin Reef
    "ribbon reef no.6": "15-032",
    "sand bank no.1": "14-045",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def normalized_name(name: str) -> str:
    return re.sub(r"\s+\([0-9]{2}-[0-9a-z]+\)$", "", str(name).strip(), flags=re.I).casefold()


def crosswalk(global_reefs: pd.DataFrame, survey_names: set[str]) -> pd.DataFrame:
    by_name: dict[str, list[int]] = {}
    by_id: dict[str, list[int]] = {}
    for idx, row in global_reefs.iterrows():
        by_name.setdefault(normalized_name(row.reefName), []).append(idx)
        by_id.setdefault(str(row.LABEL_ID).casefold(), []).append(idx)
    entries = []
    for name in sorted(survey_names):
        key = normalized_name(name)
        if key == "lizard isles":
            matches = [i for i, row in global_reefs.iterrows()
                       if str(row.reefName).startswith("Lizard Island Reef (")]
        elif key == "north direction island":
            matches = by_name.get("north direction reef", [])
        elif key in EXPLICIT_REEF_IDS:
            matches = by_id.get(EXPLICIT_REEF_IDS[key], [])
        else:
            id_match = re.fullmatch(r"reef\s+([0-9]{2}-[0-9a-z]+)", key)
            matches = by_id.get(id_match.group(1), []) if id_match else by_name.get(key, [])
        if not matches or (len(matches) != 1 and key != "lizard isles"):
            entries.append((name, "unmatched" if not matches else "ambiguous", len(matches),
                            np.nan, np.nan, ""))
            continue
        matched = global_reefs.loc[matches]
        weight = matched.Shape_Area.to_numpy(dtype=float)
        if not np.all(np.isfinite(weight)) or not np.all(weight > 0):
            raise ValueError(f"Invalid global reef area for {name}")
        lon = np.average(matched.xCentroid.to_numpy(dtype=float), weights=weight)
        lat = np.average(matched.yCentroid.to_numpy(dtype=float), weights=weight)
        entries.append((name, "matched", len(matches), lon, lat,
                        ";".join(matched.LABEL_ID.astype(str))))
    return pd.DataFrame(entries, columns=["survey_reef", "status", "global_reef_count",
                                           "longitude", "latitude", "global_reef_ids"])


def haversine_km(lon: np.ndarray, lat: np.ndarray,
                 donor_lon: np.ndarray, donor_lat: np.ndarray) -> np.ndarray:
    x1 = np.deg2rad(lon)[:, None]
    y1 = np.deg2rad(lat)[:, None]
    x2 = np.deg2rad(donor_lon)[None, :]
    y2 = np.deg2rad(donor_lat)[None, :]
    a = np.sin((y2 - y1) / 2) ** 2 + np.cos(y1) * np.cos(y2) * np.sin((x2 - x1) / 2) ** 2
    return 6371.0088 * 2 * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def idw(distances_km: np.ndarray, values: np.ndarray, precision: np.ndarray,
        radius_km: float = RADIUS_KM, min_donors: int = MIN_DONORS
        ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if distances_km.shape[1] != len(values) or len(values) != len(precision):
        raise ValueError("IDW donor dimensions differ")
    if np.any(~np.isfinite(values)) or np.any(~np.isfinite(precision)) or np.any(precision <= 0):
        raise ValueError("IDW donors must have finite values and positive precision")
    eligible = distances_km <= radius_km
    count = eligible.sum(axis=1)
    weights = np.where(eligible, precision[None, :] /
                       np.maximum(distances_km, DISTANCE_FLOOR_KM) ** 2, 0.0)
    totals = weights.sum(axis=1)
    estimate = np.divide(weights @ values, totals, out=np.full(len(count), np.nan), where=totals > 0)
    estimate[count < min_donors] = np.nan
    furthest = np.max(np.where(eligible, distances_km, 0.0), axis=1)
    furthest[count < min_donors] = np.nan
    return estimate, count, furthest


def build(output: Path) -> None:
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite {output}")
    if output.resolve() == ROOT or ROOT not in output.resolve().parents:
        raise ValueError("Output must be a new directory within the repository")
    reef_shape = gpd.read_file(GLOBAL_REEFS)[["LABEL_ID", "reefName", "Shape_Area",
                                              "xCentroid", "yCentroid"]]
    if reef_shape.LABEL_ID.duplicated().any():
        raise ValueError("Global reef IDs are not unique")
    cots = pd.read_csv(COTS)
    cots = cots[cots.year.isin(YEARS)].copy()
    if (cots.tows <= 0).any() or (cots.cots < 0).any():
        raise ValueError("Invalid COTS count or tow effort")
    pooled = cots.groupby("reef_name", as_index=False)[["cots", "tows"]].sum()
    pooled["cots_per_tow"] = pooled.cots / pooled.tows
    coral = pd.read_csv(CORAL)
    coral = coral[(coral.data_type == "manta") & (coral.domain_category == "reef") &
                  (coral.variable == "HC") & (coral.purpose == "MANTA") &
                  (coral.project_code == "LTMP") & (coral.depth == 9) &
                  coral.report_year.isin(YEARS)].copy()
    if coral.duplicated(["reef_name", "report_year"]).any():
        raise ValueError("Duplicate manta-cover reef/year records")
    if coral["median"].isna().any() or not coral["median"].between(0, 1).all():
        raise ValueError("Invalid manta-cover median")
    coral_pooled = coral.groupby("reef_name", as_index=False).agg(
        coral_manta_9m_fraction=("median", "mean"), coral_years=("report_year", "nunique"))
    xref = crosswalk(reef_shape, set(pooled.reef_name) | set(coral_pooled.reef_name))
    cots_donors = pooled.merge(xref[xref.status == "matched"], left_on="reef_name",
                               right_on="survey_reef", validate="one_to_one")
    coral_donors = coral_pooled.merge(xref[xref.status == "matched"], left_on="reef_name",
                                     right_on="survey_reef", validate="one_to_one")
    sites = gpd.read_file(SITES)[["site_id", "xCentroid", "yCentroid"]]
    mapping = pd.read_csv(SITE_MAP, dtype={"UNIQUE_ID": str, "reef_id": str})
    sites = sites.merge(mapping[["site_id", "reef_id", "area_m2"]], on="site_id",
                        validate="one_to_one")
    if len(sites) != 2895 or sites.site_id.duplicated().any():
        raise ValueError("Unexpected Lizard site set")
    predictions = pd.read_csv(PREDICTIONS)
    predictions = predictions[predictions.year == 1991].copy()
    if predictions.reefName.duplicated().any() or not predictions.outbrProb.between(0, 1).all():
        raise ValueError("Invalid 1991 outbreak-probability grid")
    predicted = reef_shape[["LABEL_ID", "reefName"]].merge(
        predictions[["reefName", "outbrProb"]], on="reefName", validate="one_to_one")
    if len(predicted) != len(reef_shape):
        raise ValueError("1991 prediction does not cover every global reef")
    sites = sites.merge(predicted[["LABEL_ID", "outbrProb"]], left_on="reef_id",
                        right_on="LABEL_ID", validate="many_to_one")
    if sites.outbrProb.isna().any():
        raise ValueError("1991 outbreak probability missing at a site")
    lon = sites.xCentroid.to_numpy(dtype=float)
    lat = sites.yCentroid.to_numpy(dtype=float)
    cots_dist = haversine_km(lon, lat, cots_donors.longitude.to_numpy(),
                             cots_donors.latitude.to_numpy())
    coral_dist = haversine_km(lon, lat, coral_donors.longitude.to_numpy(),
                              coral_donors.latitude.to_numpy())
    cots_cpue, cots_count, cots_radius = idw(
        cots_dist, cots_donors.cots_per_tow.to_numpy(dtype=float),
        cots_donors.tows.to_numpy(dtype=float))
    coral_cover, coral_count, coral_radius = idw(
        coral_dist, coral_donors.coral_manta_9m_fraction.to_numpy(dtype=float),
        np.ones(len(coral_donors)))
    def leave_one_reef_out(donors: pd.DataFrame, column: str,
                           precision: np.ndarray) -> pd.DataFrame:
        distances = haversine_km(donors.longitude.to_numpy(dtype=float),
                                 donors.latitude.to_numpy(dtype=float),
                                 donors.longitude.to_numpy(dtype=float),
                                 donors.latitude.to_numpy(dtype=float))
        np.fill_diagonal(distances, np.inf)
        observed = donors[column].to_numpy(dtype=float)
        predicted, count, _ = idw(distances, observed, precision)
        return pd.DataFrame({"survey_reef": donors.reef_name,
                             "observed": observed, "predicted": predicted,
                             "donor_reefs": count, "residual": observed - predicted})

    cots_loocv = leave_one_reef_out(cots_donors, "cots_per_tow",
                                    cots_donors.tows.to_numpy(dtype=float))
    coral_loocv = leave_one_reef_out(coral_donors, "coral_manta_9m_fraction",
                                     np.ones(len(coral_donors)))
    surface = pd.DataFrame({
        "site_id": sites.site_id, "reef_id": sites.reef_id,
        "site_area_m2": sites.area_m2, "longitude": lon, "latitude": lat,
        "cots_cpue_1989_1991_idw": cots_cpue,
        "cots_donor_reefs": cots_count, "cots_max_donor_km": cots_radius,
        "coral_manta_9m_fraction_1989_1991_idw": coral_cover,
        "coral_donor_reefs": coral_count, "coral_max_donor_km": coral_radius,
        "outbreak_probability_1991": sites.outbrProb,
    })
    output.mkdir(parents=True)
    files = {"surface.csv": surface, "crosswalk.csv": xref,
             "cots_donors.csv": cots_donors, "coral_donors.csv": coral_donors,
             "cots_loocv.csv": cots_loocv, "coral_loocv.csv": coral_loocv}
    for name, table in files.items():
        table.to_csv(output / name, index=False)
    metadata = {
        "status": "experimental_not_promoted",
        "years": list(YEARS), "model_start_year": 1991,
        "cots_units": "COTS per manta tow; sum(COTS)/sum(tows) at each donor reef",
        "coral_units": "9 m manta hard-coral cover fraction; mean of available annual medians",
        "outbreak_probability_meaning": "P(COTS/tow > 0.22), recorded only; not multiplied into CPUE",
        "interpolation": {"method": "IDW", "distance": "haversine km",
                          "radius_km": RADIUS_KM, "min_distinct_donor_reefs": MIN_DONORS,
                          "distance_floor_km": DISTANCE_FLOOR_KM,
                          "cots_weight": "tows / max(distance_km, 1)^2",
                          "coral_weight": "1 / max(distance_km, 1)^2"},
        "explicit_reef_id_crosswalk": EXPLICIT_REEF_IDS,
        "coverage": {"site_count": len(surface), "cots_donors": len(cots_donors),
                     "coral_donors": len(coral_donors),
                     "cots_missing_sites": int(np.isnan(cots_cpue).sum()),
                     "coral_missing_sites": int(np.isnan(coral_cover).sum()),
                     "cots_loocv_available": int(cots_loocv.predicted.notna().sum()),
                     "coral_loocv_available": int(coral_loocv.predicted.notna().sum()),
                     "cots_loocv_mae_per_tow": float(cots_loocv.residual.abs().mean()),
                     "coral_loocv_mae_fraction": float(coral_loocv.residual.abs().mean()),
                     "unmatched_survey_names": xref[xref.status != "matched"].survey_reef.tolist()},
        "inputs": {str(path.relative_to(ROOT)).replace("\\", "/"): sha256(path)
                   for path in (COTS, CORAL, PREDICTIONS, SITES, SITE_MAP,
                                Path(__file__), *sorted(GLOBAL_REEFS.parent.glob("gbrShapeLL.*")))
                   if path.is_file()},
        "outputs": {name: sha256(output / name) for name in files},
    }
    (output / "metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    print(json.dumps(metadata["coverage"], indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    build(args.output)
