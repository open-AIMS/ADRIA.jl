"""Audit 2021-2026 SALAD, EOTR manta-count, and cull overlap; fit no conversion."""

from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime
import hashlib
import json
from pathlib import Path

import openpyxl

ROOT = Path(__file__).resolve().parents[2]
SALAD = ROOT / "sandbox/data/COTS SALAD MASTER 2021-2026.xlsx"
EOTR = ROOT / "sandbox/data/260529-COTS-Manta-Cull-RHIS-Lawrence-CSIRO.xlsx"

# Explicit crosswalk. Lizard/South Direction pool multiple mapped reef IDs;
# matching is reef-level only and does not establish co-located survey support.
REEF_LABELS = {
    "Lizard": {"14-116a", "14-116b", "14-116c", "14-116d"},
    "MacGillivray": {"14-114"},
    "North Direction": {"14-143"},
    "South Direction": {"14-147a", "14-147b"},
    "Martin": {"14-123"},
    "McSweeney": {"11-016"},
}
LABEL_TO_REEF = {label: reef for reef, labels in REEF_LABELS.items() for label in labels}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def number(value) -> float | None:
    try:
        return float(value) if value not in (None, "") else None
    except (ValueError, TypeError):
        return None


def date(value) -> datetime | None:
    if isinstance(value, datetime):
        return value
    if value in (None, ""):
        return None
    try:
        return datetime.fromisoformat(str(value).replace("Z", "+00:00")).replace(tzinfo=None)
    except ValueError:
        return None


def records(path: Path, sheet: str):
    workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
    try:
        rows = workbook[sheet].values
        header = next(rows)
        for values in rows:
            yield dict(zip(header, values))
    finally:
        workbook.close()


def main(output: Path) -> None:
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite {output}")
    output.mkdir(parents=True)
    salad = defaultdict(lambda: {"tracks": 0, "count": 0.0, "area_m2": 0.0,
                                 "dates": set(), "density_reported": []})
    excluded = defaultdict(int)
    for row in records(SALAD, "Tracks"):
        reef = str(row.get("Reef") or "").strip()
        when = date(row.get("Date"))
        if reef not in REEF_LABELS or when is None:
            excluded["salad_other_reef_or_date"] += 1
            continue
        count = number(row.get("No. COTS"))
        area = number(row.get("AREA(5MWIDE)"))
        reported_density = number(row.get("Density(Ha)"))
        if count is None or area is None or area <= 0:
            excluded["salad_missing_count_or_area"] += 1
            continue
        group = salad[(reef, when.year)]
        group["tracks"] += 1
        group["count"] += count
        group["area_m2"] += area
        group["dates"].add(when.date())
        if reported_density is not None:
            group["density_reported"].append((reported_density, 10000 * count / area))

    manta = defaultdict(lambda: {"tows": 0, "count": 0.0, "distance_m": 0.0,
                                 "dates": set(), "labels": set()})
    for row in records(EOTR, "Manta"):
        label = str(row.get("ReefLabel") or "").strip()
        when = date(row.get("SurveyTime"))
        if label not in LABEL_TO_REEF or when is None:
            excluded["manta_other_reef_or_date"] += 1
            continue
        count = number(row.get("CrownOfThornsStarfishCount"))
        if count is None:
            excluded["manta_missing_count"] += 1
            continue
        group = manta[(LABEL_TO_REEF[label], when.year)]
        group["tows"] += 1
        group["count"] += count
        group["distance_m"] += number(row.get("TowDistance")) or 0.0
        group["dates"].add(when.date())
        group["labels"].add(label)

    cull = defaultdict(lambda: {"dives": 0, "count": 0.0, "adult25_count": 0.0,
                                "bottom_minutes": 0.0, "dates": set()})
    for row in records(EOTR, "Cull"):
        label = str(row.get("ReefLabel") or "").strip()
        when = date(row.get("SurveyDate"))
        if label not in LABEL_TO_REEF or when is None:
            excluded["cull_other_reef_or_date"] += 1
            continue
        cohorts = [number(row.get(f"Cohort{i}")) for i in range(1, 5)]
        minutes = number(row.get("Bottomtime"))
        if any(v is None for v in cohorts) or minutes is None or minutes <= 0:
            excluded["cull_missing_cohort_or_time"] += 1
            continue
        group = cull[(LABEL_TO_REEF[label], when.year)]
        group["dives"] += 1
        group["count"] += sum(cohorts)
        group["adult25_count"] += cohorts[2] + cohorts[3]
        group["bottom_minutes"] += minutes
        group["dates"].add(when.date())

    rows = []
    for key, group in sorted(salad.items()):
        reef, year = key
        tow = manta.get(key)
        dive = cull.get(key)
        if tow is None:
            excluded["salad_reef_year_without_manta"] += 1
        salad_dates = group["dates"]
        manta_dates = tow["dates"] if tow else set()
        minimum_days = min((abs((a - b).days) for a in salad_dates for b in manta_dates),
                           default=None)
        checks = group["density_reported"]
        max_density_discrepancy = max((abs(a - b) for a, b in checks), default=None)
        rows.append({
            "salad_reef": reef, "year": year,
            "manta_reef_labels": ";".join(sorted(tow["labels"])) if tow else "",
            "salad_tracks": group["tracks"],
            "salad_cots_count_all_sizes": group["count"],
            "salad_survey_area_m2": round(group["area_m2"], 3),
            "salad_cots_all_sizes_per_ha": round(10000 * group["count"] / group["area_m2"], 6),
            "salad_reported_density_max_abs_discrepancy": max_density_discrepancy,
            "manta_tows": tow["tows"] if tow else 0,
            "manta_cots_count": tow["count"] if tow else 0,
            "manta_cots_per_tow": round(tow["count"] / tow["tows"], 6) if tow else None,
            "manta_mean_tow_distance_m": round(tow["distance_m"] / tow["tows"], 2) if tow else None,
            "minimum_salad_manta_survey_day_gap": minimum_days,
            "cull_dives": dive["dives"] if dive else 0,
            "cull_all_cohorts_per_bottom_minute":
                round(dive["count"] / dive["bottom_minutes"], 6) if dive else None,
            "cull_25cm_plus_per_bottom_minute":
                round(dive["adult25_count"] / dive["bottom_minutes"], 6) if dive else None,
        })

    import csv
    csv_path = output / "reef_year_overlap.csv"
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    matched = [row for row in rows if row["manta_tows"] > 0]
    metadata = {
        "status": "overlap_audit_only_no_conversion_fitted",
        "scope": "2021-2026 reef-year summaries for six explicit named reefs; EOTR manta/cull, not LTMP",
        "limitations": [
            "SALAD density includes all recorded COTS sizes; model observation state is adults/ha",
            "reef-year overlap does not prove same habitat, site, date, observer or detection probability",
            "Lizard and South Direction each pool multiple hydrodynamic reef IDs",
            "culling is targeted and catch per bottom minute is not unbiased density",
            "tow distance varies and visible-area/detection factors are unknown",
        ],
        "matched_reef_years": len(matched),
        "matched_within_30_days": sum(row["minimum_salad_manta_survey_day_gap"] <= 30
                                      for row in matched),
        "matched_with_cull": sum(row["cull_dives"] > 0 for row in matched),
        "excluded": dict(excluded),
        "reef_crosswalk": {reef: sorted(labels) for reef, labels in REEF_LABELS.items()},
        "inputs_sha256": {str(path.relative_to(ROOT)): sha256(path) for path in (SALAD, EOTR)},
        "script_sha256": sha256(Path(__file__)),
        "output_sha256": sha256(csv_path),
    }
    with (output / "metadata.json").open("w") as stream:
        json.dump(metadata, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps({k: metadata[k] for k in
                      ("matched_reef_years", "matched_within_30_days", "matched_with_cull",
                       "excluded")}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    main(parser.parse_args().output)
