"""Screen SALAD tracks against EOTR manta segments by reef, date, and midpoint."""

from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path

import numpy as np

from audit_cots_conversion_overlap import (
    EOTR, LABEL_TO_REEF, REEF_LABELS, SALAD, date, number, records, sha256,
)


def midpoint(lat0, lon0, lat1, lon1):
    return ((lat0 + lat1) / 2, (lon0 + lon1) / 2)


def distance_km(lat0, lon0, lat1, lon1):
    a0, a1 = np.radians([lat0, lat1])
    dlat = a1 - a0
    dlon = np.radians(lon1 - lon0)
    h = np.sin(dlat / 2) ** 2 + np.cos(a0) * np.cos(a1) * np.sin(dlon / 2) ** 2
    return float(6371.0 * 2 * np.arcsin(min(1.0, np.sqrt(h))))


def main(output: Path, max_days: int, max_km: float) -> None:
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite {output}")
    if max_days < 0 or max_km <= 0:
        raise ValueError("Invalid matching bounds")
    output.mkdir(parents=True)
    manta = defaultdict(list)
    for row in records(EOTR, "Manta"):
        label = str(row.get("ReefLabel") or "").strip()
        when = date(row.get("SurveyTime"))
        count = number(row.get("CrownOfThornsStarfishCount"))
        lat0 = number(row.get("StartLatitude"))
        lon0 = number(row.get("StartLongitude"))
        lat1 = number(row.get("EndLatitude"))
        lon1 = number(row.get("EndLongitude"))
        if label not in LABEL_TO_REEF or None in (when, count, lat0, lon0):
            continue
        lat, lon = midpoint(lat0, lon0, lat1 or lat0, lon1 or lon0)
        manta[LABEL_TO_REEF[label]].append({
            "date": when.date(), "count": count, "lat": lat, "lon": lon,
            "label": label, "tow_id": row.get("CrownOfThornsStarfishSurveillanceId"),
            "distance_m": number(row.get("TowDistance")),
        })
    rows = []
    no_coord = 0
    # The SALAD sheet repeats its coordinate headers, so use positional columns.
    import openpyxl
    workbook = openpyxl.load_workbook(SALAD, read_only=True, data_only=True)
    try:
        values_iter = iter(workbook["Tracks"].values)
        next(values_iter)
        for values in values_iter:
            reef = str(values[6] or "").strip()
            when = date(values[2])
            count = number(values[13])
            area = number(values[22])
            if reef not in REEF_LABELS or when is None or count is None or not area:
                continue
            lat0, lon0, lat1, lon1 = (number(value) for value in values[17:21])
            if None in (lat0, lon0):
                no_coord += 1
                continue
            lat, lon = midpoint(lat0, lon0, lat1 or lat0, lon1 or lon0)
            matches = []
            for tow in manta[reef]:
                day_gap = abs((when.date() - tow["date"]).days)
                if day_gap > max_days:
                    continue
                km = distance_km(lat, lon, tow["lat"], tow["lon"])
                if km <= max_km:
                    matches.append((tow, day_gap, km))
            rows.append({
                "salad_reef": reef, "salad_year": when.year,
                "salad_date": when.date().isoformat(),
                "salad_track": str(values[3]),
                "salad_cots_count_all_sizes": count,
                "salad_area_m2": area,
                "salad_all_sizes_per_ha": round(10000 * count / area, 6),
                "matched_manta_tows": len(matches),
                "matched_manta_cots_count": sum(tow["count"] for tow, _, _ in matches),
                "matched_manta_cots_per_tow":
                    round(sum(tow["count"] for tow, _, _ in matches) / len(matches), 6)
                    if matches else "",
                "nearest_midpoint_km": min((km for _, _, km in matches), default=""),
                "nearest_day_gap": min((days for _, days, _ in matches), default=""),
                "matched_tow_ids": ";".join(str(tow["tow_id"]) for tow, _, _ in matches),
                "matched_reef_labels": ";".join(sorted({tow["label"] for tow, _, _ in matches})),
            })
    finally:
        workbook.close()
    path = output / "track_overlap.csv"
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    matched = [row for row in rows if row["matched_manta_tows"] > 0]
    metadata = {
        "status": "spatiotemporal_screen_only_no_conversion_fitted",
        "max_days": max_days, "max_midpoint_distance_km": max_km,
        "salad_tracks_with_coordinates": len(rows),
        "salad_tracks_without_coordinates": no_coord,
        "matched_salad_tracks": len(matched),
        "matched_tracks_with_positive_manta_count":
            sum(row["matched_manta_cots_count"] > 0 for row in matched),
        "distinct_matched_manta_tow_ids": len({identifier for row in matched for identifier
                                             in row["matched_tow_ids"].split(";") if identifier}),
        "limitations": [
            "Midpoint distance is only a screen, not tow-path/track-line intersection",
            "Multiple SALAD tracks may reuse the same manta tows; matched rows are not independent",
            "SALAD COTS includes all sizes; individual-to-track count reconciliation is separate",
            "EOTR manta observation method may not transfer to historical LTMP",
        ],
        "inputs": {"salad": sha256(SALAD), "eotr": sha256(EOTR),
                   "script": sha256(Path(__file__)),
                   "shared_reader_and_crosswalk":
                       sha256(Path(__file__).with_name("audit_cots_conversion_overlap.py"))},
        "output_sha256": sha256(path),
    }
    (output / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps({key: metadata[key] for key in (
        "salad_tracks_with_coordinates", "salad_tracks_without_coordinates",
        "matched_salad_tracks", "matched_tracks_with_positive_manta_count",
        "distinct_matched_manta_tow_ids")}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-days", type=int, default=30)
    parser.add_argument("--max-km", type=float, default=1.0)
    args = parser.parse_args()
    main(args.output, args.max_days, args.max_km)
