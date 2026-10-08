"""Check SALAD COTS size support for an adult-density observation outcome."""

from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import hashlib
import json
from pathlib import Path

import openpyxl

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "sandbox/data/COTS SALAD MASTER 2021-2026.xlsx"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main(overlap: Path, output: Path) -> None:
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite {output}")
    output.mkdir(parents=True)
    with overlap.open(newline="") as stream:
        reef_year = list(csv.DictReader(stream))
    keys = {(row["salad_reef"], int(row["year"])) for row in reef_year}
    by_key = defaultdict(lambda: {"individuals": 0, "sized": 0,
                                  "at_least_250mm": 0, "at_least_260mm": 0,
                                  "scars": 0})
    workbook = openpyxl.load_workbook(SOURCE, read_only=True, data_only=True)
    try:
        rows = workbook["COTS recorded"].values
        header = next(rows)
        for values in rows:
            row = dict(zip(header, values))
            try:
                key = (str(row["Reef"]).strip(), int(row["Year"]))
            except (TypeError, ValueError):
                continue
            if key not in keys:
                continue
            kind = str(row["CoTS/ Scar"] or "").strip().lower()
            group = by_key[key]
            if kind == "scar":
                group["scars"] += 1
            elif kind == "cots":
                group["individuals"] += 1
                try:
                    size = float(row["Size(mm)"])
                except (ValueError, TypeError):
                    continue
                group["sized"] += 1
                group["at_least_250mm"] += size >= 250
                group["at_least_260mm"] += size >= 260
    finally:
        workbook.close()
    result = []
    for row in reef_year:
        key = (row["salad_reef"], int(row["year"]))
        group = by_key[key]
        track_count = float(row["salad_cots_count_all_sizes"])
        area = float(row["salad_survey_area_m2"])
        result.append({
            "salad_reef": key[0], "year": key[1],
            "salad_track_cots_count": track_count,
            "salad_individual_cots_records": group["individuals"],
            "count_reconciliation_difference": group["individuals"] - track_count,
            "sized_individuals": group["sized"],
            "at_least_250mm": group["at_least_250mm"],
            "at_least_260mm": group["at_least_260mm"],
            "adult250_per_ha_if_reconciled":
                round(10000 * group["at_least_250mm"] / area, 6)
                if abs(group["individuals"] - track_count) < 1e-9 else "",
            "adult260_per_ha_if_reconciled":
                round(10000 * group["at_least_260mm"] / area, 6)
                if abs(group["individuals"] - track_count) < 1e-9 else "",
            "scar_records": group["scars"],
            "manta_tows": row["manta_tows"],
            "manta_cots_per_tow": row["manta_cots_per_tow"],
            "minimum_salad_manta_survey_day_gap":
                row["minimum_salad_manta_survey_day_gap"],
        })
    output_csv = output / "salad_stage_sizes.csv"
    with output_csv.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(result[0]))
        writer.writeheader()
        writer.writerows(result)
    summary = {
        "status": "stage_size_audit_only_no_density_conversion_fitted",
        "rows": len(result),
        "exact_count_reconciliations": sum(row["count_reconciliation_difference"] == 0
                                            for row in result),
        "individuals": sum(row["salad_individual_cots_records"] for row in result),
        "sized_individuals": sum(row["sized_individuals"] for row in result),
        "at_least_250mm": sum(row["at_least_250mm"] for row in result),
        "at_least_260mm": sum(row["at_least_260mm"] for row in result),
        "note": "250/260 mm are sensitivity cutoffs; blank adult densities mean track and individual counts do not exactly reconcile. No historical 1991 age information is supplied.",
        "inputs": {"salad_workbook": sha256(SOURCE), "overlap": sha256(overlap),
                   "script": sha256(Path(__file__))},
        "output_sha256": sha256(output_csv),
    }
    (output / "metadata.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: summary[k] for k in ("rows", "exact_count_reconciliations",
                                             "individuals", "sized_individuals",
                                             "at_least_250mm", "at_least_260mm")}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--overlap", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    main(args.overlap, args.output)
