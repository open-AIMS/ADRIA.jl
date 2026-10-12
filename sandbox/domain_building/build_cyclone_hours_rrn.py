"""Extract historical cyclone wave exposure for the Lizard V2 domain reefs from
the licensed GBRMPA Resilient Reefs Network workbook (outputs stay in the
git-ignored sandbox/data tree).

Variable: ``cyc_maxHrs4mw`` = maximum hours of exposure to >=4 m significant
waves during any one cyclone in an austral summer (1 Nov to 30 Apr), per
GBRMPA reef ID, summers 1984/85 to 2024/25. Each summer is assigned to the
calendar year in which it ends (e.g. 199798 -> 1998), matching the model's
annual disturbance step. Missing reef-summers are written as NaN and treated
as zero exposure downstream; coverage is recorded.

Requires openpyxl (e.g. the base miniconda Python):

    python sandbox/domain_building/build_cyclone_hours_rrn.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import geopandas as gpd
import openpyxl
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_WORKBOOK = Path(
    r"C:\Users\smatthew\Documents\GitKraken\COTS-Coral-OutcomesV2\data\raw\GBRMPA_RRN_2025_AIMS.xlsx")
DEFAULT_SITES = REPO_ROOT / "sandbox/data/Lizard_Historical_v0.2/spatial/lizard_cluster.gpkg"
DEFAULT_OUTPUT = REPO_ROOT / "sandbox/data/Lizard_Historical_v0.2/derived/cyclone_hours4m_rrn2025.csv"
SHEET = "4.IndividualPressureData"
VARIABLE = "cyc_maxHrs4mw"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main(workbook: Path, sites_path: Path, output: Path) -> None:
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite {output}")
    reefs = sorted(gpd.read_file(sites_path).GBRMPA_ID.astype(str).unique())
    book = openpyxl.load_workbook(workbook, read_only=True)
    header = None
    rows = []
    for row in book[SHEET].iter_rows(values_only=True):
        if header is None:
            header = row
            if tuple(header[:5]) != ("LABEL_ID", "Type", "Year", "variable", "value"):
                raise ValueError(f"Unexpected {SHEET} header: {header}")
            continue
        if row[1] == VARIABLE and str(row[0]) in reefs:
            rows.append((str(row[0]), int(row[2]), row[4]))
    data = pd.DataFrame(rows, columns=["reef_id", "summer", "hours_4m"])
    if data.empty:
        raise ValueError("No cyclone rows matched the domain reefs")
    if data.duplicated(["reef_id", "summer"]).any():
        raise ValueError("Duplicate reef-summer cyclone rows")
    data["year"] = data.summer // 100 + 1  # summer YYYYyy ends in calendar year YYYY+1
    if not (data.summer % 100 == data.year % 100).all():
        raise ValueError("Unexpected summer coding")
    data["hours_4m"] = pd.to_numeric(data.hours_4m, errors="coerce")
    if (data.hours_4m < 0).any():
        raise ValueError("Negative cyclone hours")
    years = range(int(data.year.min()), int(data.year.max()) + 1)
    table = (data.pivot(index="reef_id", columns="year", values="hours_4m")
             .reindex(index=reefs, columns=list(years)))
    output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output, index_label="reef_id")
    matched = sorted(set(data.reef_id))
    provenance = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": str(Path(__file__).relative_to(REPO_ROOT)),
        "script_sha256": sha256_file(Path(__file__)),
        "workbook": str(workbook), "workbook_sha256": sha256_file(workbook),
        "sheet": SHEET, "variable": VARIABLE,
        "units": "maximum hours of >=4 m significant waves in any one cyclone, Nov-Apr summer",
        "year_convention": "summer YYYY/YY+1 assigned to calendar year YY+1",
        "sites": str(sites_path), "sites_sha256": sha256_file(sites_path),
        "domain_reefs": len(reefs), "matched_reefs": len(matched),
        "unmatched_reefs": sorted(set(reefs) - set(matched)),
        "missing_reef_years": int(table.isna().sum().sum()),
        "years": [min(years), max(years)],
        "access": "derived from licensed GBRMPA RRN data; do not publish or commit outputs",
        "output": str(output), "output_sha256": sha256_file(output),
    }
    output.with_suffix(".provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(f"matched {len(matched)}/{len(reefs)} reefs; years {min(years)}-{max(years)}; "
          f"missing reef-years {provenance['missing_reef_years']}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workbook", type=Path, default=DEFAULT_WORKBOOK)
    parser.add_argument("--sites", type=Path, default=DEFAULT_SITES)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    main(args.workbook, args.sites, args.output)
