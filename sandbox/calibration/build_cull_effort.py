"""Recent cull effort per domain reef for forward effort culling (crash_mechanism_protocol.md,
2026-10-21, M1; ADRIA switch COTS_CULL_EFFORT_CSV).

Effort = summed dive bottom time (sheet Cull, `Bottomtime`, minutes -> diver-hours), the same
definition as the cull CPUE behind the SALAD-cull density anchor (build_2020s_density_targets.py).
Reporting year: survey month >= September counts to the next year (LTMP convention). Recent effort =
mean diver-hours per year over reporting years 2023-2025 (2026 is partial: data end May 2026), per
GBRMPA reef id parsed from `ReefName` "(14-116b)". Reefs outside the domain are reported, not used.

Run with base miniconda Python (pandas, openpyxl):
    python sandbox/calibration/build_cull_effort.py
Outputs (git-ignored; licensed source): sandbox/data/derived/cull_effort_lizard_v1/
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "sandbox" / "data"
CSIRO = DATA / "260529-COTS-Manta-Cull-RHIS-Lawrence-CSIRO.xlsx"
SITES = DATA / "Lizard_Historical_v0.4" / "site_to_reef.csv"
OUT = DATA / "derived" / "cull_effort_lizard_v1"
YEARS = (2023, 2024, 2025)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    if OUT.exists():
        raise FileExistsError(f"Refusing to overwrite {OUT}")
    cull = pd.read_excel(CSIRO, sheet_name="Cull")
    t = pd.to_datetime(cull.SurveyDate)
    cull["report_year"] = t.dt.year + (t.dt.month >= 9).astype(int)
    cull["reef_id"] = cull.ReefName.str.extract(r"\((\d\d-\d{3}[a-z]?)\)\s*$")[0]
    cull["hours"] = cull.Bottomtime.fillna(0.0).astype(float) / 60.0
    cull["cots"] = cull[["Cohort1", "Cohort2", "Cohort3", "Cohort4"]].fillna(0).sum(axis=1)
    domain = set(pd.read_csv(SITES, dtype={"reef_id": str}).reef_id)
    by_year = (cull.groupby(["reef_id", "report_year"]).agg(hours=("hours", "sum"), cots=("cots", "sum"),
                                                            dives=("hours", "size")).reset_index())
    by_year["in_domain"] = by_year.reef_id.isin(domain)
    recent = by_year[by_year.report_year.isin(YEARS)]
    effort = (recent[recent.in_domain].groupby("reef_id").hours.sum() / len(YEARS)).rename("hours_per_year").reset_index()
    effort = effort[effort.hours_per_year > 0].sort_values("reef_id")
    OUT.mkdir(parents=True)
    effort.to_csv(OUT / "effort_recent_2023_2025.csv", index=False)
    by_year.to_csv(OUT / "effort_by_reef_year.csv", index=False)
    outside = recent[~recent.in_domain]
    report = {
        "status": "derived_input_not_promoted",
        "protocol": "crash_mechanism_protocol.md, 2026-10-21, M1 forward effort culling",
        "definition": "diver-hours = sum(Bottomtime minutes)/60 per reef and reporting year; recent = mean over reporting years 2023-2025",
        "units": "diver-hours per year (COTS_CULL_EFFORT_CSV hours_per_year)",
        "years": list(YEARS), "domain_reefs_with_effort": int(len(effort)),
        "domain_hours_per_year": float(effort.hours_per_year.sum()),
        "domain_hours_by_year": recent[recent.in_domain].groupby("report_year").hours.sum().round(1).to_dict(),
        "domain_cots_by_year": recent[recent.in_domain].groupby("report_year").cots.sum().astype(int).to_dict(),
        "outside_domain_hours_per_year": float(outside.hours.sum() / len(YEARS)),
        "licence": "CSIRO workbook internal/licensed: outputs git-ignored",
        "inputs": {"csiro_workbook": sha256(CSIRO), "site_to_reef": sha256(SITES)},
        "script_sha256": sha256(Path(__file__)), "created": datetime.now().isoformat(timespec="seconds"),
        "outputs": {n: sha256(OUT / n) for n in ("effort_recent_2023_2025.csv", "effort_by_reef_year.csv")},
    }
    (OUT / "metadata.json").write_text(json.dumps(report, indent=2, default=str) + "\n")
    print(effort.sort_values("hours_per_year", ascending=False).head(20).round(1).to_string(index=False))
    print(json.dumps({k: report[k] for k in ("domain_reefs_with_effort", "domain_hours_per_year", "domain_hours_by_year",
                                             "domain_cots_by_year", "outside_domain_hours_per_year")}, indent=1, default=str))


if __name__ == "__main__":
    main()
