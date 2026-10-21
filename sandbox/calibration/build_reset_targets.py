"""Projection-start reset targets (crash_mechanism_protocol.md, 2026-10-21, M2; ADRIA switch
COTS_STATE_RESET_CSV).

For each model reef with a diver-anchored target (targets_reef_year.csv: recorded SALAD density or
cull-converted, COTS ha^-1), the latest target year <= 2026 gives `target_adults_ha`; every site of that
model reef (site_to_reef.csv reef_name before " (") is listed under the reef as its group. The reset is
applied by ADRIA after the table's `year` (2026), whatever the target year (2025 where 2026 is missing).

Run with base miniconda Python:
    python sandbox/calibration/build_reset_targets.py
Outputs (git-ignored; derived from licensed data): sandbox/data/derived/reset_targets_2026/
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "sandbox" / "data"
TARGETS = DATA / "derived" / "density_targets_2020s" / "targets_reef_year.csv"
SITES = DATA / "Lizard_Historical_v0.4" / "site_to_reef.csv"
OUT = DATA / "derived" / "reset_targets_2026"
RESET_YEAR = 2026


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    if OUT.exists():
        raise FileExistsError(f"Refusing to overwrite {OUT}")
    targets = pd.read_csv(TARGETS)
    targets = targets[targets.report_year <= RESET_YEAR].sort_values("report_year").groupby("target_reef").tail(1)
    sites = pd.read_csv(SITES, dtype={"site_id": str})
    sites["model_reef"] = sites.reef_name.str.split(" (", regex=False).str[0]
    rows, summary = [], []
    for t in targets.itertuples(index=False):
        s = sites[sites.model_reef == t.target_reef]
        if s.empty:
            raise ValueError(f"No domain sites for {t.target_reef}")
        rows += [{"site_id": sid, "group": t.target_reef, "year": RESET_YEAR,
                  "target_adults_ha": float(t.target_density_ha)} for sid in s.site_id]
        summary.append({"group": t.target_reef, "target_year": int(t.report_year), "target_adults_ha": float(t.target_density_ha),
                        "source": t.target_source, "sites": int(len(s))})
    OUT.mkdir(parents=True)
    pd.DataFrame(rows).to_csv(OUT / "reset_targets.csv", index=False)
    report = {"status": "derived_input_not_promoted",
              "protocol": "crash_mechanism_protocol.md, 2026-10-21, M2 projection-start reset",
              "reset_after_year": RESET_YEAR, "groups": summary,
              "units": "target adult density COTS ha^-1 of site area (area-weighted over the group's sites)",
              "inputs": {"targets": sha256(TARGETS), "site_to_reef": sha256(SITES)},
              "script_sha256": sha256(Path(__file__)), "created": datetime.now().isoformat(timespec="seconds"),
              "outputs": {"reset_targets.csv": sha256(OUT / "reset_targets.csv")}}
    (OUT / "metadata.json").write_text(json.dumps(report, indent=2) + "\n")
    print(pd.DataFrame(summary).round(2).to_string(index=False))


if __name__ == "__main__":
    main()
