"""Plausible 1991 starting fast-coral share per calibration reef (crash_mechanism_protocol.md,
2026-10-14, partial calibration of the starting fast share).

Fast share = (Acropora + Pocilloporidae) / all hard coral on LTMP 9 m photo transects
(COMPOSITION), per reporting year, over ALL available years. The plausible range for a
pre-outbreak 1991 state keeps only years in which that reef's LTMP manta hard-coral cover was
at or above its own median (healthy, non-degraded states). Levels: min, q25, q50, q75, max of
that filtered range; reefs without transects take the mean of the three reefs at each level.

    python build_initial_share_range.py <run_dir>
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "sandbox" / "data"
OUT_DIR = DATA / "derived" / "initial_fast_share_range"
REEFS = {"Lizard Island Reef": "Lizard Isles", "MacGillivray Reef": "Macgillivray Reef",
         "North Direction Reef": "North Direction Island"}
FAST = ("Acropora", "Pocilloporidae")
LEVELS = {"min": 0.0, "q25": 0.25, "q50": 0.5, "q75": 0.75, "max": 1.0}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(run: Path) -> None:
    if run.exists() or OUT_DIR.exists():
        raise FileExistsError("Refusing to overwrite outputs")
    photo = pd.read_csv(DATA / "reef_photo_transect.csv", low_memory=False)
    manta = pd.read_csv(DATA / "reef_manta.csv")
    rows, levels = [], {}
    for reef, name in REEFS.items():
        p = photo[(photo.reef_name == name) & (photo.purpose == "COMPOSITION") &
                  (photo.variable == "HARD CORAL") & (photo.depth == 9)]
        share = p.groupby("report_year").apply(
            lambda d: d[d.reefpage_category.isin(FAST)]["median"].sum() / d["median"].sum(), include_groups=False)
        m = manta[(manta.reef_name == name) & (manta.data_type == "manta") & (manta.domain_category == "reef") &
                  (manta.variable == "HC") & (manta.purpose == "MANTA") & (manta.project_code == "LTMP") &
                  (manta.depth == 9)].groupby("report_year")["median"].mean()
        healthy = m[m >= m.median()].index
        for year, value in share.items():
            rows.append({"reef_name": reef, "report_year": int(year), "fast_share": float(value),
                         "manta_hc": float(m.get(year, float("nan"))), "healthy_year": bool(year in healthy)})
        kept = share[share.index.isin(healthy)]
        if len(kept) < 5:
            raise ValueError(f"Too few healthy transect years for {reef}")
        levels[reef] = {level: float(kept.quantile(q)) for level, q in LEVELS.items()} | {
            "n_years_kept": int(len(kept)), "n_years_all": int(len(share)), "manta_median": float(m.median())}
    levels["default"] = {level: sum(levels[r][level] for r in REEFS) / len(REEFS) for level in LEVELS}
    OUT_DIR.mkdir(parents=True)
    run.mkdir(parents=True)
    pd.DataFrame(rows).to_csv(OUT_DIR / "fast_share_by_year.csv", index=False)
    (OUT_DIR / "share_levels.json").write_text(json.dumps(levels, indent=2) + "\n")
    summary = {"status": "derived_input_not_promoted", "rule": __doc__.split("\n\n")[1],
               "levels": levels,
               "inputs": {p.name: sha256(p) for p in (DATA / "reef_photo_transect.csv", DATA / "reef_manta.csv", Path(__file__))},
               "outputs": {n: sha256(OUT_DIR / n) for n in ("fast_share_by_year.csv", "share_levels.json")}}
    (run / "initial_share_range.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(levels, indent=1))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    main(parser.parse_args().run.resolve())
