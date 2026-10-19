"""Food-scale evidence check (crash_mechanism_protocol.md, 2026-10-19, twenty-ninth round, test 2).

Compares reef-wide shallow preferred-prey cover from the CSIRO RHIS snorkel surveys with the LTMP 9 m
photo-transect fast cover (Acropora + Pocilloporidae, the model's COTS prey groups 1-3) in matched
reporting years.

RHIS proxy fast cover = live coral % x (branching + bushy + plate/table) / sum of all growth-form shares.
Growth form is not taxon (bushy/branching include Porites cylindrica, Seriatopora, Isopora), so the proxy is
an upper bound on Acropora + Pocilloporidae. RHIS live coral includes soft coral, which the growth-form share
removes from the numerator only. RHIS sites are chosen by the control programme (cull and surveillance
sites), not at random. Readout only: no model change follows (protocol, test 2, item 4).

Run with base miniconda Python (pandas, openpyxl):
    python analyse_rhis_fast_cover.py
Outputs (git-ignored; licensed source): sandbox/data/derived/rhis_fast_cover_lizard/
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "sandbox" / "data"
CSIRO = DATA / "260529-COTS-Manta-Cull-RHIS-Lawrence-CSIRO.xlsx"
PHOTO = DATA / "reef_photo_transect.csv"
OUT = DATA / "derived" / "rhis_fast_cover_lizard"
# model reef -> (RHIS ReefName pattern, LTMP reef_name)
REEFS = {
    "Lizard Island Reef": (r"^Lizard Island", "Lizard Isles"),
    "MacGillivray Reef": (r"^Macgillivray|^MacGillivray", "Macgillivray Reef"),
    "North Direction Reef": (r"^North Direction", "North Direction Island"),
    "Eyrie Reef": (r"^Eyrie", "Eyrie Reef"),
    "Martin Reef": (r"^Martin Reef", "Martin Reef"),
}
FORMS = ["Soft", "Branching", "Bushy", "PlateTable", "VaseFoliose", "Encrusting", "Mushroom", "Massive"]
FAST_FORMS = ["Branching", "Bushy", "PlateTable"]
FAST_TAXA = ("Acropora", "Pocilloporidae")
RATIO_SUPPORT = 1.5


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rhis_records() -> pd.DataFrame:
    d = pd.read_excel(CSIRO, sheet_name="RHIS")
    t = pd.to_datetime(d.SurveyTime)
    # LTMP reporting-year convention: surveys from September count to the next year.
    d["report_year"] = t.dt.year + (t.dt.month >= 9).astype(int)
    shares = d[[f"Coral{f}CoverProportion" for f in FORMS]].astype(float)
    total = shares.sum(axis=1)
    d["fast_form_share"] = shares[[f"Coral{f}CoverProportion" for f in FAST_FORMS]].sum(axis=1) / total
    d["hard_form_share"] = 1.0 - shares["CoralSoftCoverProportion"] / total
    d["live_coral"] = d.BenthosLiveCoralPercentage.astype(float) / 100.0
    d["fast_proxy"] = d.live_coral * d.fast_form_share
    d["hard_proxy"] = d.live_coral * d.hard_form_share
    rows = []
    for reef, (pattern, _) in REEFS.items():
        x = d[d.ReefName.str.contains(pattern, case=False, na=False, regex=True)].copy()
        x["reef_name"] = reef
        rows.append(x)
    out = pd.concat(rows, ignore_index=True)
    return out[np.isfinite(out.fast_proxy)]


def photo_fast() -> pd.DataFrame:
    p = pd.read_csv(PHOTO, low_memory=False)
    rows = []
    for reef, (_, name) in REEFS.items():
        x = p[(p.reef_name == name) & (p.purpose == "COMPOSITION") & (p.variable == "HARD CORAL") & (p.depth == 9)]
        for year, g in x.groupby("report_year"):
            rows.append({"reef_name": reef, "report_year": int(year),
                         "photo_fast_9m": float(g[g.reefpage_category.isin(FAST_TAXA)]["median"].sum()),
                         "photo_hard_9m": float(g["median"].sum())})
    return pd.DataFrame(rows)


def main() -> None:
    if OUT.exists():
        raise FileExistsError(f"Refusing to overwrite {OUT}")
    rhis = rhis_records()
    photo = photo_fast()
    by_habitat = (rhis.groupby(["reef_name", "report_year", "SurveyHabitat"])
                  .agg(n=("fast_proxy", "size"), rhis_fast=("fast_proxy", "mean"), rhis_hard=("hard_proxy", "mean"),
                       depth_mean=("Depth", "mean")).reset_index())
    by_year = (rhis.groupby(["reef_name", "report_year"])
               .agg(n=("fast_proxy", "size"), rhis_fast=("fast_proxy", "mean"),
                    rhis_fast_median=("fast_proxy", "median"), rhis_hard=("hard_proxy", "mean"),
                    rhis_fast_form_share=("fast_form_share", "mean")).reset_index())
    matched = by_year.merge(photo, on=["reef_name", "report_year"], how="inner")
    matched["ratio_fast"] = matched.rhis_fast / matched.photo_fast_9m
    matched["ratio_hard"] = matched.rhis_hard / matched.photo_hard_9m
    summary = {}
    for reef, g in matched.groupby("reef_name"):
        g = g[g.n >= 5]
        if g.empty:
            continue
        pooled = float(g.rhis_fast.mean() / g.photo_fast_9m.mean())
        summary[reef] = {"years": g.report_year.astype(int).tolist(), "records": int(g.n.sum()),
                         "rhis_fast_mean": float(g.rhis_fast.mean()), "photo_fast_9m_mean": float(g.photo_fast_9m.mean()),
                         "ratio_fast_pooled": pooled, "ratio_fast_by_year": g.ratio_fast.round(3).tolist(),
                         "ratio_hard_pooled": float(g.rhis_hard.mean() / g.photo_hard_9m.mean()),
                         "supports_broader_prey": pooled >= RATIO_SUPPORT}
    OUT.mkdir(parents=True)
    by_habitat.to_csv(OUT / "rhis_by_habitat.csv", index=False)
    by_year.to_csv(OUT / "rhis_by_year.csv", index=False)
    matched.to_csv(OUT / "rhis_vs_photo_matched.csv", index=False)
    report = {
        "status": "evidence_check_not_a_model_change",
        "protocol": "crash_mechanism_protocol.md, 2026-10-19, twenty-ninth round, test 2",
        "proxy": "RHIS live coral x (branching+bushy+plate/table)/sum(growth-form shares); upper proxy for Acropora+Pocilloporidae",
        "comparison": "LTMP 9 m photo-transect Acropora+Pocilloporidae (sum of COMPOSITION medians), same reporting year",
        "units": "fraction of benthos (RHIS visual estimate) vs fraction of transect points (photo)",
        "reporting_year": "RHIS survey month >= September counts to the next year (LTMP convention)",
        "matched_year_min_records": 5, "support_threshold_ratio": RATIO_SUPPORT,
        "summary": summary,
        "licence": "CSIRO workbook internal/licensed: outputs git-ignored, summary statistics only in the protocol",
        "inputs": {"csiro_workbook": sha256(CSIRO), "photo_transect": sha256(PHOTO)},
        "script_sha256": sha256(Path(__file__)), "created": datetime.now().isoformat(timespec="seconds"),
        "outputs": {n: sha256(OUT / n) for n in ("rhis_by_habitat.csv", "rhis_by_year.csv", "rhis_vs_photo_matched.csv")},
    }
    (OUT / "rhis_fast_cover.json").write_text(json.dumps(report, indent=2) + "\n")
    pd.set_option("display.width", 220)
    print(matched.round(3).to_string(index=False))
    print(by_habitat[by_habitat.reef_name == "Lizard Island Reef"].round(3).to_string(index=False))
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
