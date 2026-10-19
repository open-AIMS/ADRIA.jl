"""Juvenile-to-adult pipeline diagnostic, 2010-26 (lead's decision 3, 2026-10-19, twenty-first round).

Where does the model lose the cohort that became the observed 2020s adults?

Model: W-next arms (wc_b100_q25_..._aprmak_cv2_pnd, corr_w3, run 20261016T_03, five seeds; seed means)
at Lizard, Eyrie, MacGillivray and North Direction, from mechanisms.csv / flows.csv:
settled recruits (age 0), subadults (age 1, N2), adults (age 2+), maturation, juvenile predator-pit
removals, food survival, fast coral. "Maturation yield" = adults matured in year t / subadults at the
end of year t-1 (reef aggregates; it combines background juvenile survival, the pit and food survival).

Observed (summary statistics only; the source workbooks are licensed and stay unpublished):
- culls (cots_cull_removals_lizard_v0_3_v2/cull_reef_year.csv): culled per cohort (Cohort1 <15 cm,
  Cohort2 15-<25, Cohort3 25-40, Cohort4 >40 cm) and catch per hour of bottom time; Lizard = 14-116b/c/d;
- SALAD diver tracks ('COTS recorded'): n, median and 10th percentile diameter, share < 150 / < 250 mm.
Size-at-age is approximate: COTS are herbivorous at < 1-2 cm for at least 6 months, then corallivorous,
reaching about 15-25 cm at age 2 and > 25 cm at age 3+ (e.g. Pratchett et al. 2014 review; Deaker et al.
2020 show juveniles can stay at 16-18 mm for years without coral). So adults > 25 cm in year y
settled in about y-3 or earlier.

    python analyse_juvenile_pipeline.py [--run <run_dir>]
Outputs: sandbox/data/derived/juvenile_pipeline/ (git-ignored; never overwritten).
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "sandbox" / "data"
CULL = DATA / "derived" / "cots_cull_removals_lizard_v0_3_v2" / "cull_reef_year.csv"
SALAD = DATA / "COTS SALAD MASTER 2021-2026.xlsx"
OUT = DATA / "derived" / "juvenile_pipeline"
DEFAULT_RUN = ROOT / "sandbox/calibration/runs/20261016T_03_pit_q_w3_1991_v1"
ARM = "wc_b100_q25_eb_h10_obs015_s{k}_t60_e50_aprmak_cv2_pnd"
REEFS = {  # model reef: (cull reef_id prefix, SALAD reef name)
    "Lizard Island Reef": ("14-116", "Lizard"),
    "Eyrie Reef": ("14-118", None),
    "MacGillivray Reef": ("14-114", "MacGillivray"),
    "North Direction Reef": ("14-143", "North Direction"),
}
COHORTS = ("Cohort1", "Cohort2", "Cohort3", "Cohort4")
YEARS = range(2010, 2027)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def model_table(run: Path) -> pd.DataFrame:
    arms = [ARM.format(k=k) for k in range(1, 6)]
    mech = pd.read_csv(run / "mechanisms.csv")
    flows = pd.read_csv(run / "flows.csv")
    mech = mech[mech.treatment.isin(arms) & mech.reef_name.isin(REEFS)]
    flows = flows[flows.treatment.isin(arms) & flows.reef_name.isin(REEFS)]
    if mech.treatment.nunique() != 5:
        raise ValueError("Expected five W-next seed arms")
    m = mech.merge(flows[["treatment", "reef_name", "year", "settled_recruits", "internal_immigration", "external_immigration"]],
                   on=["treatment", "reef_name", "year"])
    m = m.sort_values(["treatment", "reef_name", "year"])
    m["subadults_prev"] = m.groupby(["treatment", "reef_name"]).subadults_ha.shift(1)
    m["maturation_yield"] = m.maturation_ha / m.subadults_prev
    cols = ["settled_recruits", "internal_immigration", "external_immigration", "recruits_ha", "subadults_ha", "adults_ha",
            "maturation_ha", "juvenile_predation_ha", "maturation_yield", "food_survival", "fast_coral", "total_coral"]
    return m[m.year.isin(YEARS)].groupby(["reef_name", "year"])[cols].mean().reset_index()


def observed_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    cull = pd.read_csv(CULL)
    rows = []
    for reef, (prefix, _) in REEFS.items():
        c = cull[cull.reef_id.str.startswith(prefix)].groupby("calendar_year")[list(COHORTS) + ["culled", "bottomtime_min"]].sum()
        for y, r in c.iterrows():
            n = max(r.culled, 1)
            rows.append({"reef_name": reef, "year": int(y), "culled": int(r.culled), "bottom_hours": r.bottomtime_min / 60,
                         "catch_per_hour": r.culled / (r.bottomtime_min / 60) if r.bottomtime_min > 0 else np.nan,
                         **{f"share_{c.lower()}": r[c] / n for c in COHORTS},
                         **{f"catch_per_hour_{c.lower()}": r[c] / (r.bottomtime_min / 60) if r.bottomtime_min > 0 else np.nan
                            for c in COHORTS}})
    cull_t = pd.DataFrame(rows)
    d = pd.read_excel(SALAD, sheet_name="COTS recorded")
    d = d[d["CoTS/ Scar"].astype(str).str.strip().str.lower() == "cots"]
    d["size_mm"] = pd.to_numeric(d["Size(mm)"], errors="coerce")
    rows = []
    for reef, (_, salad) in REEFS.items():
        if salad is None:
            continue
        for y, g in d[d.Reef == salad].groupby("Year"):
            s = g.size_mm.dropna()
            if len(s) == 0:
                continue
            rows.append({"reef_name": reef, "year": int(y), "n_sized": len(s), "median_mm": float(s.median()),
                         "p10_mm": float(s.quantile(0.1)), "share_lt150": float((s < 150).mean()),
                         "share_lt250": float((s < 250).mean())})
    return cull_t, pd.DataFrame(rows)


def main(run: Path) -> None:
    if OUT.exists():
        raise FileExistsError(f"Refusing to overwrite {OUT}")
    model = model_table(run)
    cull, salad = observed_tables()
    OUT.mkdir(parents=True)
    paths = {n: OUT / n for n in ("pipeline_model.csv", "pipeline_cull.csv", "pipeline_salad_sizes.csv", "pipeline.png")}
    model.to_csv(paths["pipeline_model.csv"], index=False)
    cull.to_csv(paths["pipeline_cull.csv"], index=False)
    salad.to_csv(paths["pipeline_salad_sizes.csv"], index=False)

    fig, axes = plt.subplots(2, len(REEFS), figsize=(16, 7.5), sharex=True)
    for col, reef in enumerate(REEFS):
        g = model[model.reef_name == reef].set_index("year")
        ax = axes[0, col]
        for name, color in (("settled_recruits", "#86b6ef"), ("subadults_ha", "#3987e5"), ("adults_ha", "#0d366b")):
            ax.plot(g.index, g[name], color=color, linewidth=1.9, label={"settled_recruits": "settled recruits (age 0)",
                    "subadults_ha": "subadults (age 1)", "adults_ha": "adults (age 2+)"}[name])
        ax.set_yscale("log")
        ax.set_ylim(0.1, 400)
        ax.set_title(reef, loc="left", fontsize=10)
        ax2 = ax.twinx()
        ax2.plot(g.index, g.food_survival, color="#eb6834", linewidth=1.2, linestyle="--", label="food survival")
        ax2.set_ylim(0, 1.05)
        ax2.spines[["top"]].set_visible(False)
        if col == len(REEFS) - 1:
            ax2.set_ylabel("Food survival (dashed)")
        ax = axes[1, col]
        c = cull[cull.reef_name == reef].set_index("year")
        bottom = np.zeros(len(c))
        for cohort, color in zip(COHORTS, ("#cfe3fa", "#86b6ef", "#3987e5", "#0d366b")):
            v = c[f"catch_per_hour_{cohort.lower()}"].to_numpy()
            ax.bar(c.index, v, bottom=bottom, color=color, width=0.7, edgecolor="#fcfcfb", linewidth=1,
                   label={"Cohort1": "< 15 cm", "Cohort2": "15-25 cm", "Cohort3": "25-40 cm", "Cohort4": "> 40 cm"}[cohort])
            bottom += np.nan_to_num(v)
        s = salad[salad.reef_name == reef].set_index("year")
        if len(s):
            ax3 = ax.twinx()
            ax3.plot(s.index, s.median_mm, "o-", color="#eb6834", markersize=5, linewidth=1.2, label="SALAD median diameter")
            ax3.set_ylim(0, 550)
            if col == len(REEFS) - 1:
                ax3.set_ylabel("SALAD median diameter (mm)")
        for a in axes[:, col]:
            a.grid(color="#e4e3df", linewidth=.8)
            a.spines[["top", "right"]].set_visible(False)
            a.axvspan(2016, 2017.5, color="#f3d9cc", alpha=0.5, linewidth=0)
    axes[0, 0].set_ylabel("Model COTS ha$^{-1}$ (log)")
    axes[1, 0].set_ylabel("Culled per hour of bottom time")
    axes[0, 0].legend(frameon=False, fontsize=8, loc="lower left")
    axes[1, 0].legend(frameon=False, fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(paths["pipeline.png"], dpi=150, facecolor="#fcfcfb")
    plt.close(fig)

    summary = {}
    for reef in REEFS:
        g = model[model.reef_name == reef].set_index("year")
        summary[reef] = {
            "model_maturation_yield_2012_14": float(g.loc[2012:2014].maturation_yield.mean()),
            "model_maturation_yield_2015_18": float(g.loc[2015:2018].maturation_yield.mean()),
            "model_maturation_yield_2020_24": float(g.loc[2020:2024].maturation_yield.mean()),
            "model_food_survival_min_2015_19": float(g.loc[2015:2019].food_survival.min()),
            "model_subadults_2015_17": float(g.loc[2015:2017].subadults_ha.mean()),
            "model_matured_2015_18_total": float(g.loc[2015:2018].maturation_ha.sum()),
            "model_pit_removed_2015_18_total": float(g.loc[2015:2018].juvenile_predation_ha.sum()),
            "model_first_year_adults_ge_5": next((int(y) for y, v in g.loc[2019:].adults_ha.items() if v >= 5), None)}
    meta = {"status": "diagnostic", "summary": summary,
            "inputs": {"mechanisms.csv": sha256(run / "mechanisms.csv"), "flows.csv": sha256(run / "flows.csv"),
                       "cull_reef_year.csv": sha256(CULL), "salad_workbook": sha256(SALAD)},
            "run": str(run.relative_to(ROOT)), "arms": ARM, "script_sha256": sha256(Path(__file__)),
            "outputs": {n: sha256(p) for n, p in paths.items()}}
    (OUT / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
    pd.set_option("display.width", 250)
    print(json.dumps(summary, indent=1))
    print(cull[cull.year >= 2018].round(2).to_string(index=False))
    print(salad.round(2).to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, default=DEFAULT_RUN)
    main(parser.parse_args().run.resolve())
