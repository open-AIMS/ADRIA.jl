"""Audit the tolerance-spread x low-DHW threshold screen (crash_mechanism_protocol.md,
"Preregistered tolerance-spread x low-DHW threshold screen", 2026-10-15).

Bleaching-attributable relative loss in year y (Lizard, MacGillivray, North Direction mean):
    L_y = 1 - cover_y / cover_y(bleaching off, same mu)

    python audit_1991_coral_spread.py <run>
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from audit_1991_coral_calibration import observed_photo  # noqa: E402
from audit_1991_growth_starvation import PHOTO  # noqa: E402
from audit_1991_initiation import exact_replay  # noqa: E402
from audit_1991_joint_coral import objective, reef_errors  # noqa: E402
from audit_1991_switching_repair import REEFS, sha256  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / "sandbox/calibration/runs"
REPLAYS = {"spread_s100_t40_m120": ("20261015T_01_coral_thermal_1991_v1", "therm_b100_d000_h030_m120_p100"),
           "spread_off_t40_m120": ("20261015T_02_coral_thermal_nobleach_1991_v1", "therm_b000_d000_h030_m120_p100")}
RECOVERY_REEFS = ("Lizard Island Reef", "MacGillivray Reef", "North Direction Reef")
EVENT_YEARS = (2016, 2017, 2020, 2022, 2024)
REFERENCE = "spread_s100_t40_m120"


def main(run: Path) -> None:
    paths = {n: run / n for n in ("coral_spread_audit.csv", "coral_spread_audit.json")}
    if any(p.exists() for p in paths.values()):
        raise FileExistsError("Refusing to overwrite coral-spread audit")
    tables = {n: pd.read_csv(run / f"{n}.csv") for n in ("trajectories", "mechanisms", "failures")}
    if len(tables["failures"]):
        raise ValueError("Failed evaluations in run")
    mech = tables["mechanisms"]
    replay = {}
    for current, (prior_run, prior_arm) in REPLAYS.items():
        for name in ("trajectories", "mechanisms"):
            prior = pd.read_csv(RUNS / prior_run / f"{name}.csv")
            columns = [c for c in prior if c not in ("treatment", "reef_name", "year")]
            replay[f"{current}/{name}"] = exact_replay(tables[name][tables[name].treatment == current],
                                                       prior[prior.treatment == prior_arm], columns)
    photo = observed_photo()
    manta_raw = pd.read_csv(ROOT / "sandbox/data/reef_manta.csv")
    manta_raw = manta_raw[(manta_raw.data_type == "manta") & (manta_raw.domain_category == "reef") &
                          (manta_raw.variable == "HC") & (manta_raw.purpose == "MANTA") &
                          (manta_raw.project_code == "LTMP") & (manta_raw.depth == 9)]
    manta = {r: manta_raw[(manta_raw.reef_name == REEFS[r]) & manta_raw.report_year.between(1992, 2024)]
             .groupby("report_year")["median"].mean() for r in REEFS}
    errs = {"all": [], "pre": [], "post": []}
    windows = {"all": (1992, 2024), "pre": (1992, 2016), "post": (2018, 2024)}
    for (name, reef), m in mech.groupby(["treatment", "reef_name"]):
        m = m.set_index("year").assign(slow_coral=lambda x: x.total_coral - x.fast_coral)
        for k, w in windows.items():
            errs[k].append({"treatment": name, "reef_name": reef} | reef_errors(m, reef, photo, manta[reef], w))
    j = {k: objective(pd.DataFrame(v)) for k, v in errs.items()}
    post = pd.DataFrame(errs["post"])
    post = post[post.reef_name.isin(RECOVERY_REEFS)].groupby("treatment")[["rmse_total_vs_manta", "rmse_fast_vs_photo"]].mean().mean(axis=1)

    sub = mech[mech.reef_name.isin(RECOVERY_REEFS)].groupby(["treatment", "year"])[["fast_coral", "total_coral"]].mean()
    rows = []
    for arm in sorted(set(mech.treatment)):
        p = arm.split("_")
        off = arm.startswith("spread_off")
        mu = int(p[3][1:]) / 1000
        cf = f"spread_off_t40_m{p[3][1:]}"
        row = {"treatment": arm, "s": np.nan if off else int(p[1][1:]) / 100, "threshold": int(p[2][1:]) / 10,
               "mu": mu, "bleaching_off": off, "J_all": float(j["all"][arm]), "J_pre2017": float(j["pre"][arm]),
               "post_recovery_error": float(post[arm])}
        for y in EVENT_YEARS:
            for c in ("fast_coral", "total_coral"):
                row[f"loss_{c.split('_')[0]}_{y}"] = float(1 - sub.loc[(arm, y), c] / sub.loc[(cf, y), c])
        row["loss_fast_2016_17"] = max(row["loss_fast_2016"], row["loss_fast_2017"])
        row["loss_fast_2020_24"] = float(np.mean([row["loss_fast_2020"], row["loss_fast_2024"]]))
        rows.append(row)
    result = pd.DataFrame(rows)
    ref = result.set_index("treatment").loc[REFERENCE]
    cand = result[~result.bleaching_off & (result.loss_fast_2016_17 >= 0.5 * ref.loss_fast_2016_17)]
    best = cand.sort_values(["loss_fast_2020_24", "J_all"]).iloc[0].treatment if len(cand) else None
    result["keeps_2016_17"] = result.treatment.isin(cand.treatment)
    result.to_csv(paths["coral_spread_audit.csv"], index=False)
    report = {"status": "bounded_diagnostic_not_promoted", "reference_replay_max_rel_delta": replay,
              "selection": {"rule": "keep >= half the reference 2016-17 attributable fast loss; min mean 2020/2024 attributable fast loss; then J",
                            "eligible": cand.treatment.tolist(), "selected": best},
              "summary": json.loads(result.to_json(orient="records")),
              "inputs": {n: sha256(run / n) for n in ("metadata.toml", "trajectories.csv", "mechanisms.csv")}
                        | {"reef_photo_transect.csv": sha256(PHOTO)},
              "script_sha256": sha256(Path(__file__)),
              "outputs": {"coral_spread_audit.csv": sha256(paths["coral_spread_audit.csv"])}}
    paths["coral_spread_audit.json"].write_text(json.dumps(report, indent=2, sort_keys=True, default=str) + "\n")
    pd.set_option("display.width", 300)
    pd.set_option("display.max_columns", 40)
    print(result[["treatment", "mu", "J_all", "J_pre2017", "post_recovery_error", "loss_fast_2016", "loss_fast_2017",
                  "loss_fast_2020", "loss_fast_2022", "loss_fast_2024", "loss_total_2024", "keeps_2016_17"]]
          .sort_values(["mu", "treatment"]).round(3).to_string(index=False))
    print("Replay:", replay, "| selected:", best)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    main(parser.parse_args().run.resolve())
