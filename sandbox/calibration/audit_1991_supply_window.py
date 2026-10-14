"""Audit the 2017-2024 external-supply counterfactual (crash_mechanism_protocol.md,
"Preregistered diagnostic: 2017-2024 external-supply counterfactual", 2026-10-14) against the
diver-anchored 2020s density targets (build_2020s_density_targets.py; git-ignored).

    python audit_1991_supply_window.py <run> <rescore_q015> <gate_q015>
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from audit_1991_initiation import exact_replay  # noqa: E402
from audit_1991_switching_repair import REEFS, sha256  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
PRIOR = ROOT / "sandbox/calibration/runs/20261014T_37_candidates_1991_v1"
TARGETS = ROOT / "sandbox/data/derived/density_targets_2020s/targets_reef_year.csv"
BASE = "wc_b150_q20_eb_lin_obs015_s1"
REFERENCE, PRIOR_REFERENCE = f"{BASE}_sw2017_2024_x001_cv2_pnd", f"{BASE}_cv2_pnd"
WINDOW = (2020, 2024)
OFFSET = 0.5  # adults/ha added before the log ratio (keeps near-zero years finite)


def main(run: Path, rescore: Path, gate: Path) -> None:
    paths = {n: run / n for n in ("supply_window_audit.csv", "supply_window_audit.json", "supply_window.png")}
    if any(p.exists() for p in paths.values()):
        raise FileExistsError("Refusing to overwrite supply-window audit")
    tables = {n: pd.read_csv(run / f"{n}.csv") for n in ("trajectories", "flows", "mechanisms", "failures")}
    if len(tables["failures"]):
        raise ValueError("Failed evaluations in run")
    trajectories, flows, mechanisms = (tables[k] for k in ("trajectories", "flows", "mechanisms"))
    replay = {}
    for name in ("trajectories", "flows", "mechanisms"):
        prior = pd.read_csv(PRIOR / f"{name}.csv")
        columns = [c for c in prior if c not in ("treatment", "reef_name", "year")]
        replay[name] = exact_replay(tables[name][tables[name].treatment == REFERENCE],
                                    prior[prior.treatment == PRIOR_REFERENCE], columns)
    targets = pd.read_csv(TARGETS)
    targets = targets[targets.report_year.between(*WINDOW) & targets.target_reef.isin(REEFS)]
    scores = pd.read_csv(rescore / "detectability_scores.csv")
    g2 = pd.read_csv(gate / "gate_g2prime.csv").set_index(["treatment", "a_half"])
    rows = []
    for arm in sorted(set(trajectories.treatment)):
        s = scores[(scores.run == run.name) & (scores.treatment == arm) & (scores.a_half == 0.0)]
        primary = s[s.score_v2_primary.astype(str).str.lower() == "true"]
        traj = trajectories[trajectories.treatment == arm]
        joined = targets.merge(traj, left_on=["target_reef", "report_year"], right_on=["reef_name", "year"])
        log_err = np.log((joined.adults_ha + OFFSET) / (joined.target_density_ha + OFFSET))
        row = {"treatment": arm, "factor": int(arm.split("_x")[1][:3]), "window": arm.split("_sw")[1][:9],
               "summed_loss_v2": float(primary.loss_v2.sum()),
               "g2_pass_ahalf30": bool(g2.loc[(arm, 30.0)].passes_g2prime),
               "g2_pass_ahalf10": bool(g2.loc[(arm, 10.0)].passes_g2prime),
               "target_reef_years": int(len(joined)),
               "mean_log_ratio_2020s": float(log_err.mean()),
               "mean_abs_log_ratio_2020s": float(log_err.abs().mean())}
        for reef in REEFS:
            short = {"Lizard Island Reef": "liz", "MacGillivray Reef": "macg", "North Direction Reef": "nd", "Eyrie Reef": "eyr"}[reef]
            t = traj[traj.reef_name == reef].set_index("year")
            f = flows[(flows.treatment == arm) & (flows.reef_name == reef)].set_index("year")
            m = mechanisms[(mechanisms.treatment == arm) & (mechanisms.reef_name == reef)].set_index("year")
            above = t.loc[2017:2024][t.loc[2017:2024].adults_ha >= 3.0]
            row |= {f"{short}_adults_2020_2024": ";".join(f"{v:.1f}" for v in t.loc[2020:2024].adults_ha),
                    f"{short}_first_year_ge3": int(above.index.min()) if len(above) else np.nan,
                    f"{short}_settled_2017_2024": float(f.loc[2017:2024].settled_recruits.sum()),
                    f"{short}_juv_predation_2017_2024": float(m.loc[2017:2024].juvenile_predation_ha.sum()),
                    f"{short}_maturation_2017_2024": float(m.loc[2017:2024].maturation_ha.sum()),
                    f"{short}_food_survival_mean_2017_2024": float(m.loc[2017:2024].food_survival.mean()),
                    f"{short}_total_coral_2024": float(m.loc[2024].total_coral)}
        rows.append(row)
    result = pd.DataFrame(rows).sort_values(["window", "factor"]).reset_index(drop=True)
    result.to_csv(paths["supply_window_audit.csv"], index=False)

    fig, axes = plt.subplots(4, 1, figsize=(10, 12), sharex=True)
    colors = {1: "#52514e", 3: "#86b6ef", 10: "#3987e5", 30: "#0d366b"}
    all_targets = pd.read_csv(TARGETS)
    for i, reef in enumerate(REEFS):
        axis = axes[i]
        for arm in result.treatment:
            factor = int(arm.split("_x")[1][:3])
            late = "2019_2024" in arm
            x = trajectories[(trajectories.treatment == arm) & (trajectories.reef_name == reef) & (trajectories.year >= 2008)]
            axis.plot(x.year, x.adults_ha, color=colors[factor], linestyle=":" if late else "-", linewidth=1.6,
                      label=(f"x{factor}" + (" (2019-24)" if late else " (2017-24)")) if i == 0 else None)
        tg = all_targets[all_targets.target_reef == reef]
        for source, marker in (("SALAD", "o"), ("cull_converted", "^")):
            z = tg[tg.target_source == source]
            axis.scatter(z.report_year, z.target_density_ha, marker=marker, s=28, color="#eb6834", zorder=5,
                         label=("SALAD" if source == "SALAD" else "cull, converted") if i == 0 else None)
        axis.set_ylabel(f"{reef}\nadults / ha")
        axis.set_yscale("symlog", linthresh=1.0)
        axis.grid(color="#e4e3df", linewidth=.8)
    axes[0].set_title("External-supply counterfactual vs diver-anchored 2020s densities (labelled diagnostic)",
                      loc="left", fontsize=10)
    axes[0].legend(fontsize=8, frameon=False, ncol=4)
    axes[-1].set_xlim(2008, 2026.5)
    fig.tight_layout()
    fig.savefig(paths["supply_window.png"], dpi=160, facecolor="#fcfcfb")
    plt.close(fig)

    report = {"status": "labelled_counterfactual_diagnostic", "reference_replay_max_rel_delta": replay,
              "score_2020s": f"mean log((model + {OFFSET}) / (target + {OFFSET})) over reef-years {WINDOW}",
              "summary": json.loads(result.to_json(orient="records")),
              "inputs": {n: sha256(run / n) for n in ("metadata.toml", "trajectories.csv", "flows.csv", "mechanisms.csv")}
                        | {"targets_reef_year.csv": sha256(TARGETS),
                           "detectability_scores.csv": sha256(rescore / "detectability_scores.csv"),
                           "gate_g2prime.csv": sha256(gate / "gate_g2prime.csv")},
              "script_sha256": sha256(Path(__file__)),
              "outputs": {n: sha256(p) for n, p in paths.items() if n != "supply_window_audit.json"}}
    paths["supply_window_audit.json"].write_text(json.dumps(report, indent=2, sort_keys=True, default=str) + "\n")
    pd.set_option("display.width", 320)
    pd.set_option("display.max_columns", 60)
    print(result.T.to_string())
    print("Replay:", replay)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("rescore", type=Path)
    parser.add_argument("gate", type=Path)
    a = parser.parse_args()
    main(a.run.resolve(), a.rescore.resolve(), a.gate.resolve())
