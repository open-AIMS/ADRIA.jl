"""Audit the predator-pit capacity x correction weight test (crash_mechanism_protocol.md,
"Preregistered predator-pit capacity x correction weight", 2026-10-16, eighteenth round):
W + h10 x Q {10, 15, 20, 25} x five seed pairs, on the face-value (w 1), corr_w2 and corr_w3 layers.

    python audit_1991_pit_q.py <run_w1> <rescore_w1> <gate_w1> <run_w2> <rescore_w2> <gate_w2> \\
        <run_w3> <rescore_w3> <gate_w3>
Outputs go to <run_w1>.
"""

from __future__ import annotations

import argparse
import json
import re
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
RUNS = ROOT / "sandbox/calibration/runs"
TARGETS = ROOT / "sandbox/data/derived/density_targets_2020s/targets_reef_year.csv"
PRIOR = "20261015T_42_facevalue_h10_1991_v1"
ARM = re.compile(r"^wc_b100_q(\d+)_eb_h10_obs015_s(\d)_t60_e50_aprmak_cv2_pnd$")
SHORT = {"Lizard Island Reef": "liz", "MacGillivray Reef": "macg", "North Direction Reef": "nd", "Eyrie Reef": "eyr"}
BASE = (1, 20)  # current configuration (w, Q)
YEARS = (2022, 2023, 2024, 2025, 2026)


def main(sets: dict[int, tuple[Path, Path, Path]]) -> None:
    out = sets[1][0]
    paths = {n: out / n for n in ("pit_q_audit.csv", "pit_q_summary.csv", "pit_q_audit.json", "pit_q.png")}
    if any(p.exists() for p in paths.values()):
        raise FileExistsError("Refusing to overwrite pit_q audit")
    targets = pd.read_csv(TARGETS)
    rows, series, replay = [], {}, {}
    for w, (run, rescore, gate) in sets.items():
        tables = {n: pd.read_csv(run / f"{n}.csv") for n in ("trajectories", "flows", "mechanisms", "failures")}
        if len(tables["failures"]):
            raise ValueError(f"Failed evaluations in {run.name}")
        arms = sorted(set(tables["trajectories"].treatment))
        if len(arms) != 20 or not all(ARM.match(a) for a in arms):
            raise ValueError(f"Unexpected arm set in {run.name}")
        if w == 1:
            for name in ("trajectories", "flows", "mechanisms"):
                prior = pd.read_csv(RUNS / PRIOR / f"{name}.csv")
                columns = [c for c in prior if c not in ("treatment", "reef_name", "year")]
                for arm in (a for a in arms if ARM.match(a).group(1) == "20"):
                    replay[f"{arm}/{name}"] = exact_replay(tables[name][tables[name].treatment == arm],
                                                           prior[prior.treatment == arm], columns)
        scores = pd.read_csv(rescore / "detectability_scores.csv")
        g2 = pd.read_csv(gate / "gate_g2prime.csv").set_index(["treatment", "a_half"])
        traj, mech = tables["trajectories"], tables["mechanisms"]
        for arm in arms:
            q, seed = int(ARM.match(arm).group(1)), int(ARM.match(arm).group(2))
            s = scores[(scores.treatment == arm) & (scores.a_half == 0.0)]
            primary = s[s.score_v2_primary.astype(str).str.lower() == "true"].set_index("reef_name")
            gr = g2.loc[(arm, 30.0)]
            row = {"w": w, "Q": q, "seed": seed, "treatment": arm, "loss_v3": float(primary.loss_v3.sum()),
                   "loss_v2_manta": float(primary.loss_v2_manta.sum()), "diver_error": float(primary.diver_error.sum()),
                   "g2_pass_ahalf30": bool(gr.passes_g2prime),
                   "liz_trough_cpue": float(gr.liz_trough_model_mean_cpue), "liz_trough_bound": float(gr.liz_trough_observed_upper_95),
                   "macg_trough_cpue": float(gr.macg_trough_model_mean_cpue), "macg_trough_bound": float(gr.macg_trough_observed_upper_95)}
            for reef, short in SHORT.items():
                a = traj[(traj.treatment == arm) & (traj.reef_name == reef)].set_index("year").adults_ha
                series[(w, q, seed, reef)] = a
                row[f"{short}_diver_error"] = float(primary.loc[reef, "diver_error"])
                tg = targets[targets.target_reef == reef].set_index("report_year").target_density_ha
                for y in YEARS:
                    row[f"{short}_{y}_model"] = float(a.get(y, np.nan))
                    row[f"{short}_{y}_target"] = float(tg.get(y, np.nan))
            m = mech[(mech.treatment == arm) & (mech.reef_name == "Lizard Island Reef")].set_index("year")
            sub_prev = m.subadults_ha.shift(1)
            row["liz_pit_share_2020_23"] = float((m.juvenile_predation_ha / sub_prev).loc[2020:2023].mean())
            rows.append(row)
    result = pd.DataFrame(rows).sort_values(["w", "Q", "seed"]).reset_index(drop=True)
    result.to_csv(paths["pit_q_audit.csv"], index=False)
    g = result.groupby(["w", "Q"])
    summary = g.agg(v3_mean=("loss_v3", "mean"), v3_min=("loss_v3", "min"), v3_max=("loss_v3", "max"),
                    manta_mean=("loss_v2_manta", "mean"), diver_mean=("diver_error", "mean"),
                    g2_passes_ahalf30=("g2_pass_ahalf30", "sum"), liz_trough_cpue=("liz_trough_cpue", "mean"),
                    macg_trough_cpue=("macg_trough_cpue", "mean"), liz_pit_share_2020_23=("liz_pit_share_2020_23", "mean"))
    for short in SHORT.values():
        summary[f"{short}_diver_error_mean"] = g[f"{short}_diver_error"].mean()
        for y in YEARS:
            summary[f"{short}_{y}_model_mean"] = g[f"{short}_{y}_model"].mean()
            ratio = np.log((result[f"{short}_{y}_model"] + 0.5) / (result[f"{short}_{y}_target"] + 0.5))
            summary[f"{short}_{y}_model_over_target_gmean"] = np.exp(ratio.groupby([result.w, result.Q]).mean())
    summary = summary.reset_index()
    base = summary[(summary.w == BASE[0]) & (summary.Q == BASE[1])].iloc[0]
    summary["qualifies"] = (summary.g2_passes_ahalf30 >= 4) & (summary.v3_mean < base.v3_mean) & (summary.diver_mean < base.diver_mean)
    summary.to_csv(paths["pit_q_summary.csv"], index=False)

    colors = {10: "#0d366b", 15: "#3987e5", 20: "#86b6ef", 25: "#b6b5b0"}
    fig, axes = plt.subplots(4, 3, figsize=(16, 13), sharex=True)
    for col, w in enumerate(sorted(sets)):
        for r, reef in enumerate(REEFS):
            ax = axes[r, col]
            for q in sorted(colors):
                m = pd.concat([series[(w, q, s, reef)] for s in range(1, 6)], axis=1).loc[2012:]
                ax.plot(m.index, m.mean(axis=1), color=colors[q], linewidth=2.2 if q == 20 else 1.6,
                        label=f"pit capacity {q}" if r == 0 else None)
            tg = targets[(targets.target_reef == reef) & (targets.report_year >= 2012)]
            ax.scatter(tg.report_year, tg.target_density_ha, color="#eb6834", s=22, zorder=5,
                       label="diver/cull target" if r == 0 else None)
            ax.set_yscale("symlog", linthresh=1.0)
            ax.set_ylim(bottom=0)
            ax.grid(color="#e4e3df", linewidth=.8)
            ax.spines[["top", "right"]].set_visible(False)
            if col == 0:
                ax.set_ylabel(f"{reef}\nadults / ha")
        axes[0, col].set_title(f"correction weight w = {w} (seed means)", loc="left", fontsize=10)
        axes[0, col].legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(paths["pit_q.png"], dpi=150, facecolor="#fcfcfb")
    plt.close(fig)

    report = {"status": "diagnostic_not_adopted", "reference_replay_max_rel_delta": replay,
              "selection_rule": "G2' A_half 30 >= 4/5 seeds; seed-mean v3 and diver error both below (w 1, Q 20)",
              "summary": json.loads(summary.to_json(orient="records")),
              "inputs": {f"w{w}/{p.name}/{n}": sha256(p / n) for w, (run, rescore, gate) in sets.items()
                         for p, n in ((run, "metadata.toml"), (run, "trajectories.csv"), (run, "mechanisms.csv"),
                                      (rescore, "detectability_scores.csv"), (gate, "gate_g2prime.csv"))}
                        | {"targets_reef_year.csv": sha256(TARGETS)},
              "script_sha256": sha256(Path(__file__)),
              "outputs": {n: sha256(p) for n, p in paths.items() if n != "pit_q_audit.json"}}
    paths["pit_q_audit.json"].write_text(json.dumps(report, indent=2, sort_keys=True, default=str) + "\n")
    pd.set_option("display.width", 250)
    print(summary.set_index(["w", "Q"]).T.round(3).to_string())
    print(result.pivot_table(index="seed", columns=["w", "Q"], values="loss_v3").round(2).to_string())
    print("Replay max:", max(replay.values()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for w in (1, 2, 3):
        for k in ("run", "rescore", "gate"):
            parser.add_argument(f"{k}_w{w}", type=Path)
    a = vars(parser.parse_args())
    main({w: tuple(a[f"{k}_w{w}"].resolve() for k in ("run", "rescore", "gate")) for w in (1, 2, 3)})
