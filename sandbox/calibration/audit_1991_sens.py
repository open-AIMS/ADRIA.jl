"""Audit promotion gate 5, the +/-25% sensitivity (crash_mechanism_protocol.md, "Preregistered
promotion checks", 2026-10-19): W + h10, Q 20, Ext x5, corr_w2, one parameter at a time at
x0.75 / x1.25, five seed pairs. Q 15 / 25 come from the pit_q run (20261016T_02, re-score 06,
gate 07).
W-next (2026-10-19, twenty-first round): --base-q 25 audits the corr_w3 run with Q 25 x f in the same
run (_sxq<f>_); its base arms replay the q25 arms of 20261016T_03 (pit_q, corr_w3).

Pass (preregistered): each perturbation passes G2' (A_half 30, trough on Lizard and MacGillivray)
in >= 4/5 seeds and its seed-mean v3 is no more than 1.0 above the base.

    python audit_1991_sens.py <run> <rescore> <gate> [--base-q 25]
Outputs go to <run>.
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
from audit_1991_switching_repair import sha256  # noqa: E402

RUNS = Path(__file__).resolve().parent / "runs"
PITQ = ("20261016T_02_pit_q_w2_1991_v1", "20261016T_06_rescore_pit_q_w2", "20261016T_07_gate_pit_q_w2_lizmacg")
PRIOR = {20: PITQ[0], 25: "20261016T_03_pit_q_w3_1991_v1"}
TARGETS = Path(__file__).resolve().parents[2] / "sandbox/data/derived/density_targets_2020s/targets_reef_year.csv"
SHORT = {"Lizard Island Reef": "liz", "MacGillivray Reef": "macg", "North Direction Reef": "nd", "Eyrie Reef": "eyr"}
V3_TOLERANCE, MIN_PASSES = 1.0, 4


def perturbation(arm: str, base_q: int = 20) -> tuple[str, float]:
    m = re.search(r"_sx(ext|mu|low|chlb)(\d{3})_", arm)
    if m:
        return m.group(1), int(m.group(2)) / 100
    m = re.search(r"_sxq(\d{3})_", arm)
    if m:
        return "Q", int(m.group(1)) / 100
    m = re.search(r"_sh(\d{3})$", arm)
    if m:
        return "share", int(m.group(1)) / 100
    q = int(re.search(r"_q(\d+)_", arm).group(1))
    return ("Q", q / base_q) if q != base_q else ("base", 1.0)


def rows_for(run: Path, rescore: Path, gate: Path, keep, base_q: int = 20) -> list[dict]:
    traj = pd.read_csv(run / "trajectories.csv")
    scores = pd.read_csv(rescore / "detectability_scores.csv")
    g2 = pd.read_csv(gate / "gate_g2prime.csv").set_index(["treatment", "a_half"])
    out = []
    for arm in sorted(set(traj.treatment)):
        if not keep(arm):
            continue
        s = scores[(scores.treatment == arm) & (scores.a_half == 0.0)]
        p = s[s.score_v2_primary.astype(str).str.lower() == "true"].set_index("reef_name")
        par, f = perturbation(arm, base_q)
        row = {"param": par, "factor": f, "seed": int(re.search(r"_s(\d)_", arm).group(1)), "treatment": arm,
               "loss_v3": float(p.loss_v3.sum()), "loss_v2_manta": float(p.loss_v2_manta.sum()),
               "diver_error": float(p.diver_error.sum()), "g2_pass_ahalf30": bool(g2.loc[(arm, 30.0)].passes_g2prime)}
        for reef, short in SHORT.items():
            a = traj[(traj.treatment == arm) & (traj.reef_name == reef)].set_index("year").adults_ha
            row[f"{short}_2024_adults"] = float(a.get(2024, np.nan))
            row[f"{short}_2026_adults"] = float(a.get(2026, np.nan))
            w3 = a.loc[1993:2003]
            row[f"{short}_w3_peak_year"] = int(w3.idxmax())
        out.append(row)
    return out


def main(run: Path, rescore: Path, gate: Path, base_q: int = 20) -> None:
    paths = {n: run / n for n in ("sens_audit.csv", "sens_summary.csv", "sens_audit.json", "sens.png")}
    if any(p.exists() for p in paths.values()):
        raise FileExistsError("Refusing to overwrite sensitivity audit")
    if len(pd.read_csv(run / "failures.csv")):
        raise ValueError("Failed evaluations")
    rows = rows_for(run, rescore, gate, lambda a: True, base_q)
    if base_q == 20:
        rows += rows_for(*(RUNS / p for p in PITQ), lambda a: re.search(r"_q(15|25)_", a) is not None)
    result = pd.DataFrame(rows).sort_values(["param", "factor", "seed"]).reset_index(drop=True)
    # Base replays the pit_q run on the same layer (q20: 20261016T_02; q25: 20261016T_03).
    replay = {}
    cur = {n: pd.read_csv(run / f"{n}.csv") for n in ("trajectories", "flows", "mechanisms")}
    for name in cur:
        prior = pd.read_csv(RUNS / PRIOR[base_q] / f"{name}.csv")
        cols = [c for c in prior if c not in ("treatment", "reef_name", "year")]
        for arm in result[result.param == "base"].treatment:
            replay[f"{arm}/{name}"] = exact_replay(cur[name][cur[name].treatment == arm], prior[prior.treatment == arm], cols)
    result.to_csv(paths["sens_audit.csv"], index=False)
    g = result.groupby(["param", "factor"])
    summary = g.agg(v3_mean=("loss_v3", "mean"), v3_min=("loss_v3", "min"), v3_max=("loss_v3", "max"),
                    manta_mean=("loss_v2_manta", "mean"), diver_mean=("diver_error", "mean"),
                    g2_passes_ahalf30=("g2_pass_ahalf30", "sum"),
                    **{f"{s}_2024_adults": (f"{s}_2024_adults", "mean") for s in SHORT.values()},
                    **{f"{s}_2026_adults": (f"{s}_2026_adults", "mean") for s in SHORT.values()},
                    **{f"{s}_w3_peak_year": (f"{s}_w3_peak_year", "mean") for s in ("liz", "macg", "nd")}).reset_index()
    base = float(summary[summary.param == "base"].v3_mean.iloc[0])
    summary["d_v3"] = summary.v3_mean - base
    summary["pass"] = (summary.g2_passes_ahalf30 >= MIN_PASSES) & (summary.d_v3 <= V3_TOLERANCE)
    summary.loc[summary.param == "base", "pass"] = summary.loc[summary.param == "base", "g2_passes_ahalf30"] >= MIN_PASSES
    summary.to_csv(paths["sens_summary.csv"], index=False)

    s = summary[summary.param != "base"].copy()
    s["label"] = s.param + " x" + s.factor.map(lambda f: f"{f:g}")
    fig, ax = plt.subplots(figsize=(9, 4.8))
    y = np.arange(len(s))
    ax.barh(y, s.d_v3, color=np.where(s["pass"], "#3987e5", "#eb6834"), height=0.6)
    ax.axvline(0, color="#52514e", linewidth=0.8)
    ax.axvline(V3_TOLERANCE, color="#eb6834", linewidth=0.8, linestyle="--")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{l}  (G2' {int(p)}/5)" for l, p in zip(s.label, s.g2_passes_ahalf30)])
    ax.set_xlabel(f"Change in seed-mean v3 vs base ({base:.2f}); dashed = +{V3_TOLERANCE} limit")
    ax.grid(axis="x", color="#e4e3df", linewidth=.8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(paths["sens.png"], dpi=150, facecolor="#fcfcfb")
    plt.close(fig)
    report = {"status": "promotion_gate_5", "base_q": base_q, "replay_prior": PRIOR[base_q], "base_v3": base, "gate5_pass": bool(summary["pass"].all()),
              "reference_replay_max_rel_delta": replay, "summary": json.loads(summary.to_json(orient="records")),
              "inputs": {**{f"{p.name}/{n}": sha256(p / n) for p, n in ((run, "metadata.toml"), (run, "trajectories.csv"),
                                                                           (rescore, "detectability_scores.csv"), (gate, "gate_g2prime.csv"))},
                         "targets": sha256(TARGETS)},
              "script_sha256": sha256(Path(__file__)),
              "outputs": {n: sha256(p) for n, p in paths.items() if n != "sens_audit.json"}}
    paths["sens_audit.json"].write_text(json.dumps(report, indent=2, default=str) + "\n")
    pd.set_option("display.width", 250)
    print(summary.round(2).to_string(index=False))
    print("Gate 5 pass:", report["gate5_pass"], "| replay max:", max(replay.values()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("run", "rescore", "gate"):
        parser.add_argument(name, type=Path)
    parser.add_argument("--base-q", type=int, default=20, choices=sorted(PRIOR))
    a = parser.parse_args()
    main(a.run.resolve(), a.rescore.resolve(), a.gate.resolve(), a.base_q)
