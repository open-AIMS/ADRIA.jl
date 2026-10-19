"""Audit the graded outbreak-layer correction (crash_mechanism_protocol.md, "Preregistered: graded
outbreak-layer correction", 2026-10-19, twenty-fifth round): W-wave on owen_global_2018_2023_corr_w3_graded
against W-wave on corr_w3 (the _wv025 arms of run 20261019T_15, re-score 16, G2' 17), five seeds.

- Identity: 1991-2017 trajectories, flows and mechanisms must equal the control exactly (the layers
  differ only from 2018).
- Failure (preregistered): diver error not lower, dv3 > +0.5, or G2' (A_half 30) < 4/5.

    python audit_1991_graded.py <run> <rescore> <gate>
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
CONTROL = ("20261019T_15_wave_w3_q25_1991_v1", "20261019T_16_rescore_wave_w3_q25", "20261019T_17_gate_wave_w3_q25_lizmacg")
TARGETS = Path(__file__).resolve().parents[2] / "sandbox/data/derived/density_targets_2020s/targets_reef_year.csv"
REEFS = ("Lizard Island Reef", "MacGillivray Reef", "North Direction Reef", "Eyrie Reef")
RIBBONS = ("Yonge Reef", "Carter Reef", "No Name Reef")
Q, A_HALF = 0.015, 30.0


def load(run: Path, rescore: Path, gate: Path, label: str) -> tuple[list[dict], dict]:
    traj = pd.read_csv(run / "trajectories.csv")
    traj = traj[traj.treatment.str.endswith("_wv025")]
    flows = pd.read_csv(run / "flows.csv")
    extra = pd.read_csv(run / "trajectories_extra.csv")
    scores = pd.read_csv(rescore / "detectability_scores.csv")
    g2 = pd.read_csv(gate / "gate_g2prime.csv").set_index(["treatment", "a_half"])
    rows = []
    for arm in sorted(set(traj.treatment)):
        s = scores[(scores.treatment == arm) & (scores.a_half == 0.0)]
        p = s[s.score_v2_primary.astype(str).str.lower() == "true"].set_index("reef_name")
        row = {"layer": label, "treatment": arm, "seed": int(re.search(r"_s(\d)_", arm).group(1)),
               "loss_v3": float(p.loss_v3.sum()), "loss_v2_manta": float(p.loss_v2_manta.sum()),
               "diver_error": float(p.diver_error.sum()), "g2_pass_ahalf30": bool(g2.loc[(arm, 30.0)].passes_g2prime)}
        for reef in REEFS:
            key = reef.split()[0].lower()
            row[f"{key}_v3"] = float(p.loss_v3.get(reef, np.nan))
            row[f"{key}_diver"] = float(p.diver_error.get(reef, np.nan))
            a = traj[(traj.treatment == arm) & (traj.reef_name == reef)].set_index("year").adults_ha
            for y in range(2018, 2027):
                row[f"{key}_{y}"] = float(a.get(y, np.nan))
            hit = a.loc[2019:][a.loc[2019:] >= 5.0]
            row[f"{key}_first_ge5"] = int(hit.index.min()) if len(hit) else 2027
            f = flows[(flows.treatment == arm) & (flows.reef_name == reef)].set_index("year")
            for y in (2018, 2019, 2020, 2021, 2022):
                row[f"{key}_ext_{y}"] = float(f.external_immigration.get(y, np.nan))
        for reef in RIBBONS:
            a = extra[(extra.treatment == arm) & (extra.reef_name == reef)].set_index("year").adults_ha.loc[2018:]
            row[f"{reef.split()[0].lower()}_max_cpue_2018_26"] = float((Q * a ** 2 / (a + A_HALF)).max())
        rows.append(row)
    return rows, {"traj": traj}


def main(run: Path, rescore: Path, gate: Path) -> None:
    paths = {n: run / n for n in ("graded_audit.csv", "graded_summary.csv", "graded_audit.json", "graded.png")}
    if any(p.exists() for p in paths.values()):
        raise FileExistsError("Refusing to overwrite graded audit")
    if len(pd.read_csv(run / "failures.csv")):
        raise ValueError("Failed evaluations")
    ctrl = [RUNS / p for p in CONTROL]
    identity = {}
    for name in ("trajectories", "flows", "mechanisms"):
        cur, prior = pd.read_csv(run / f"{name}.csv"), pd.read_csv(ctrl[0] / f"{name}.csv")
        cols = [c for c in prior if c not in ("treatment", "reef_name", "year")]
        for arm in sorted(set(cur.treatment)):
            identity[f"{arm}/{name}"] = exact_replay(cur[(cur.treatment == arm) & (cur.year <= 2017)],
                                                     prior[(prior.treatment == arm) & (prior.year <= 2017)], cols)
    rows_g, tg = load(run, rescore, gate, "graded")
    rows_c, tc = load(*ctrl, "binary")
    d = pd.DataFrame(rows_c + rows_g)
    d.to_csv(paths["graded_audit.csv"], index=False)
    num = [c for c in d if c not in ("layer", "treatment", "seed")]
    summary = d.groupby("layer")[num].mean()
    summary["g2_passes_ahalf30"] = d.groupby("layer").g2_pass_ahalf30.sum()
    b, g = summary.loc["binary"], summary.loc["graded"]
    verdict = {"d_v3": float(g.loss_v3 - b.loss_v3), "d_manta": float(g.loss_v2_manta - b.loss_v2_manta),
               "d_diver_error": float(g.diver_error - b.diver_error), "g2_passes": int(g.g2_passes_ahalf30),
               "first_year_ge5": {r: (float(b[f"{r.split()[0].lower()}_first_ge5"]), float(g[f"{r.split()[0].lower()}_first_ge5"])) for r in REEFS}}
    verdict["pass"] = bool(verdict["d_diver_error"] < 0 and verdict["d_v3"] <= 0.5 and verdict["g2_passes"] >= 4)
    summary = summary.reset_index()
    summary.to_csv(paths["graded_summary.csv"], index=False)

    targets = pd.read_csv(TARGETS)
    fig, axes = plt.subplots(1, 4, figsize=(17, 4.2), sharex=True)
    for ax, reef in zip(axes, REEFS):
        for t, color, label in ((tc["traj"], "#0d366b", "binary correction (W-wave)"), (tg["traj"], "#eb6834", "graded correction")):
            a = t[t.reef_name == reef].pivot(index="year", columns="treatment", values="adults_ha")
            ax.plot(a.index, a.mean(axis=1), color=color, linewidth=1.8, label=label)
        tt = targets[(targets.target_reef == reef) & (targets.target_density_ha > 0)]
        ax.scatter(tt.report_year, tt.target_density_ha, s=16, color="#52514e", zorder=5, label="2020s diver/cull target")
        ax.set_yscale("log")
        ax.set_title(reef, loc="left", fontsize=10)
        ax.set_xlim(2010, 2027)
        ax.grid(color="#e4e3df", linewidth=.8)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Adults ha$^{-1}$ (log)")
    axes[0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(paths["graded.png"], dpi=150, facecolor="#fcfcfb")
    plt.close(fig)

    report = {"status": "layer_screen", "verdict": verdict, "identity_1991_2017_max_rel_delta": max(identity.values()),
              "summary": json.loads(summary.to_json(orient="records")),
              "inputs": {f"{p.name}/{n}": sha256(p / n) for p, n in ((run, "metadata.toml"), (run, "trajectories.csv"),
                                                                       (rescore, "detectability_scores.csv"), (gate, "gate_g2prime.csv"),
                                                                       (ctrl[0], "trajectories.csv"))},
              "script_sha256": sha256(Path(__file__)),
              "outputs": {n: sha256(p) for n, p in paths.items() if n != "graded_audit.json"}}
    paths["graded_audit.json"].write_text(json.dumps(report, indent=2, default=str) + "\n")
    pd.set_option("display.width", 260)
    print(summary[["layer", "loss_v3", "loss_v2_manta", "diver_error", "g2_passes_ahalf30"] +
                  [f"{r}_{c}" for r in ("lizard", "macgillivray", "north", "eyrie") for c in ("v3", "diver", "first_ge5")]].round(2).to_string(index=False))
    for r in ("lizard", "eyrie", "macgillivray", "north"):
        print(summary[["layer"] + [f"{r}_{y}" for y in range(2018, 2027)] + [f"{r}_ext_{y}" for y in (2018, 2019, 2020, 2021, 2022)]].round(1).to_string(index=False))
    print(summary[["layer"] + [c for c in summary if c.endswith("_max_cpue_2018_26")]].round(3).to_string(index=False))
    print(json.dumps(verdict, indent=1))
    print("Identity 1991-2017 max:", report["identity_1991_2017_max_rel_delta"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("run", "rescore", "gate"):
        parser.add_argument(name, type=Path)
    a = parser.parse_args()
    main(a.run.resolve(), a.rescore.resolve(), a.gate.resolve())
