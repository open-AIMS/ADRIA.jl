"""Audit the Lizard 1991 fast-coral test (crash_mechanism_protocol.md, "Lead's request (2026-10-19,
twenty-eighth round) and preregistered Lizard 1991 fast-coral test"): W-wave x Lizard fast increment
{0, +0.02, +0.04}, five seeds.

Replay: 0 arms must equal the _wv025 arms of 20261019T_15. "Fit better": seed-mean v3 below W-wave with
G2' (A_half 30) >= 4/5. Also reported: Lizard's own v3 term, peak height against Pratchett (40-55 adults/ha)
and manta (1.92 COTS/tow), fast cover against photo transects, troughs and neighbouring reefs.

    python audit_1991_lizfast.py <run> <rescore> <gate>
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
PRIOR = "20261019T_15_wave_w3_q25_1991_v1"
DATA = Path(__file__).resolve().parents[2] / "sandbox" / "data"
REEFS = ("Lizard Island Reef", "MacGillivray Reef", "North Direction Reef", "Eyrie Reef")
Q = 0.015
COLORS = {0.0: "#0d366b", 0.02: "#3987e5", 0.04: "#eb6834"}


def delta(arm: str) -> float:
    m = re.search(r"_lf(\d{2})_", arm)
    return 0.0 if m is None else int(m.group(1)) / 100


def main(run: Path, rescore: Path, gate: Path) -> None:
    paths = {n: run / n for n in ("lizfast_audit.csv", "lizfast_summary.csv", "lizfast_audit.json", "lizfast.png")}
    if any(p.exists() for p in paths.values()):
        raise FileExistsError("Refusing to overwrite audit")
    if len(pd.read_csv(run / "failures.csv")):
        raise ValueError("Failed evaluations")
    traj, mech = pd.read_csv(run / "trajectories.csv"), pd.read_csv(run / "mechanisms.csv")
    scores = pd.read_csv(rescore / "detectability_scores.csv")
    g2 = pd.read_csv(gate / "gate_g2prime.csv").set_index(["treatment", "a_half"])
    arms = sorted(set(traj.treatment))
    replay = {}
    for name in ("trajectories", "flows", "mechanisms"):
        cur, prior = pd.read_csv(run / f"{name}.csv"), pd.read_csv(RUNS / PRIOR / f"{name}.csv")
        cols = [c for c in prior if c not in ("treatment", "reef_name", "year")]
        for arm in (a for a in arms if delta(a) == 0.0):
            replay[f"{arm}/{name}"] = exact_replay(cur[cur.treatment == arm], prior[prior.treatment == arm], cols)
    rows = []
    for arm in arms:
        s = scores[(scores.treatment == arm) & (scores.a_half == 0.0)]
        p = s[s.score_v2_primary.astype(str).str.lower() == "true"].set_index("reef_name")
        row = {"treatment": arm, "delta": delta(arm), "seed": int(re.search(r"_s(\d)_", arm).group(1)),
               "loss_v3": float(p.loss_v3.sum()), "loss_v2_manta": float(p.loss_v2_manta.sum()),
               "diver_error": float(p.diver_error.sum()), "g2_pass_ahalf30": bool(g2.loc[(arm, 30.0)].passes_g2prime)}
        for reef in REEFS:
            key = reef.split()[0].lower()
            row[f"{key}_v3"] = float(p.loss_v3.get(reef, np.nan))
            row[f"{key}_manta"] = float(p.loss_v2_manta.get(reef, np.nan))
            a = traj[(traj.treatment == arm) & (traj.reef_name == reef)].set_index("year")
            w3 = a.adults_ha.loc[1993:2003]
            row[f"{key}_w3_peak_year"], row[f"{key}_w3_peak_adults"] = int(w3.idxmax()), float(w3.max())
            w4 = a.adults_ha.loc[2008:2016]
            row[f"{key}_w4_peak_adults"] = float(w4.max())
            row[f"{key}_trough_cpue_2000_2010"] = float((Q * a.adults_ha.loc[2000:2010]).mean())
            row[f"{key}_coral_1999"] = float(a.coral_cover.get(1999, np.nan))
        m = mech[(mech.treatment == arm) & (mech.reef_name == "Lizard Island Reef")].set_index("year")
        for y in (1991, 1993, 1995, 1997):
            row[f"lizard_fast_{y}"] = float(m.fast_coral.get(y, np.nan))
        row["lizard_min_food_survival_1995_99"] = float(m.food_survival.loc[1995:1999].min())
        row["lizard_max_condition_1993_99"] = float(m.body_condition.loc[1993:1999].max())
        rows.append(row)
    d = pd.DataFrame(rows)
    d.to_csv(paths["lizfast_audit.csv"], index=False)
    num = [c for c in d if c not in ("treatment", "delta", "seed")]
    summary = d.groupby("delta")[num].mean()
    summary["g2_passes_ahalf30"] = d.groupby("delta").g2_pass_ahalf30.sum()
    summary["d_v3"] = summary.loss_v3 - summary.loc[0.0, "loss_v3"]
    summary["fits_better"] = (summary.d_v3 < 0) & (summary.g2_passes_ahalf30 >= 4)
    summary = summary.reset_index()
    summary.to_csv(paths["lizfast_summary.csv"], index=False)

    photo = pd.read_csv(DATA / "reef_photo_transect.csv", low_memory=False)
    ph = photo[(photo.variable == "HARD CORAL") & (photo.project_code == "LTMP") & (photo.domain_name == "Lizard Isles")]
    fast_obs = ph[(ph.purpose == "COMPOSITION") & ph.reefpage_category.isin(("Acropora", "Pocilloporidae"))].groupby("report_year")["median"].sum()
    obs = pd.read_csv(DATA / "reef_cots.csv")
    o = obs[obs.reef_name == "Lizard Isles"].groupby("year").cotsptow.mean()
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.3))
    for dl, color in COLORS.items():
        sel = [a for a in arms if delta(a) == dl]
        lab = "W-wave" if dl == 0 else f"Lizard 1991 fast +{dl:.2f}"
        a = traj[traj.treatment.isin(sel) & (traj.reef_name == "Lizard Island Reef")].pivot(index="year", columns="treatment", values="adults_ha")
        axes[0].plot(a.index, a.mean(axis=1), color=color, linewidth=1.9, label=lab)
        axes[0].fill_between(a.index, a.min(axis=1), a.max(axis=1), color=color, alpha=0.15, linewidth=0)
        m = mech[mech.treatment.isin(sel) & (mech.reef_name == "Lizard Island Reef")].pivot(index="year", columns="treatment", values="fast_coral")
        axes[1].plot(m.index, m.mean(axis=1), color=color, linewidth=1.9, label=lab)
        c = traj[traj.treatment.isin(sel) & (traj.reef_name == "Lizard Island Reef")].pivot(index="year", columns="treatment", values="coral_cover")
        axes[2].plot(c.index, c.mean(axis=1), color=color, linewidth=1.9, label=lab)
    axes[0].scatter(o.index, o.values / Q, s=14, color="#1f1e1d", zorder=5, label="LTMP manta ÷ 0.015")
    axes[0].axhspan(40, 55, color="#b6b5b0", alpha=0.35, linewidth=0, label="Pratchett 2005 diver density")
    axes[0].set_ylabel("Lizard adults ha$^{-1}$")
    axes[0].set_ylim(0, 140)
    axes[1].scatter(fast_obs.index, fast_obs.values, s=16, color="#1f1e1d", zorder=5, label="LTMP 9 m photo fast cover")
    axes[1].set_ylabel("Fast coral (fraction)")
    axes[2].set_ylabel("Total coral (fraction of habitable area)")
    axes[2].scatter([1996, 1999], [0.322, 0.229], s=30, marker="D", color="#eb6834", zorder=5, label="Pratchett 2005 reef-wide")
    for ax, title in zip(axes, ("Lizard adults", "Lizard fast coral", "Lizard total coral")):
        ax.set_title(title, loc="left", fontsize=10)
        ax.set_xlim(1991, 2027)
        ax.grid(color="#e4e3df", linewidth=.8)
        ax.spines[["top", "right"]].set_visible(False)
        ax.legend(frameon=False, fontsize=7.5)
    fig.tight_layout()
    fig.savefig(paths["lizfast.png"], dpi=150, facecolor="#fcfcfb")
    plt.close(fig)

    report = {"status": "lead_request_screen", "replay_max_rel_delta": max(replay.values()),
              "summary": json.loads(summary.to_json(orient="records")),
              "inputs": {f"{p.name}/{n}": sha256(p / n) for p, n in ((run, "metadata.toml"), (run, "trajectories.csv"),
                                                                       (rescore, "detectability_scores.csv"), (gate, "gate_g2prime.csv"))},
              "script_sha256": sha256(Path(__file__)),
              "outputs": {n: sha256(p) for n, p in paths.items() if n != "lizfast_audit.json"}}
    paths["lizfast_audit.json"].write_text(json.dumps(report, indent=2, default=str) + "\n")
    pd.set_option("display.width", 260)
    print(summary[["delta", "loss_v3", "d_v3", "loss_v2_manta", "diver_error", "g2_passes_ahalf30", "fits_better"] +
                  [f"{r}_v3" for r in ("lizard", "macgillivray", "north", "eyrie")]].round(3).to_string(index=False))
    print(summary[["delta"] + [c for c in summary if c.startswith("lizard_")]].round(3).T.to_string())
    print(summary[["delta"] + [f"{r}_w3_peak_adults" for r in ("macgillivray", "north", "eyrie")] +
                  [f"{r}_trough_cpue_2000_2010" for r in ("lizard", "macgillivray")]].round(3).to_string(index=False))
    print("Replay max:", report["replay_max_rel_delta"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("run", "rescore", "gate"):
        parser.add_argument(name, type=Path)
    a = parser.parse_args()
    main(a.run.resolve(), a.rescore.resolve(), a.gate.resolve())
