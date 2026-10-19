"""Audit chl survival on external supply (crash_mechanism_protocol.md, "Preregistered: chl survival on
external supply", 2026-10-19, twenty-fourth round): W-wave x {off, _xchl_}, five seeds.

Replay: off arms must equal the _wv025 arms of 20261019T_15. Failure (preregistered): dv3 > +0.5,
G2' (A_half 30) < 4/5, a 1990s peak moved by > 1 y, or diver error not lower.

    python audit_1991_ext_chl_yr.py <run> <rescore> <gate>
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
TARGETS = Path(__file__).resolve().parents[2] / "sandbox/data/derived/density_targets_2020s/targets_reef_year.csv"
REEFS = ("Lizard Island Reef", "MacGillivray Reef", "North Direction Reef", "Eyrie Reef")
W3_REEFS = REEFS[:3]
WAVES = {"w3": (1993, 2003), "w4": (2008, 2016)}


def arm_on(arm: str) -> bool:
    return "_xchl_" in arm


def main(run: Path, rescore: Path, gate: Path) -> None:
    paths = {n: run / n for n in ("ext_chl_yr_audit.csv", "ext_chl_yr_summary.csv", "ext_chl_yr_audit.json", "ext_chl_yr.png")}
    if any(p.exists() for p in paths.values()):
        raise FileExistsError("Refusing to overwrite audit")
    if len(pd.read_csv(run / "failures.csv")):
        raise ValueError("Failed evaluations")
    traj = pd.read_csv(run / "trajectories.csv")
    flows = pd.read_csv(run / "flows.csv")
    scores = pd.read_csv(rescore / "detectability_scores.csv")
    g2 = pd.read_csv(gate / "gate_g2prime.csv").set_index(["treatment", "a_half"])
    targets = pd.read_csv(TARGETS)
    arms = sorted(set(traj.treatment))

    replay = {}
    for name in ("trajectories", "flows", "mechanisms"):
        cur, prior = pd.read_csv(run / f"{name}.csv"), pd.read_csv(RUNS / PRIOR / f"{name}.csv")
        cols = [c for c in prior if c not in ("treatment", "reef_name", "year")]
        for arm in (a for a in arms if not arm_on(a)):
            replay[f"{arm}/{name}"] = exact_replay(cur[cur.treatment == arm], prior[prior.treatment == arm], cols)

    rows = []
    for arm in arms:
        s = scores[(scores.treatment == arm) & (scores.a_half == 0.0)]
        p = s[s.score_v2_primary.astype(str).str.lower() == "true"].set_index("reef_name")
        row = {"treatment": arm, "on": arm_on(arm), "seed": int(re.search(r"_s(\d)_", arm).group(1)),
               "loss_v3": float(p.loss_v3.sum()), "loss_v2_manta": float(p.loss_v2_manta.sum()),
               "diver_error": float(p.diver_error.sum()), "g2_pass_ahalf30": bool(g2.loc[(arm, 30.0)].passes_g2prime)}
        for reef in REEFS:
            key = reef.split()[0].lower()
            row[f"{key}_v3"] = float(p.loss_v3.get(reef, np.nan))
            row[f"{key}_diver"] = float(p.diver_error.get(reef, np.nan))
            a = traj[(traj.treatment == arm) & (traj.reef_name == reef)].set_index("year").adults_ha
            for w, (y0, y1) in WAVES.items():
                row[f"{key}_{w}_peak_year"] = int(a.loc[y0:y1].idxmax())
                row[f"{key}_{w}_peak_adults"] = float(a.loc[y0:y1].max())
            for y in range(2019, 2027):
                row[f"{key}_{y}"] = float(a.get(y, np.nan))
            f = flows[(flows.treatment == arm) & (flows.reef_name == reef)].set_index("year")
            row[f"{key}_ext_2014_18"] = float(f.external_immigration.loc[2014:2018].mean())
            row[f"{key}_ext_2019"] = float(f.external_immigration.get(2019, np.nan))
        rows.append(row)
    d = pd.DataFrame(rows)
    d.to_csv(paths["ext_chl_yr_audit.csv"], index=False)
    num = [c for c in d if c not in ("treatment", "on", "seed")]
    summary = d.groupby("on")[num].mean()
    summary["g2_passes_ahalf30"] = d.groupby("on").g2_pass_ahalf30.sum()
    off, on = summary.loc[False], summary.loc[True]
    shifts = {r: float(on[f"{r.split()[0].lower()}_w3_peak_year"] - off[f"{r.split()[0].lower()}_w3_peak_year"]) for r in W3_REEFS}
    verdict = {"d_v3": float(on.loss_v3 - off.loss_v3), "d_diver_error": float(on.diver_error - off.diver_error),
               "d_manta": float(on.loss_v2_manta - off.loss_v2_manta), "g2_passes": int(on.g2_passes_ahalf30),
               "w3_peak_shift_years": shifts}
    verdict["pass"] = bool(verdict["d_v3"] <= 0.5 and verdict["g2_passes"] >= 4 and
                           max(abs(v) for v in shifts.values()) <= 1.0 and verdict["d_diver_error"] < 0)
    summary = summary.reset_index()
    summary.to_csv(paths["ext_chl_yr_summary.csv"], index=False)

    fig, axes = plt.subplots(1, 4, figsize=(17, 4.2), sharex=True)
    for ax, reef in zip(axes, REEFS):
        for flag, color, label in ((False, "#0d366b", "W-wave"), (True, "#eb6834", "+ chl on external supply")):
            sel = [a for a in arms if arm_on(a) == flag]
            a = traj[traj.treatment.isin(sel) & (traj.reef_name == reef)].pivot(index="year", columns="treatment", values="adults_ha")
            ax.plot(a.index, a.mean(axis=1), color=color, linewidth=1.8, label=label)
        t = targets[(targets.target_reef == reef) & (targets.target_density_ha > 0)]
        ax.scatter(t.report_year, t.target_density_ha, s=16, color="#52514e", zorder=5, label="2020s diver/cull target")
        ax.set_yscale("log")
        ax.set_title(reef, loc="left", fontsize=10)
        ax.set_xlim(1991, 2027)
        ax.grid(color="#e4e3df", linewidth=.8)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Adults ha$^{-1}$ (log)")
    axes[0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(paths["ext_chl_yr.png"], dpi=150, facecolor="#fcfcfb")
    plt.close(fig)

    report = {"status": "mechanism_screen", "verdict": verdict, "replay_max_rel_delta": max(replay.values()),
              "summary": json.loads(summary.to_json(orient="records")),
              "inputs": {f"{p.name}/{n}": sha256(p / n) for p, n in ((run, "metadata.toml"), (run, "trajectories.csv"),
                                                                       (rescore, "detectability_scores.csv"), (gate, "gate_g2prime.csv"))},
              "script_sha256": sha256(Path(__file__)),
              "outputs": {n: sha256(p) for n, p in paths.items() if n != "ext_chl_yr_audit.json"}}
    paths["ext_chl_yr_audit.json"].write_text(json.dumps(report, indent=2, default=str) + "\n")
    pd.set_option("display.width", 250)
    print(summary[["on", "loss_v3", "loss_v2_manta", "diver_error", "g2_passes_ahalf30"] +
                  [f"{r}_{c}" for r in ("lizard", "macgillivray", "north", "eyrie") for c in ("v3", "diver")]].round(2).to_string(index=False))
    print(summary[["on"] + [c for c in summary if re.search(r"_w[34]_peak_", c)]].round(1).to_string(index=False))
    print(summary[["on"] + [c for c in summary if re.search(r"_(2020|2022|2024|2026|ext_2014_18|ext_2019)$", c)]].round(1).to_string(index=False))
    print(json.dumps(verdict, indent=1))
    print("Replay max:", report["replay_max_rel_delta"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("run", "rescore", "gate"):
        parser.add_argument(name, type=Path)
    a = parser.parse_args()
    main(a.run.resolve(), a.rescore.resolve(), a.gate.resolve())
