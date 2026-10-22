"""Compare outside-reef cull scenarios for the 2020s build-up study (crash_mechanism_protocol.md, 2026-10-22,
thirty-third round). Uses the per-scenario audits (audit_buildup_2020s.py) and the run trajectories.

Scenarios: none (no outside-cull discount; runs 20261022T_01-03), B (observed outside culls to 2025), K (B + recent
outside effort carried on 2026-30), D (outside effort doubled 2020-30), L (observed to 2025, doubled 2026-30).

Figures (same layouts as the build-up study):
- outside1_trajectories_2018_2030: rows = outbreak layer, columns = reef; one line per outside scenario, no extra
  inside effort (observed inside culls only); seed mean.
- outside2_success_vs_inside_effort: rows = early / late inside effort (targeted), columns = outbreak layer;
  share of seeds with no breach of 14.7/ha against the inside effort multiple, one line per outside scenario.
Table: outside_cull_comparison.csv (required inside effort and peak adults per scenario).

    python compare_outside_cull_2020s.py <report_dir> --audits none=<dir> B=<dir> K=<dir> D=<dir> L=<dir>
        --runs none=<C>,<G>,<O> B=<C>,<G>,<O> ...
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REEFS = ["Lizard Island Reef", "MacGillivray Reef", "North Direction Reef", "Eyrie Reef", "Martin Reef"]
CONFIG_NAME = {"C": "current (binary correction)", "G": "graded correction 0.5", "O": "original layer (no correction)"}
SCEN = {"none": "no outside-cull discount", "B": "outside culls observed to 2025",
        "K": "outside: recent effort carried on 2026-30", "D": "outside: effort doubled 2020-30",
        "L": "outside: effort doubled 2026-30"}
COLORS = {"none": "#2a78d6", "B": "#eb6834", "K": "#1baf7a", "D": "#eda100", "L": "#e87ba4"}
INK, MUTED, GRID = "#1f1e1d", "#6f6e69", "#e4e3df"
OUTBREAK = 0.22 / 0.015


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def style(ax):
    ax.grid(color=GRID, linewidth=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=8)


def main(out: Path, audits: dict[str, Path], runs: dict[str, dict[str, Path]]) -> None:
    if out.exists():
        raise FileExistsError(f"Refusing to overwrite {out}")
    out.mkdir(parents=True)
    plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold", "text.color": INK,
                         "axes.labelcolor": INK})
    files = {}

    def save(fig, stem):
        for ext in ("png", "pdf"):
            path = out / f"{stem}.{ext}"
            fig.savefig(path, dpi=200, facecolor="#ffffff", bbox_inches="tight")
            files[path.name] = sha256(path)
        plt.close(fig)

    # Control trajectories (no extra inside effort) per scenario.
    traj = []
    for sc, cfgs in runs.items():
        for cfg, run in cfgs.items():
            t = pd.concat([pd.read_csv(run / "trajectories.csv"), pd.read_csv(run / "trajectories_extra.csv")])
            t = t[t.reef_name.isin(REEFS) & t.year.between(2018, 2030) & ~t.treatment.str.contains("_xe")]
            traj.append(t.assign(scenario=sc, config=cfg))
    traj = pd.concat(traj, ignore_index=True)
    fig, axes = plt.subplots(3, len(REEFS), figsize=(18, 9.5), sharex=True)
    for i, cfg in enumerate(CONFIG_NAME):
        for j, reef in enumerate(REEFS):
            ax = axes[i, j]
            for sc, color in COLORS.items():
                g = traj[(traj.config == cfg) & (traj.reef_name == reef) & (traj.scenario == sc)].groupby("year").adults_ha
                ax.plot(g.mean().index, g.mean().values, color=color, linewidth=1.7, label=SCEN[sc])
            ax.axhline(OUTBREAK, color=MUTED, linewidth=0.9, linestyle="--")
            ax.axhline(3.0, color=MUTED, linewidth=0.9, linestyle=":")
            ax.axvspan(2019.5, 2025.5, color="#f1d9c8", alpha=0.25, linewidth=0)
            ax.axvspan(2025.5, 2030.5, color="#cdeee0", alpha=0.25, linewidth=0)
            style(ax)
            if i == 0:
                ax.set_title(reef, loc="left")
            if j == 0:
                ax.set_ylabel(f"{CONFIG_NAME[cfg]}\nadults ha$^{{-1}}$")
    h, lab = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, lab, loc="lower center", ncol=5, frameon=False, fontsize=9)
    fig.suptitle("2020s build-up under outside-reef culling scenarios, no extra culling inside the domain (seed mean; "
                 "dashed = 14.7/ha outbreak, dotted = 3/ha)", x=0.01, ha="left", fontsize=11)
    fig.tight_layout(rect=(0, 0.04, 1, 0.96))
    save(fig, "outside1_trajectories_2018_2030")

    # Success against inside effort (targeted), per scenario.
    metrics = pd.concat([pd.read_csv(a / "buildup_arm_metrics.csv").assign(scenario=sc) for sc, a in audits.items()],
                        ignore_index=True)
    mults = sorted(metrics.mult.unique())
    fig, axes = plt.subplots(2, 3, figsize=(16, 7.5), sharex=True, sharey=True)
    for i, (wkey, wname) in enumerate((("E", "inside effort 2020-25 (judged 2020-26)"), ("L", "inside effort 2026-30 (judged 2027-30)"))):
        for j, cfg in enumerate(CONFIG_NAME):
            ax = axes[i, j]
            for sc, color in COLORS.items():
                d = metrics[(metrics.scenario == sc) & (metrics.config == cfg) & (metrics.window == wkey) &
                            metrics["mode"].isin(["none", "T"])]
                f = d.groupby("mult")[f"success_{OUTBREAK:g}"].mean().reindex(mults)
                x = np.array([0.25 if v == 0 else v for v in f.index])
                ax.plot(x, f.values, marker="o", markersize=4, color=color, linewidth=1.6, label=SCEN[sc])
            ax.set_xscale("log", base=2)
            ax.set_xticks([0.25] + [v for v in mults if v > 0], ["0"] + [f"{v:g}" for v in mults if v > 0])
            ax.axhline(0.8, color=MUTED, linewidth=0.8, linestyle=":")
            style(ax)
            if i == 0:
                ax.set_title(CONFIG_NAME[cfg], loc="left")
            if j == 0:
                ax.set_ylabel(f"Share of seeds with no breach of 14.7/ha\n{wname}")
            if i == 1:
                ax.set_xlabel("Extra targeted effort inside the domain, multiple of 3,733 diver-hours/yr")
    h, lab = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, lab, loc="lower center", ncol=5, frameon=False, fontsize=9)
    fig.suptitle("Does culling outside the domain reduce the inside effort needed? (five reefs; >= 4 of 5 seeds = 0.8 line)",
                 x=0.01, ha="left", fontsize=11)
    fig.tight_layout(rect=(0, 0.07, 1, 0.95))
    save(fig, "outside2_success_vs_inside_effort")

    # Table.
    rows = []
    for sc, a in audits.items():
        req = pd.read_csv(a / "required_effort.csv")
        req = req[req.threshold.str.startswith("14.7")]
        ctrl = traj[traj.scenario == sc]
        for cfg in CONFIG_NAME:
            c = ctrl[ctrl.config == cfg]
            peak = lambda y0, y1: float(c[c.year.between(y0, y1)].groupby("treatment").adults_ha.max().median())
            row = {"scenario": SCEN[sc], "config": CONFIG_NAME[cfg], "control_max_adults_2020_25": peak(2020, 2025),
                   "control_max_adults_2026": peak(2026, 2026), "control_max_adults_2027_30": peak(2027, 2030)}
            for r in req[req.config == cfg].itertuples():
                key = ("early" if r.window.startswith("early") else "late") + ("_targeted" if "targeted" in r.deployment else "_recent")
                row[f"required_inside_multiplier_{key}"] = r.required_multiplier
            rows.append(row)
    table = pd.DataFrame(rows)
    table.to_csv(out / "outside_cull_comparison.csv", index=False)
    files["outside_cull_comparison.csv"] = sha256(out / "outside_cull_comparison.csv")
    (out / "metadata.json").write_text(json.dumps({
        "status": "scenario_study_not_promoted", "protocol": "crash_mechanism_protocol.md, 2026-10-22, thirty-third round",
        "scenarios": SCEN, "audits": {k: str(v) for k, v in audits.items()},
        "runs": {k: {c: r.name for c, r in v.items()} for k, v in runs.items()},
        "script_sha256": sha256(Path(__file__)), "outputs": files}, indent=2) + "\n")
    pd.set_option("display.width", 250)
    print(table.round(1).to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out", type=Path)
    parser.add_argument("--audits", nargs="+", required=True)
    parser.add_argument("--runs", nargs="+", required=True)
    a = parser.parse_args()
    audits = {k: Path(v).resolve() for k, v in (x.split("=", 1) for x in a.audits)}
    runs = {k: dict(zip(("C", "G", "O"), (Path(p).resolve() for p in v.split(",")))) for k, v in (x.split("=", 1) for x in a.runs)}
    main(a.out.resolve(), audits, runs)
