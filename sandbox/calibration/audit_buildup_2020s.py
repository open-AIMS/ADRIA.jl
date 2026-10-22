"""Audit and report figures for the 2020s build-up culling study (crash_mechanism_protocol.md, "Lead's request
(2026-10-22, thirty-second round) and preregistered 2020s build-up culling study").

Three outbreak layers (C current corr_w3, G graded 0.5, O original), observed culls in every arm, extra forward
effort in an early window (E, 2020-2025) or a late window (L, 2026-2030), deployed as each reef's recent effort x m
(F) or the domain total x m to reefs at >= 3 adults/ha (T), m in {0.5 ... 32}; 5 seeds; years to 2030.

Checks: no failures; C controls equal the culling control (run 20261021T_01) over 1991-2026; E arms equal their
control before 2020 and L arms before 2026; deployed effort is zero outside each window.
Success: no breach (adults >= threshold) on Lizard, MacGillivray, North Direction, Eyrie and Martin in the window
(E: 2020-2026; L: 2027-2030), thresholds 14.7/ha (0.22 COTS/tow) and 3/ha. Required effort = smallest m with >= 4
of 5 seeds succeeding.

    python audit_buildup_2020s.py <out_dir> --C <run> --G <run> --O <run> [--scenario-label "..."] [--control-replay-before 2002]

Outside-reef cull scenarios (2026-10-22, thirty-third round): runs with OWEN_1991_PROJ_OUTSIDE_CULL carry _oc<code>
in arm names; their supply differs from the culling control from the first discounted year, so the C-control replay
is checked only before --control-replay-before (default 2027). --scenario-label is added to figure titles.
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

HERE = Path(__file__).resolve().parent
CONTROL = HERE / "runs" / "20261021T_01_lizfast_cull_control_1991_v1"
REEFS = ["Lizard Island Reef", "MacGillivray Reef", "North Direction Reef", "Eyrie Reef", "Martin Reef"]
THRESHOLDS = {"14.7/ha (0.22 COTS/tow)": 0.22 / 0.015, "3/ha": 3.0}
WINDOWS = {"E": ((2020, 2025), (2020, 2026), "early: extra effort 2020-25"),
           "L": ((2026, 2030), (2027, 2030), "late: extra effort 2026-30")}
CONFIG_NAME = {"C": "current (binary correction)", "G": "graded correction 0.5", "O": "original layer (no correction)"}
LINE = {("E", "F"): "#2a78d6", ("E", "T"): "#eb6834", ("L", "F"): "#1baf7a", ("L", "T"): "#eda100"}
MODE_NAME = {"F": "recent pattern x m", "T": "targeted x m"}
TRAJ_COLORS = {"control": "#2a78d6", "E": "#eb6834", "L": "#1baf7a"}
INK, MUTED, GRID = "#1f1e1d", "#6f6e69", "#e4e3df"
RECENT_HOURS = 3732.79  # domain diver-hours per year, 2023-25 (build_cull_effort.py metadata)
SHOW_M = 8.0


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parts(arm: str) -> dict:
    m = re.search(r"_xe(E|L)(F|T)(\d{3})_", arm)
    seed = int(re.search(r"_s(\d+)_", arm).group(1))
    return {"seed": seed, "window": "control" if m is None else m.group(1),
            "mode": None if m is None else m.group(2), "mult": 0.0 if m is None else int(m.group(3)) / 10}


def style(ax):
    ax.grid(color=GRID, linewidth=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=8)


def max_diff(a, b, cols):
    a, b = a.sort_values(["reef_name", "year"]), b.sort_values(["reef_name", "year"])
    if len(a) != len(b) or len(a) == 0:
        raise ValueError("Replay rows differ")
    return float(np.nanmax(np.abs(a[cols].to_numpy(float) - b[cols].to_numpy(float))))


def main(out: Path, runs: dict[str, Path], scenario_label: str = "", control_before: int = 2027) -> None:
    if out.exists():
        raise FileExistsError(f"Refusing to overwrite {out}")
    checks, rows, traj_all = {"failures": {}, "replay": {}, "effort_outside_window": {}}, [], []
    ctrl = pd.read_csv(CONTROL / "trajectories.csv")
    for cfg, run in runs.items():
        checks["failures"][cfg] = int(len(pd.read_csv(run / "failures.csv")))
        t = pd.concat([pd.read_csv(run / "trajectories.csv"), pd.read_csv(run / "trajectories_extra.csv")], ignore_index=True)
        t = t[t.reef_name.isin(REEFS)]
        main_t = pd.read_csv(run / "trajectories.csv")
        cols = [c for c in main_t if c not in ("treatment", "reef_name", "year")]
        plog = pd.read_csv(run / "projection_log.csv")
        culls = pd.read_csv(run / "culls.csv")
        culls = culls[culls.scope == "domain"]
        arms = sorted(t.treatment.unique())
        meta = {a: parts(a) for a in arms}
        controls = {meta[a]["seed"]: a for a in arms if meta[a]["window"] == "control"}
        rep = 0.0
        for a, p in meta.items():
            c = controls[p["seed"]]
            x, y = main_t[main_t.treatment == a], main_t[main_t.treatment == c]
            cut = {"control": None, "E": 2019, "L": 2025}[p["window"]]
            if cut is not None:
                rep = max(rep, max_diff(x[x.year <= cut], y[y.year <= cut], cols))
            if p["window"] != "control":
                (w0, w1), _, _ = WINDOWS[p["window"]]
                h = plog[plog.treatment == a].set_index("year").effort_hours
                checks["effort_outside_window"][a] = float(h[(h.index < w0) | (h.index > w1)].abs().max())
        if cfg == "C":
            for s, a in controls.items():
                base = re.sub(r"_pj\d{3}g\d{2}W(_oc[BKDL])?_cullreef_xc_aprmak_", "_cullreef_aprmak_", a)
                x = main_t[(main_t.treatment == a) & (main_t.year < min(2027, control_before))]
                y = ctrl[(ctrl.treatment == base) & (ctrl.year < min(2027, control_before))]
                rep = max(rep, max_diff(x, y, cols))
        checks["replay"][cfg] = rep
        t = t.assign(config=cfg)
        traj_all.append(t)
        for a, p in meta.items():
            x = t[t.treatment == a]
            h = plog[plog.treatment == a].set_index("year").effort_hours
            r = culls[culls.treatment == a].set_index("year")
            removed = r.removed_pre_adults + r.removed_post_adults
            windows = ["E", "L"] if p["window"] == "control" else [p["window"]]
            for wkey in windows:
                (e0, e1), (s0, s1), _ = WINDOWS[wkey]
                peak = x[x.year.between(s0, s1)].groupby("reef_name").adults_ha.max()
                row = {"config": cfg, "treatment": a, "seed": p["seed"], "window": wkey,
                       "mode": p["mode"] or "none", "mult": p["mult"],
                       "effort_hours_window": float(h.loc[e0:e1].sum()) if p["window"] != "control" else 0.0,
                       "removed_window": float(removed.loc[e0:e1].sum()),
                       "max_adults_5reefs": float(peak.max()),
                       "max_adults_2027_2030": float(x[x.year.between(2027, 2030)].adults_ha.max()),
                       "max_adults_2026": float(x[x.year == 2026].adults_ha.max())}
                row |= {f"peak_{reef.split()[0].lower()}": float(peak.get(reef, np.nan)) for reef in REEFS}
                for tname, tv in THRESHOLDS.items():
                    row[f"success_{tv:g}"] = bool(peak.max() < tv)
                rows.append(row)
    m = pd.DataFrame(rows)
    out.mkdir(parents=True)
    m.to_csv(out / "buildup_arm_metrics.csv", index=False)
    traj = pd.concat(traj_all, ignore_index=True)

    # Required effort: smallest m with >= 4/5 seeds succeeding (controls give m = 0).
    req = []
    for (cfg, wkey), d in m.groupby(["config", "window"]):
        base = d[d["mode"] == "none"]
        for mode in ("F", "T"):
            dm = pd.concat([base, d[d["mode"] == mode]])
            for tname, tv in THRESHOLDS.items():
                frac = dm.groupby("mult")[f"success_{tv:g}"].mean()
                ok = frac[frac >= 0.8]
                mreq = float(ok.index.min()) if len(ok) else np.nan
                sel = dm[dm.mult == mreq] if np.isfinite(mreq) else dm[dm.mult == dm.mult.max()]
                years = WINDOWS[wkey][0][1] - WINDOWS[wkey][0][0] + 1
                req.append({"config": cfg, "config_name": CONFIG_NAME[cfg], "window": WINDOWS[wkey][2],
                            "deployment": MODE_NAME[mode], "threshold": tname, "required_multiplier": mreq,
                            "nominal_extra_hours_per_year": mreq * RECENT_HOURS if np.isfinite(mreq) else np.nan,
                            "deployed_extra_hours_window": float(sel.effort_hours_window.mean()),
                            "deployed_extra_hours_per_year": float(sel.effort_hours_window.mean() / years),
                            "removed_window_total": float(sel.removed_window.mean()),
                            "success_fraction_at_max_tested": float(frac.loc[frac.index.max()]),
                            "best_max_adults_at_max_tested": float(dm[dm.mult == dm.mult.max()].max_adults_5reefs.median())})
    req = pd.DataFrame(req)
    req.to_csv(out / "required_effort.csv", index=False)

    plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold", "text.color": INK,
                         "axes.labelcolor": INK})
    files = {}

    def save(fig, stem):
        for ext in ("png", "pdf"):
            path = out / f"{stem}.{ext}"
            fig.savefig(path, dpi=200, facecolor="#ffffff", bbox_inches="tight")
            files[path.name] = sha256(path)
        plt.close(fig)

    # Figure 1: trajectories 2018-2030, control vs early and late targeted effort at x SHOW_M.
    fig, axes = plt.subplots(3, len(REEFS), figsize=(18, 9.5), sharex=True)
    for i, cfg in enumerate(runs):
        for j, reef in enumerate(REEFS):
            ax = axes[i, j]
            d = traj[(traj.config == cfg) & (traj.reef_name == reef) & traj.year.between(2018, 2030)]
            d = d.assign(**pd.DataFrame([parts(a) for a in d.treatment], index=d.index))
            for key, color in TRAJ_COLORS.items():
                sel = d[d.window == "control"] if key == "control" else d[(d.window == key) & (d["mode"] == "T") & (d.mult == SHOW_M)]
                g = sel.groupby("year").adults_ha
                label = {"control": "observed culls only",
                         "E": f"+ targeted x{SHOW_M:g} effort 2020-25",
                         "L": f"+ targeted x{SHOW_M:g} effort 2026-30"}[key]
                ax.plot(g.mean().index, g.mean().values, color=color, linewidth=1.8, label=label)
                ax.fill_between(g.mean().index, g.min().values, g.max().values, color=color, alpha=0.12, linewidth=0)
            for tv, ls in ((THRESHOLDS["14.7/ha (0.22 COTS/tow)"], "--"), (3.0, ":")):
                ax.axhline(tv, color=MUTED, linewidth=0.9, linestyle=ls)
            ax.axvspan(2019.5, 2025.5, color="#f1d9c8", alpha=0.25, linewidth=0)
            ax.axvspan(2025.5, 2030.5, color="#cdeee0", alpha=0.25, linewidth=0)
            style(ax)
            if i == 0:
                ax.set_title(reef, loc="left")
            if j == 0:
                ax.set_ylabel(f"{CONFIG_NAME[cfg]}\nadults ha$^{{-1}}$")
    h, lab = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, lab, loc="lower center", ncol=3, frameon=False, fontsize=9)
    fig.suptitle(f"2020s build-up under extra culling effort{scenario_label} (seed mean, ribbon = seed range; "
                 "dashed = 14.7/ha outbreak, dotted = 3/ha; shaded = early and late windows)", x=0.01, ha="left", fontsize=11)
    fig.tight_layout(rect=(0, 0.04, 1, 0.96))
    save(fig, "buildup1_trajectories_2018_2030")

    # Figure 2: success fraction against extra effort.
    mults = sorted(m.mult.unique())
    fig, axes = plt.subplots(2, 3, figsize=(16, 7.5), sharex=True, sharey=True)
    for i, (tname, tv) in enumerate(THRESHOLDS.items()):
        for j, cfg in enumerate(runs):
            ax = axes[i, j]
            for (wkey, mode), color in LINE.items():
                d = m[(m.config == cfg) & (m.window == wkey) & (m["mode"].isin(["none", mode]))]
                f = d.groupby("mult")[f"success_{tv:g}"].mean().reindex(mults)
                x = np.array([0.25 if v == 0 else v for v in f.index])
                ax.plot(x, f.values, marker="o", markersize=4, color=color, linewidth=1.6,
                        label=f"{WINDOWS[wkey][2]}, {MODE_NAME[mode]}")
            ax.set_xscale("log", base=2)
            ax.set_xticks([0.25] + [v for v in mults if v > 0], ["0"] + [f"{v:g}" for v in mults if v > 0])
            ax.axhline(0.8, color=MUTED, linewidth=0.8, linestyle=":")
            style(ax)
            if i == 0:
                ax.set_title(CONFIG_NAME[cfg], loc="left")
            if j == 0:
                ax.set_ylabel(f"Share of seeds with no breach\n(threshold {tname})")
            if i == 1:
                ax.set_xlabel(f"Extra effort, multiple of recent effort ({RECENT_HOURS:,.0f} diver-h/yr)")
    h, lab = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, lab, loc="lower center", ncol=2, frameon=False, fontsize=9)
    fig.suptitle(f"Extra effort needed to keep Lizard, MacGillivray, North Direction, Eyrie and Martin below the threshold"
                 f"{scenario_label} (early: 2020-26 judged; late: 2027-30 judged)", x=0.01, ha="left", fontsize=11)
    fig.tight_layout(rect=(0, 0.08, 1, 0.95))
    save(fig, "buildup2_success_vs_effort")

    # Figure 3: peak adults on the five reefs against deployed extra effort.
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6), sharey=True)
    for j, cfg in enumerate(runs):
        ax = axes[j]
        for (wkey, mode), color in LINE.items():
            d = m[(m.config == cfg) & (m.window == wkey) & (m["mode"] == mode)]
            g = d.groupby("mult").agg(h=("effort_hours_window", "mean"), p=("max_adults_5reefs", "median"))
            ax.plot(g.h / 1e3, g.p, marker="o", markersize=4, color=color, linewidth=1.6,
                    label=f"{WINDOWS[wkey][2]}, {MODE_NAME[mode]}")
        for tv, ls in ((THRESHOLDS["14.7/ha (0.22 COTS/tow)"], "--"), (3.0, ":")):
            ax.axhline(tv, color=MUTED, linewidth=0.9, linestyle=ls)
        ax.set_xscale("log")
        ax.set_title(CONFIG_NAME[cfg], loc="left")
        ax.set_xlabel("Extra diver-hours deployed in the window (thousands, log)")
        style(ax)
    axes[0].set_ylabel("Highest adults/ha on the five reefs\nin the judged years (seed median)")
    if scenario_label:
        fig.suptitle(f"Peak adults against extra effort{scenario_label}", x=0.01, ha="left", fontsize=11)
    h, lab = axes[0].get_legend_handles_labels()
    fig.legend(h, lab, loc="lower center", ncol=2, frameon=False, fontsize=9)
    fig.tight_layout(rect=(0, 0.14, 1, 1))
    save(fig, "buildup3_peak_vs_effort")

    report = {"status": "scenario_study_not_promoted",
              "protocol": "crash_mechanism_protocol.md, 2026-10-22, thirty-second round",
              "runs": {k: v.name for k, v in runs.items()}, "scenario_label": scenario_label,
              "control_replay_before": control_before, "checks": {
                  "failures": checks["failures"], "replay_max": checks["replay"],
                  "effort_outside_window_max": max(checks["effort_outside_window"].values())},
              "thresholds_ha": THRESHOLDS, "reefs": REEFS, "recent_domain_hours_per_year": RECENT_HOURS,
              "required_effort": json.loads(req.to_json(orient="records")),
              "notes": ["Extra effort is additional to observed culls; success needs >= 4 of 5 seeds.",
                        "2027-30 uses CNRM-ESM2-1 SSP2-4.5, resampled supply scaled by domain coral; juveniles are not culled; "
                        "outside reefs are not culled.",
                        "Cull effort and 2020s densities derive from licensed data: check permissions before publication."],
              "inputs": {f"{r.name}/{n}": sha256(r / n) for r in runs.values()
                         for n in ("trajectories.csv", "trajectories_extra.csv", "culls.csv", "projection_log.csv")},
              "script_sha256": sha256(Path(__file__)),
              "outputs": files | {n: sha256(out / n) for n in ("buildup_arm_metrics.csv", "required_effort.csv")}}
    (out / "metadata.json").write_text(json.dumps(report, indent=2, default=str) + "\n")
    pd.set_option("display.width", 250)
    print(json.dumps(report["checks"], indent=1))
    print(req[["config", "window", "deployment", "threshold", "required_multiplier", "nominal_extra_hours_per_year",
               "deployed_extra_hours_per_year", "removed_window_total", "success_fraction_at_max_tested",
               "best_max_adults_at_max_tested"]].round(2).to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out", type=Path)
    for c in ("C", "G", "O"):
        parser.add_argument(f"--{c}", type=Path, required=True)
    parser.add_argument("--scenario-label", default="")
    parser.add_argument("--control-replay-before", type=int, default=2027)
    a = parser.parse_args()
    main(a.out.resolve(), {"C": a.C.resolve(), "G": a.G.resolve(), "O": a.O.resolve()},
         f" - {a.scenario_label}" if a.scenario_label else "", a.control_replay_before)
