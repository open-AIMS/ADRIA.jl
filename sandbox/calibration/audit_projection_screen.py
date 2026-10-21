"""Audit stage A, the projection scenario screen (crash_mechanism_protocol.md, "Preregistered:
projection mechanisms and scenarios", 2026-10-21).

Replays: (1) the control arms (no extra tokens) equal pilot run 20261019T_58 over 1991-2100 for the SSPs
it holds; (2) every _cullreef arm equals the lizfast culling control (OWEN_1991_LIZFAST_CULL) over
1991-2026; (3) all _cullreef arms are identical to each other before 2027.

Readouts per scenario (seed means) against base (_cullreef_fcF100_xc), per reef and period:
coral mean, years with coral < 0.05, years in outbreak (adults >= 14.7/ha), 2027-31 peak adults,
forward-cull removals (COTS), reset factors, coral-coupling factor.

    python audit_projection_screen.py <screen run> <cull control run>
Outputs go to <screen run>: screen_metrics.csv, screen_effects.csv, screen_audit.json, screen_<ssp>_<set>.png
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
from audit_1991_switching_repair import sha256  # noqa: E402
from audit_projection_pilot import CAL_REEFS, OUTBREAK_HA, PERIODS, episodes  # noqa: E402

HERE = Path(__file__).resolve().parent
PILOT = HERE / "runs" / "20261019T_58_projection_pilot_1991_v1"
TABLES = ("trajectories", "flows", "mechanisms", "trajectories_extra", "flows_extra")
BASE = "_cullreef_fcF100_xc"
LABELS = {"": "control (no culling, uncoupled)", BASE: "base: cull recent effort, coupled supply",
          "_cullreef_fcF100": "uncoupled supply", BASE + "_rsAJ": "2026 reset (adults + juveniles)",
          BASE + "_rsAO": "2026 reset (adults only)", "_cullreef_fcN_xc": "no forward culling",
          "_cullreef_fcT100_xc": "targeted culling x1", "_cullreef_fcF300_xc": "recent effort x3",
          BASE + "_ta05": "thermal acclimation 0.05/yr", BASE + "_ta12": "thermal acclimation 0.12/yr",
          BASE + "_cl10": "outside coral larvae 10%", BASE + "_cl50": "outside coral larvae 50%"}
SETS = {"cots": ["_cullreef_fcF100", BASE, BASE + "_rsAJ", BASE + "_rsAO", "_cullreef_fcN_xc",
                 "_cullreef_fcT100_xc", "_cullreef_fcF300_xc"],
        "coral": [BASE, BASE + "_ta05", BASE + "_ta12", BASE + "_cl10", BASE + "_cl50"]}
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7"]


def extra(arm: str) -> str:
    m = re.search(r"_pj\d{3}g\d{2}[WMCR](.*?)_aprmak_", arm)
    return m.group(1)


def meta(arm: str) -> dict:
    m = re.search(r"_pj(\d{3})g(\d{2})([WMCR])", arm)
    return {"ssp": "ssp" + m.group(1), "seed": int(re.search(r"_s(\d+)_", arm).group(1)), "extra": extra(arm)}


def max_diff(a: pd.DataFrame, b: pd.DataFrame, cols: list[str]) -> float:
    a, b = a.sort_values(["reef_name", "year"]), b.sort_values(["reef_name", "year"])
    if len(a) != len(b) or len(a) == 0:
        raise ValueError("Replay rows differ")
    return float(np.nanmax(np.abs(a[cols].to_numpy(float) - b[cols].to_numpy(float))))


def main(run: Path, control: Path) -> None:
    out = {n: run / n for n in ("screen_metrics.csv", "screen_effects.csv", "screen_audit.json")}
    if any(p.exists() for p in out.values()):
        raise FileExistsError("Refusing to overwrite audit")
    failures = pd.read_csv(run / "failures.csv")
    data = {n: pd.read_csv(run / f"{n}.csv") for n in TABLES}
    arms = sorted(data["trajectories"].treatment.unique())
    replay = {"control_vs_pilot": {}, "cull_vs_control": {}, "cull_arms_pre2027": {}}
    pilot = {n: pd.read_csv(PILOT / f"{n}.csv") for n in TABLES}
    ctrl = {n: pd.read_csv(control / f"{n}.csv") for n in TABLES}
    for name in TABLES:
        cols = [c for c in data[name] if c not in ("treatment", "reef_name", "year")]
        for arm in arms:
            x = data[name][data[name].treatment == arm]
            e = extra(arm)
            if e == "" and arm in set(pilot[name].treatment):
                replay["control_vs_pilot"][f"{arm}/{name}"] = max_diff(x, pilot[name][pilot[name].treatment == arm], cols)
            if e.startswith("_cullreef"):
                base_name = re.sub(r"_pj\d{3}g\d{2}[WMCR].*?_aprmak_", "_cullreef_aprmak_", arm)
                c = ctrl[name][ctrl[name].treatment == base_name]
                replay["cull_vs_control"][f"{arm}/{name}"] = max_diff(x[x.year <= 2026], c, cols)
                ref = arm.replace(e, BASE)
                y = data[name][data[name].treatment == ref]
                replay["cull_arms_pre2027"][f"{arm}/{name}"] = max_diff(x[x.year <= 2026], y[y.year <= 2026], cols)
    traj = pd.concat([data["trajectories"], data["trajectories_extra"]], ignore_index=True)
    culls = pd.read_csv(run / "culls.csv")
    plog = pd.read_csv(run / "projection_log.csv")
    resets = pd.read_csv(run / "projection_reset.csv")
    rows = []
    for arm in arms:
        m = meta(arm)
        for reef in CAL_REEFS:
            x = traj[(traj.treatment == arm) & (traj.reef_name == reef)].set_index("year").sort_index()
            c = culls[(culls.treatment == arm) & (culls.scope == reef)].set_index("year")
            for period, (y0, y1) in PERIODS.items():
                w = x.loc[y0:y1]
                eps = episodes(w.index.to_numpy(), w.adults_ha.to_numpy())
                cw = c.loc[y0:y1]
                rows.append(m | {"treatment": arm, "reef_name": reef, "period": period,
                                 "coral_mean": float(w.coral_cover.mean()),
                                 "years_coral_below_005": int((w.coral_cover < 0.05).sum()),
                                 "years_in_outbreak": int((w.adults_ha >= OUTBREAK_HA).sum()),
                                 "outbreak_episodes": len(eps),
                                 "peak_2027_31_ha": float(x.adults_ha.loc[2027:2031].max()),
                                 "adults_2026_ha": float(x.adults_ha.loc[2026]),
                                 "adults_2027_ha": float(x.adults_ha.loc[2027]),
                                 "removed_adults": float(cw.removed_pre_adults.sum() + cw.removed_post_adults.sum())
                                 if len(cw) else 0.0})
    metrics = pd.DataFrame(rows)
    metrics.to_csv(out["screen_metrics.csv"], index=False)
    num = [c for c in metrics if c not in ("ssp", "seed", "extra", "treatment", "reef_name", "period")]
    mean = metrics.groupby(["ssp", "extra", "reef_name", "period"])[num].mean().reset_index()
    base = mean[mean.extra == BASE].set_index(["ssp", "reef_name", "period"])[num]
    eff = mean.set_index(["ssp", "reef_name", "period"])
    for col in ("coral_mean", "years_coral_below_005", "years_in_outbreak", "peak_2027_31_ha"):
        eff[f"d_{col}"] = eff[col] - base[col].reindex(eff.index).values
    eff = eff.reset_index()
    eff["scenario"] = eff.extra.map(LABELS)
    eff.to_csv(out["screen_effects.csv"], index=False)

    coupling = plog.merge(pd.DataFrame({"treatment": arms, "extra": [extra(a) for a in arms],
                                        "ssp": [meta(a)["ssp"] for a in arms]}), on="treatment")
    coupling = coupling[coupling.extra.str.contains("_xc")].dropna(subset=["external_coral_factor"])
    g = coupling.groupby(["ssp", pd.cut(coupling.year, [2026, 2050, 2075, 2100])], observed=True).external_coral_factor.mean()
    reset = resets.groupby(["group", "juveniles"]).factor.agg(["mean", "min", "max"]).round(3)
    dom = culls[(culls.scope == "domain") & (culls.year >= 2027)].merge(
        pd.DataFrame({"treatment": arms, "extra": [extra(a) for a in arms], "ssp": [meta(a)["ssp"] for a in arms]}),
        on="treatment")
    dom["removed"] = dom.removed_pre_adults + dom.removed_post_adults
    removals = dom.groupby(["ssp", "extra", pd.cut(dom.year, [2026, 2031, 2050, 2100])], observed=True).removed.sum() / 5

    years = np.arange(1991, 2101)
    figs = {}
    for ssp in sorted(metrics.ssp.unique()):
        for set_name, extras in SETS.items():
            fig, axes = plt.subplots(2, 4, figsize=(17, 7), sharex=True)
            for j, reef in enumerate(CAL_REEFS):
                for i, var in enumerate(("adults_ha", "coral_cover")):
                    ax = axes[i, j]
                    for e, color in zip(extras, PALETTE):
                        sel = [a for a in arms if extra(a) == e and meta(a)["ssp"] == ssp]
                        piv = traj[traj.treatment.isin(sel) & (traj.reef_name == reef)].pivot(
                            index="year", columns="treatment", values=var).reindex(years)
                        ax.plot(piv.index, piv.mean(axis=1), color=color, linewidth=1.5, label=LABELS[e])
                    ax.axvline(2026.5, color="#52514e", linewidth=.8, linestyle=":")
                    ax.grid(color="#e4e3df", linewidth=.8)
                    ax.spines[["top", "right"]].set_visible(False)
                    ax.set_xlim(2015, 2100)
                axes[0, j].set_title(reef, loc="left", fontsize=10)
                axes[0, j].axhline(OUTBREAK_HA, color="#b6b5b0", linewidth=.8)
            axes[0, 0].set_ylabel("Adult COTS ha$^{-1}$ (seed mean)")
            axes[1, 0].set_ylabel("Coral cover (fraction of habitable area)")
            h, lab = axes[0, 0].get_legend_handles_labels()
            fig.legend(h, lab, loc="lower center", ncol=4, frameon=False, fontsize=8)
            fig.suptitle(f"Stage A screen, CNRM-ESM2-1 {ssp.upper()}: {set_name} scenarios (seed means; dotted = projection start)",
                         x=0.01, ha="left", fontsize=10)
            fig.tight_layout(rect=(0, 0.08, 1, 0.97))
            path = run / f"screen_{ssp}_{set_name}.png"
            fig.savefig(path, dpi=150, facecolor="#fcfcfb")
            plt.close(fig)
            figs[path.name] = sha256(path)

    report = {"status": "projection_scenario_screen_not_promoted",
              "protocol": "crash_mechanism_protocol.md, Preregistered: projection mechanisms and scenarios (2026-10-21)",
              "failures": int(len(failures)),
              "replay_max": {k: (max(v.values()) if v else None) for k, v in replay.items()},
              "replay_counts": {k: len(v) for k, v in replay.items()},
              "coupling_factor_mean": {f"{k[0]} {k[1]}": float(v) for k, v in g.items()},
              "reset_factors": json.loads(reset.reset_index().to_json(orient="records")),
              "forward_removals_domain_per_seed": {f"{k[0]} {k[1]} {k[2]}": float(v) for k, v in removals.items()},
              "inputs": {n: sha256(run / f"{n}.csv") for n in TABLES + ("culls", "projection_log", "projection_reset")} |
                        {"control/trajectories.csv": sha256(control / "trajectories.csv")},
              "script_sha256": sha256(Path(__file__)),
              "outputs": {n: sha256(p) for n, p in out.items() if n != "screen_audit.json"} | figs}
    out["screen_audit.json"].write_text(json.dumps(report, indent=2, default=str) + "\n")
    pd.set_option("display.width", 260)
    show = eff[eff.extra != BASE][["ssp", "scenario", "reef_name", "period", "d_coral_mean", "d_years_coral_below_005",
                                   "d_years_in_outbreak", "d_peak_2027_31_ha"]]
    print(show[show.period != "2076-2100"].round(2).to_string(index=False))
    print(json.dumps({k: report[k] for k in ("failures", "replay_max", "replay_counts", "coupling_factor_mean")}, indent=1))
    print(reset.to_string())
    print(removals.round(0).to_string())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("control", type=Path)
    a = parser.parse_args()
    main(a.run.resolve(), a.control.resolve())
