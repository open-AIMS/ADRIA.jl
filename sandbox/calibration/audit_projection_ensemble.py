"""Audit stage B, the projection ensemble (crash_mechanism_protocol.md, "Preregistered: projection
mechanisms and scenarios", 2026-10-21): 5 GCMs x SSP2-4.5 / SSP3-7.0 x 5 seeds x reset {off, on} x
forward cull {none; fixed x 1, 1.5, 2, 3; targeted x 1, 1.5, 2, 3}, coupled outside supply.

Checks: no failures; every arm equals the lizfast culling control (run 20261021T_01) over 1991-2026.
Metrics per arm x calibration reef x period (audit_projection_pilot definitions); summaries over
GCM x seed (median, 10th and 90th percentiles) per SSP x reset x cull level, averaged over the four
calibration reefs; culling benefit = change against no forward culling in the same GCM, seed and reset,
with coral-years gained per 1,000 COTS removed (2027-2050, domain removals).

    python audit_projection_ensemble.py <out_dir> <run> [<run> ...]
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
CONTROL = HERE / "runs" / "20261021T_01_lizfast_cull_control_1991_v1"
PROJ = HERE.parents[1] / "sandbox" / "data" / "derived" / "projection_inputs_v0.5"
TABLES = ("trajectories", "flows", "mechanisms")
CULL_ORDER = ["fcN", "fcF100", "fcF150", "fcF200", "fcF300", "fcT100", "fcT150", "fcT200", "fcT300"]
CULL_LABEL = {"fcN": "none", "fcF100": "recent x1", "fcF150": "recent x1.5", "fcF200": "recent x2", "fcF300": "recent x3",
              "fcT100": "targeted x1", "fcT150": "targeted x1.5", "fcT200": "targeted x2", "fcT300": "targeted x3"}
RESET_COLORS = {"none": "#2a78d6", "2026 reset": "#eb6834"}


def parts(arm: str) -> dict:
    m = re.search(r"_pj(\d{3})g(\d{2})W_cullreef_(fc[NFT]\d{0,3})_xc(_rsAJ)?_aprmak_", arm)
    if m is None:
        raise ValueError(f"Not an ensembleB arm: {arm}")
    return {"ssp": "ssp" + m.group(1), "gcm": int(m.group(2)), "cull": m.group(3),
            "reset": "2026 reset" if m.group(4) else "none", "seed": int(re.search(r"_s(\d+)_", arm).group(1))}


def q(x: pd.Series, p: float) -> float:
    return float(np.nanquantile(x, p))


def main(out: Path, runs: list[Path]) -> None:
    if out.exists():
        raise FileExistsError(f"Refusing to overwrite {out}")
    names = pd.read_csv(PROJ / "dhw" / "dhw_scenarios.csv")
    ctrl = {n: pd.read_csv(CONTROL / f"{n}.csv") for n in TABLES}
    rows, replay, failures, removals = [], 0.0, 0, []
    for run in runs:
        failures += len(pd.read_csv(run / "failures.csv"))
        data = {n: pd.read_csv(run / f"{n}.csv") for n in TABLES}
        culls = pd.read_csv(run / "culls.csv")
        traj = data["trajectories"]
        for arm in sorted(traj.treatment.unique()):
            p = parts(arm)
            base = re.sub(r"_pj\d{3}g\d{2}W_cullreef_.*?_aprmak_", "_cullreef_aprmak_", arm)
            for n in TABLES:
                cols = [c for c in data[n] if c not in ("treatment", "reef_name", "year")]
                a = data[n][(data[n].treatment == arm) & (data[n].year <= 2026)].sort_values(["reef_name", "year"])
                b = ctrl[n][ctrl[n].treatment == base].sort_values(["reef_name", "year"])
                if len(a) != len(b) or len(b) == 0:
                    raise ValueError(f"Replay rows differ for {arm}/{n}")
                replay = max(replay, float(np.nanmax(np.abs(a[cols].to_numpy(float) - b[cols].to_numpy(float)))))
            c = culls[(culls.treatment == arm) & (culls.scope == "domain")].set_index("year")
            rem = c.removed_pre_adults + c.removed_post_adults
            for period, (y0, y1) in PERIODS.items():
                removals.append(p | {"period": period, "removed_domain": float(rem.loc[y0:y1].sum())})
            for reef in CAL_REEFS:
                x = traj[(traj.treatment == arm) & (traj.reef_name == reef)].set_index("year").sort_index()
                for period, (y0, y1) in PERIODS.items():
                    w = x.loc[y0:y1]
                    rows.append(p | {"reef_name": reef, "period": period,
                                     "coral_mean": float(w.coral_cover.mean()),
                                     "years_coral_below_005": int((w.coral_cover < 0.05).sum()),
                                     "years_in_outbreak": int((w.adults_ha >= OUTBREAK_HA).sum()),
                                     "outbreak_episodes": len(episodes(w.index.to_numpy(), w.adults_ha.to_numpy())),
                                     "peak_2027_31_ha": float(x.adults_ha.loc[2027:2031].max())})
    m = pd.DataFrame(rows)
    rem = pd.DataFrame(removals)
    out.mkdir(parents=True)
    m.to_csv(out / "ensemble_metrics.csv", index=False)
    key = ["ssp", "gcm", "seed", "reset", "cull", "period"]
    reefmean = m.groupby(key)[["coral_mean", "years_coral_below_005", "years_in_outbreak", "peak_2027_31_ha"]].mean().reset_index()
    reefmean = reefmean.merge(rem, on=key)
    none = reefmean[reefmean.cull == "fcN"].set_index(["ssp", "gcm", "seed", "reset", "period"])
    idx = reefmean.set_index(["ssp", "gcm", "seed", "reset", "period"]).index
    for col in ("coral_mean", "years_coral_below_005", "years_in_outbreak", "peak_2027_31_ha"):
        reefmean[f"d_{col}"] = reefmean[col].values - none[col].reindex(idx).values
    span = {p: y1 - y0 + 1 for p, (y0, y1) in PERIODS.items()}
    reefmean["coral_years_per_1000_removed"] = np.where(
        reefmean.removed_domain > 0, 1000 * reefmean.d_coral_mean * reefmean.period.map(span) / reefmean.removed_domain, np.nan)
    reefmean.to_csv(out / "ensemble_reefmean.csv", index=False)
    stats = []
    for g, d in reefmean.groupby(["ssp", "reset", "cull", "period"]):
        row = dict(zip(["ssp", "reset", "cull", "period"], g)) | {"n": len(d)}
        for col in ("coral_mean", "years_coral_below_005", "years_in_outbreak", "peak_2027_31_ha", "removed_domain",
                    "d_coral_mean", "d_years_in_outbreak", "d_peak_2027_31_ha", "coral_years_per_1000_removed"):
            row |= {f"{col}_med": q(d[col], .5), f"{col}_p10": q(d[col], .1), f"{col}_p90": q(d[col], .9)}
        stats.append(row)
    summary = pd.DataFrame(stats)
    summary["cull_label"] = summary.cull.map(CULL_LABEL)
    summary.to_csv(out / "ensemble_summary.csv", index=False)

    panels = [("years_in_outbreak", "2027-2050", "Years in outbreak 2027-50 (reef mean)"),
              ("peak_2027_31_ha", "2027-2050", "Peak adults 2027-31 (adults/ha)"),
              ("coral_mean", "2027-2050", "Mean coral 2027-50"),
              ("years_coral_below_005", "2051-2075", "Years coral < 0.05, 2051-75")]
    figs = {}
    for ssp in sorted(summary.ssp.unique()):
        fig, axes = plt.subplots(1, 4, figsize=(18, 4.6))
        for ax, (col, period, title) in zip(axes, panels):
            for k, (reset, color) in enumerate(RESET_COLORS.items()):
                d = summary[(summary.ssp == ssp) & (summary.reset == reset) & (summary.period == period)]
                d = d.set_index("cull").reindex(CULL_ORDER)
                x = np.arange(len(CULL_ORDER)) + (k - 0.5) * 0.22
                ax.errorbar(x, d[f"{col}_med"], yerr=[d[f"{col}_med"] - d[f"{col}_p10"], d[f"{col}_p90"] - d[f"{col}_med"]],
                            fmt="o", color=color, markersize=5, linewidth=1.5, capsize=0, label=f"start: {reset}")
            ax.set_xticks(range(len(CULL_ORDER)), [CULL_LABEL[c] for c in CULL_ORDER], rotation=45, ha="right", fontsize=8)
            ax.set_title(title, loc="left", fontsize=10)
            ax.grid(axis="y", color="#e4e3df", linewidth=.8)
            ax.spines[["top", "right"]].set_visible(False)
        axes[0].legend(frameon=False, fontsize=8)
        fig.suptitle(f"Ensemble {ssp.upper()}: 5 GCMs x 5 seeds, culling level (forward from 2027); median and 10-90% range",
                     x=0.01, ha="left", fontsize=10)
        fig.tight_layout(rect=(0, 0, 1, 0.95))
        path = out / f"ensemble_{ssp}.png"
        fig.savefig(path, dpi=150, facecolor="#fcfcfb")
        plt.close(fig)
        figs[path.name] = sha256(path)

    gcm_names = {int(r.scenario): r.gcm for r in names[names.ssp == "ssp245"].itertuples()}
    report = {"status": "projection_ensemble_not_promoted",
              "protocol": "crash_mechanism_protocol.md, Preregistered: projection mechanisms and scenarios (2026-10-21), stage B",
              "runs": [r.name for r in runs], "arms": int(len(m) / len(CAL_REEFS) / len(PERIODS)),
              "gcms": {int(g): gcm_names[int(g)] for g in sorted(m.gcm.unique())},
              "failures": failures, "replay_max_vs_cull_control": replay,
              "inputs": {f"{r.name}/{n}.csv": sha256(r / f"{n}.csv") for r in runs for n in TABLES + ("culls",)},
              "script_sha256": sha256(Path(__file__)),
              "outputs": {n: sha256(out / n) for n in ("ensemble_metrics.csv", "ensemble_reefmean.csv", "ensemble_summary.csv")} | figs}
    (out / "ensemble_audit.json").write_text(json.dumps(report, indent=2, default=str) + "\n")
    pd.set_option("display.width", 250)
    show = summary[summary.period.isin(["2027-2050", "2051-2075"])]
    print(show[["ssp", "reset", "cull_label", "period", "years_in_outbreak_med", "peak_2027_31_ha_med", "coral_mean_med",
                "coral_mean_p10", "coral_mean_p90", "years_coral_below_005_med", "removed_domain_med",
                "coral_years_per_1000_removed_med"]].round(3).to_string(index=False))
    print(json.dumps({k: report[k] for k in ("arms", "gcms", "failures", "replay_max_vs_cull_control")}, indent=1))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out", type=Path)
    parser.add_argument("runs", type=Path, nargs="+")
    a = parser.parse_args()
    main(a.out.resolve(), [r.resolve() for r in a.runs])
