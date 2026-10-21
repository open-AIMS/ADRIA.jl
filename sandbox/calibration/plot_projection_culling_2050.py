"""Report figures: projected COTS and coral to 2050 under forward culling strategies (W-wave-lf04 projection
ensemble, stage B; crash_mechanism_protocol.md, "Stage B ensemble: results", 2026-10-21).

Ensemble: 5 GCMs (MPI-ESM1-2-HR, MIROC6, CNRM-ESM2-1, EC-Earth3-Veg, UKESM1-0-LL; bias-corrected CMIP6 DHW)
x 5 seeds per SSP; observed culls to 2026, then forward effort culling from 2027; outside supply resampled from
1991-2026 waves and scaled by projected domain coral; no 2026 reset (continuous run from 1991).
Strategies shown: no forward culling; recent effort (2023-25 diver-hours per reef) x1 and x3; targeted
(the same total effort x1 moved each half-year to reefs at >= 3 adults/ha).

Figures (central line = median over GCM x seed; ribbon = 25-75% for no culling and recent effort x3):
- proj1_cots_coral_2050_<ssp>: adults/ha and coral cover, 2015-2050, five reefs (four calibration reefs + Martin,
  the most-culled domain reef);
- proj2_culling_effects_2027_2050: effect of every culling level on years in outbreak, peak adults 2027-31,
  mean coral and COTS removed (domain), both SSPs, median and 10-90% over GCM x seed.
Table: culling_effects_2027_2050.csv (per reef, strategy and SSP).

    python plot_projection_culling_2050.py <out_dir> <stage B run> [<stage B run> ...]
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

INK, MUTED, GRID = "#1f1e1d", "#6f6e69", "#e4e3df"
OUTBREAK = 0.22 / 0.015
REEFS = ["Lizard Island Reef", "MacGillivray Reef", "North Direction Reef", "Eyrie Reef", "Martin Reef"]
SHOW = {"fcN": ("no forward culling", "#2a78d6"), "fcF100": ("recent effort x1", "#eb6834"),
        "fcF300": ("recent effort x3", "#1baf7a"), "fcT100": ("targeted x1", "#eda100")}
RIBBON = ("fcN", "fcF300")
ORDER = ["fcN", "fcF100", "fcF150", "fcF200", "fcF300", "fcT100", "fcT150", "fcT200", "fcT300"]
LABEL = {"fcN": "none", "fcF100": "recent x1", "fcF150": "recent x1.5", "fcF200": "recent x2", "fcF300": "recent x3",
         "fcT100": "targeted x1", "fcT150": "targeted x1.5", "fcT200": "targeted x2", "fcT300": "targeted x3"}
SSP_COLORS = {"ssp245": "#2a78d6", "ssp370": "#eb6834"}
SSP_NAME = {"ssp245": "SSP2-4.5", "ssp370": "SSP3-7.0"}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parts(arm: str):
    m = re.search(r"_pj(\d{3})g(\d{2})W_cullreef_(fc[NFT]\d{0,3})_xc(_rsAJ)?_aprmak_", arm)
    return None if m is None else {"ssp": "ssp" + m.group(1), "gcm": int(m.group(2)), "cull": m.group(3),
                                   "reset": bool(m.group(4)), "seed": int(re.search(r"_s(\d+)_", arm).group(1))}


def style(ax):
    ax.grid(color=GRID, linewidth=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=8)


def main(out: Path, runs: list[Path]) -> None:
    if out.exists():
        raise FileExistsError(f"Refusing to overwrite {out}")
    frames, culls = [], []
    for run in runs:
        if len(pd.read_csv(run / "failures.csv")):
            raise ValueError(f"{run.name} has failed arms")
        t = pd.concat([pd.read_csv(run / "trajectories.csv"), pd.read_csv(run / "trajectories_extra.csv")])
        frames.append(t[t.reef_name.isin(REEFS) & t.year.between(2015, 2050)])
        c = pd.read_csv(run / "culls.csv")
        culls.append(c[(c.scope == "domain") & c.year.between(2027, 2050)])
    traj = pd.concat(frames, ignore_index=True)
    meta = pd.DataFrame([{"treatment": a} | parts(a) for a in traj.treatment.unique() if parts(a)])
    meta = meta[~meta.reset]
    traj = traj.merge(meta, on="treatment")
    cull = pd.concat(culls).merge(meta, on="treatment")
    cull["removed"] = cull.removed_pre_adults + cull.removed_post_adults
    if meta.groupby(["ssp", "cull"]).size().nunique() != 1:
        raise ValueError("Unbalanced ensemble")
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

    for ssp in ("ssp245", "ssp370"):
        fig, axes = plt.subplots(2, len(REEFS), figsize=(18, 6.6), sharex=True)
        for j, reef in enumerate(REEFS):
            for i, var in enumerate(("adults_ha", "coral_cover")):
                ax = axes[i, j]
                for code, (label, color) in SHOW.items():
                    d = traj[(traj.ssp == ssp) & (traj.cull == code) & (traj.reef_name == reef)]
                    g = d.groupby("year")[var]
                    ax.plot(g.median().index, g.median().values, color=color, linewidth=1.8, label=label)
                    if code in RIBBON:
                        ax.fill_between(g.median().index, g.quantile(.25).values, g.quantile(.75).values,
                                        color=color, alpha=0.15, linewidth=0)
                ax.axvline(2026.5, color=MUTED, linewidth=0.8, linestyle=":")
                style(ax)
                ax.set_xlim(2015, 2050)
            axes[0, j].axhline(OUTBREAK, color="#b6b5b0", linewidth=0.9)
            axes[0, j].set_title(reef, loc="left")
        axes[0, 0].set_ylabel("Adult COTS ha$^{-1}$")
        axes[1, 0].set_ylabel("Coral cover (fraction of habitable area)")
        axes[0, 0].text(2015.5, OUTBREAK, " outbreak (0.22 COTS/tow)", color=MUTED, fontsize=7, va="bottom")
        axes[0, 1].text(2026.8, axes[0, 1].get_ylim()[1] * 0.95, "projection →", color=MUTED, fontsize=7, va="top")
        h, lab = axes[0, 0].get_legend_handles_labels()
        fig.legend(h, lab, loc="lower center", ncol=4, frameon=False, fontsize=9)
        fig.suptitle(f"Projected COTS and coral to 2050 under forward culling, {SSP_NAME[ssp]} "
                     "(median of 5 climate models x 5 seeds; ribbons 25-75%)", x=0.01, ha="left", fontsize=11)
        fig.tight_layout(rect=(0, 0.05, 1, 0.96))
        save(fig, f"proj1_cots_coral_2050_{ssp}")

    rows = []
    for (ssp, code, gcm, seed, arm), d in traj.groupby(["ssp", "cull", "gcm", "seed", "treatment"]):
        removed = float(cull[cull.treatment == arm].removed.sum())
        for reef, x in d.groupby("reef_name"):
            x = x.set_index("year").sort_index()
            p = x.loc[2027:2050]
            rows.append({"ssp": ssp, "cull": code, "gcm": gcm, "seed": seed, "reef_name": reef,
                         "years_in_outbreak": int((p.adults_ha >= OUTBREAK).sum()),
                         "peak_2027_31_ha": float(x.adults_ha.loc[2027:2031].max()),
                         "coral_mean": float(p.coral_cover.mean()), "removed_domain": removed})
    m = pd.DataFrame(rows)
    table = (m.groupby(["ssp", "cull", "reef_name"])[["years_in_outbreak", "peak_2027_31_ha", "coral_mean"]]
             .median().reset_index())
    table["strategy"] = table.cull.map(LABEL)
    table.to_csv(out / "culling_effects_2027_2050.csv", index=False)
    files["culling_effects_2027_2050.csv"] = sha256(out / "culling_effects_2027_2050.csv")

    reefmean = m.groupby(["ssp", "cull", "gcm", "seed"]).agg(
        years_in_outbreak=("years_in_outbreak", "mean"), peak_2027_31_ha=("peak_2027_31_ha", "mean"),
        coral_mean=("coral_mean", "mean"), removed_domain=("removed_domain", "first")).reset_index()
    reefmean["removed_domain"] /= 1e3
    panels = [("years_in_outbreak", "Years in outbreak, 2027-50"), ("peak_2027_31_ha", "Peak adults 2027-31 (ha$^{-1}$)"),
              ("coral_mean", "Mean coral cover, 2027-50"), ("removed_domain", "COTS removed 2027-50 (thousands, domain)")]
    fig, axes = plt.subplots(1, 4, figsize=(18, 4.4))
    for ax, (col, title) in zip(axes, panels):
        for k, ssp in enumerate(("ssp245", "ssp370")):
            d = reefmean[reefmean.ssp == ssp].groupby("cull")[col]
            med, lo, hi = (d.median().reindex(ORDER), d.quantile(.1).reindex(ORDER), d.quantile(.9).reindex(ORDER))
            x = np.arange(len(ORDER)) + (k - 0.5) * 0.24
            ax.errorbar(x, med, yerr=[med - lo, hi - med], fmt="o", color=SSP_COLORS[ssp], markersize=5, linewidth=1.5,
                        capsize=0, label=SSP_NAME[ssp])
        ax.axvline(4.5, color=GRID, linewidth=1.2)
        ax.set_xticks(range(len(ORDER)), [LABEL[c] for c in ORDER], rotation=45, ha="right", fontsize=8)
        ax.set_title(title, loc="left")
        style(ax)
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle("Effect of forward culling, 2027-2050 (mean of Lizard, MacGillivray, North Direction, Eyrie and Martin; "
                 "median and 10-90% over 5 climate models x 5 seeds)", x=0.01, ha="left", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    save(fig, "proj2_culling_effects_2027_2050")

    report = {"status": "report_figures_projection_not_promoted",
              "protocol": "crash_mechanism_protocol.md, Stage B ensemble: results (2026-10-21)",
              "model": "W-wave-lf04; observed culls to 2026, forward effort culling from 2027; coupled outside supply; no reset",
              "runs": [r.name for r in runs], "arms_used": int(len(meta)), "reefs": REEFS,
              "notes": ["Projections are scenario-conditional, not forecasts: outside supply is resampled history scaled by "
                        "domain coral; no thermal adaptation; the model's 2026 overshoot carries into 2027-31.",
                        "Recent effort = 2023-25 diver-hours per reef from the CSIRO cull sheet (licensed); check permissions "
                        "before publication."],
              "inputs": {f"{r.name}/{n}": sha256(r / n) for r in runs for n in ("trajectories.csv", "trajectories_extra.csv", "culls.csv")},
              "script_sha256": sha256(Path(__file__)), "outputs": files}
    (out / "metadata.json").write_text(json.dumps(report, indent=2) + "\n")
    pd.set_option("display.width", 220)
    print(table[table.cull.isin(SHOW)].pivot_table(index=["ssp", "reef_name"], columns="strategy",
                                                    values=["years_in_outbreak", "coral_mean"]).round(2).to_string())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out", type=Path)
    parser.add_argument("runs", type=Path, nargs="+")
    a = parser.parse_args()
    main(a.out.resolve(), [r.resolve() for r in a.runs])
