"""Technical-report figures for the W-wave calibration (2026-10-19): how well the calibration reefs and the
additional reported reefs are represented, with uncertainty ribbons.

Model: W-wave (W + h10 start, Q 25, Ext x5, correction weight 3 layer owen_global_2018_2023_corr_w3, new coral,
chl aprmak, wave-exposure settlement floor 0.25; culling off), run `sens` with OWEN_1991_SENS_WAVE=true:
- central line: mean of the 5 base seed pairs (hydrodynamic-year sequence + ADRIA RNG);
- inner ribbon: min-max over those 5 seeds (stochastic uncertainty at the calibrated parameters);
- outer ribbon: full range (min-max) over the 12 one-at-a-time +/-25% parameter perturbations x 5 seeds
  (external multiplier, chronic coral mortality, low-DHW bleaching scale, chl slope, pit capacity Q,
  initial fast-coral shares) plus the base arms. This is a sensitivity envelope, not a probability interval.

Observations:
- AIMS LTMP manta tows (reef_cots.csv; COTS per tow by reporting year); simulated CPUE = 0.015 x adults/ha
  (the frozen observation model, fixed_q0015); 0.22 COTS/tow outbreak line;
- 2020s diver/cull density targets (density_targets_2020s; SALAD diver densities and cull catch converted to
  COTS/ha). Derived from unpublished SALAD / CSIRO cull data: check permissions before publication;
- AIMS LTMP manta hard coral (reef_manta.csv, modelled median and 95% interval) and 9 m photo-transect hard
  coral (reef_photo_transect.csv). Model coral is a fraction of habitable area, so levels differ by support.

With --cull-run (a W-wave sens run with OWEN_1991_SENS_CULL=true: observed culls, reef-proportional allocation),
the culling-on ensemble is overlaid in purple (dashed seed-mean line, seed range and parameter envelope).

    python plot_calibration_report.py --run <sens W-wave run> [--cull-run <sens W-wave culling run>] [--output <dir>]
Outputs (git-ignored): sandbox/data/derived/calibration_report_figures/ (PNG 200 dpi + PDF, summary CSVs).
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "sandbox" / "data"
OBS_COTS = DATA / "reef_cots.csv"
OBS_MANTA = DATA / "reef_manta.csv"
OBS_PHOTO = DATA / "reef_photo_transect.csv"
TARGETS = DATA / "derived" / "density_targets_2020s" / "targets_reef_year.csv"
DEFAULT_OUT = DATA / "derived" / "calibration_report_figures"
Q, OUTBREAK = 0.015, 0.22
CALIB = {"Lizard Island Reef": "Lizard Isles", "MacGillivray Reef": "Macgillivray Reef",
         "North Direction Reef": "North Direction Island", "Eyrie Reef": "Eyrie Reef"}
EXTRA = {"Martin Reef": "Martin Reef", "Linnet Reef": "Linnet Reef", "Yonge Reef": "Yonge Reef",
         "Carter Reef": "Carter Reef", "No Name Reef": "No Name Reef", "Ribbon No 10 Reef": None}
EXTRA_NOTE = {"Martin Reef": "outbreak reef", "Linnet Reef": "outbreak reef", "Yonge Reef": "never-outbreak (ribbon)",
              "Carter Reef": "never-outbreak (ribbon)", "No Name Reef": "never-outbreak (ribbon)",
              "Ribbon No 10 Reef": "ribbon; no LTMP surveys"}
INK, MUTED, GRID, SURFACE = "#1f1e1d", "#6f6e69", "#e4e3df", "#ffffff"
LINE, INNER, OUTER, OBS, TARGET = "#0d366b", "#86b6ef", "#d6e6f9", "#1f1e1d", "#eb6834"
CLINE, CINNER, COUTER = "#5b2a86", "#a77fd6", "#d9c8f0"
YEARS = (1992, 2026)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def style(ax):
    ax.grid(color=GRID, linewidth=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=8)
    ax.set_xlim(YEARS[0] - 0.5, YEARS[1] + 0.5)


def bands(df: pd.DataFrame, value: str, base: list[str], all_arms: list[str], scale: float = 1.0) -> pd.DataFrame:
    p = df.pivot(index="year", columns="treatment", values=value) * scale
    b = p[base]
    a = p[all_arms]
    return pd.DataFrame({"mean": b.mean(axis=1), "lo": b.min(axis=1), "hi": b.max(axis=1),
                         "olo": a.min(axis=1), "ohi": a.max(axis=1)})


PREFIX = {"": ""}


def draw(ax, b: pd.DataFrame, label_prefix: str = ""):
    b = b.loc[YEARS[0]:YEARS[1]]
    pre = PREFIX[""]
    ax.fill_between(b.index, b.olo, b.ohi, color=OUTER, linewidth=0, label=f"{pre}±25% parameter envelope (full range)")
    ax.fill_between(b.index, b.lo, b.hi, color=INNER, linewidth=0, label=f"{pre}seed range (5 seeds)")
    ax.plot(b.index, b["mean"], color=LINE, linewidth=2.0, label=f"{pre}model (seed mean)")


def draw_cull(ax, b: pd.DataFrame, lo_clip: float | None = None, years=YEARS):
    b = b.loc[years[0]:years[1]]
    c = (lambda x: x.clip(lower=lo_clip)) if lo_clip is not None else (lambda x: x)
    ax.fill_between(b.index, c(b.olo), c(b.ohi), color=COUTER, alpha=0.55, linewidth=0, label="culling on: parameter envelope")
    ax.fill_between(b.index, c(b.lo), c(b.hi), color=CINNER, alpha=0.55, linewidth=0, label="culling on: seed range")
    ax.plot(b.index, c(b["mean"]), color=CLINE, linewidth=2.0, linestyle=(0, (5, 2)), label="culling on: model (seed mean)")


def fit_stats(model: pd.Series, obs: pd.Series, threshold: float | None = None, log_offset: float | None = None) -> dict:
    """Agreement between the model seed mean and observations in matching years (2026-10-22): n, bias
    (model - obs), MAE, Pearson r, Spearman rho (the four annotated); also relative bias (%), Nash-Sutcliffe
    efficiency, MAE of log(x + offset) and, with a threshold, outbreak-state skill (true skill statistic)."""
    pair = pd.concat([model.rename("m"), obs.rename("o")], axis=1, join="inner").dropna()
    n = len(pair)
    out = {"n": n}
    if n < 3:
        return out | {k: np.nan for k in ("bias", "mae", "pearson", "spearman", "rel_bias_pct", "nse", "log_mae", "tss")}
    m, o = pair.m.to_numpy(float), pair.o.to_numpy(float)
    const = np.ptp(m) == 0 or np.ptp(o) == 0
    out |= {"bias": float(np.mean(m - o)), "mae": float(np.mean(np.abs(m - o))),
            "pearson": np.nan if const else float(pearsonr(m, o)[0]),
            "spearman": np.nan if const else float(spearmanr(m, o)[0]),
            "rel_bias_pct": float(100 * np.mean(m - o) / np.mean(o)) if np.mean(o) > 0 else np.nan,
            "nse": float(1 - np.sum((m - o) ** 2) / np.sum((o - o.mean()) ** 2)) if np.ptp(o) > 0 else np.nan,
            "log_mae": float(np.mean(np.abs(np.log(m + log_offset) - np.log(o + log_offset)))) if log_offset else np.nan}
    if threshold is not None:
        mo, oo = m >= threshold, o >= threshold
        hits, misses = np.sum(mo & oo), np.sum(~mo & oo)
        fa, cn = np.sum(mo & ~oo), np.sum(~mo & ~oo)
        out["tss"] = (float(hits / (hits + misses) - fa / (fa + cn))
                      if (hits + misses) > 0 and (fa + cn) > 0 else np.nan)
    else:
        out["tss"] = np.nan
    return out


def stat_text(st: dict, compact: bool = False) -> str:
    if st.get("n", 0) < 3:
        return f"n={st.get('n', 0)}"
    signed = lambda v: "-" if not np.isfinite(v) else f"{v:+.2f}"
    plain = lambda v: "-" if not np.isfinite(v) else f"{v:.2f}"
    if compact:
        return (f"n={st['n']} bias {signed(st['bias'])} MAE {plain(st['mae'])}\n"
                f"     r {plain(st['pearson'])}  \u03c1 {plain(st['spearman'])}")
    return (f"n={st['n']:<2d} bias {signed(st['bias'])}  MAE {plain(st['mae'])}  "
            f"r {plain(st['pearson'])}  \u03c1 {plain(st['spearman'])}")


def annotate(ax, lines: list[str], loc: str = "upper left"):
    """Small statistics box drawn in ink colour (text never takes a series colour)."""
    y0, va = (0.97, "top") if loc.startswith("upper") else (0.03, "bottom")
    x0, ha = (0.02, "left") if loc.endswith("left") else (0.98, "right")
    ax.text(x0, y0, "\n".join(lines), transform=ax.transAxes, ha=ha, va=va, fontsize=6.6, color=INK,
            family="monospace", zorder=10,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="#ffffff", edgecolor=GRID, alpha=0.9))


def peak_table(model: pd.Series, obs: pd.Series, windows) -> list[dict]:
    rows = []
    for name, (y0, y1) in windows.items():
        m, o = model.loc[y0:y1], obs.loc[y0:y1] if len(obs) else obs
        rows.append({"wave": name, "model_peak_year": int(m.idxmax()), "model_peak": float(m.max()),
                     "obs_peak_year": int(o.idxmax()) if len(o) and o.max() > 0 else None,
                     "obs_peak": float(o.max()) if len(o) else np.nan})
    return rows


def main(run: Path, out: Path, cull_run: Path | None = None, model_label: str = "W-wave") -> None:
    if out.exists():
        raise FileExistsError(f"Refusing to overwrite {out}")
    traj = pd.read_csv(run / "trajectories.csv")
    extra = pd.read_csv(run / "trajectories_extra.csv")
    if len(pd.read_csv(run / "failures.csv")):
        raise ValueError("Run has failed arms")
    arms = sorted(set(traj.treatment))
    base = [a for a in arms if "_sx" not in a and "_sh" not in a]
    if len(base) != 5 or not all(a.endswith("_wv025") for a in arms) or any("_cullreef_" in a for a in arms):
        raise ValueError("Expected the W-wave sens run (5 base arms, all _wv025, no culling)")
    if cull_run is not None:
        ctraj = pd.read_csv(cull_run / "trajectories.csv")
        cextra = pd.read_csv(cull_run / "trajectories_extra.csv")
        if len(pd.read_csv(cull_run / "failures.csv")):
            raise ValueError("Culling run has failed arms")
        carms = sorted(set(ctraj.treatment))
        cbase = [a for a in carms if "_sx" not in a and "_sh" not in a]
        if len(cbase) != 5 or not all("_cullreef_" in a and a.endswith("_wv025") for a in carms) or len(carms) != len(arms):
            raise ValueError("Expected the W-wave culling sens run (all arms _cullreef_ and _wv025)")
    cull_bands = (lambda d, value, scale=1.0: bands(d, value, cbase, carms, scale)) if cull_run is not None else None
    PREFIX[""] = "culling off: " if cull_run is not None else ""
    obs = pd.read_csv(OBS_COTS)
    manta = pd.read_csv(OBS_MANTA)
    photo = pd.read_csv(OBS_PHOTO, low_memory=False)
    photo = photo[(photo.variable == "HARD CORAL") & (photo.purpose == "GROUP_LEVEL") & (photo.project_code == "LTMP")]
    targets = pd.read_csv(TARGETS)
    out.mkdir(parents=True)
    plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold", "text.color": INK,
                         "axes.labelcolor": INK})
    files = []
    windows = {"1990s": (1993, 2003), "2010s": (2008, 2016), "2020s": (2019, 2026)}
    peaks = []
    stats = []

    def add_stats(figure, reef, observation, model, obsv, cmodel=None, threshold=None, log_offset=None, label="",
                  compact=False):
        """Fit statistics for culling off (and on, when a culling run is given); returns the annotation lines."""
        st = fit_stats(model, obsv, threshold, log_offset)
        stats.append({"figure": figure, "reef": reef, "observation": observation, "culling": "off", **st})
        prefix = (label + " ") if label else ""
        if cmodel is None:
            return [prefix + stat_text(st, compact)]
        cst = fit_stats(cmodel, obsv, threshold, log_offset)
        stats.append({"figure": figure, "reef": reef, "observation": observation, "culling": "on", **cst})
        if compact:
            return [f"off {stat_text(st, True)}", f"on  {stat_text(cst, True)}"]
        return [f"{prefix}cull off: {stat_text(st)}", f"{prefix}cull on:  {stat_text(cst)}"]

    def save(fig, stem):
        for ext in ("png", "pdf"):
            p = out / f"{stem}.{ext}"
            fig.savefig(p, dpi=200, facecolor=SURFACE, bbox_inches="tight")
            files.append(p)
        plt.close(fig)

    # Figure 1: calibration reefs, manta CPUE.
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True)
    for ax, (reef, oname) in zip(axes.flat, CALIB.items()):
        d = traj[traj.reef_name == reef]
        b = bands(d, "adults_ha", base, arms, Q)
        draw(ax, b)
        if cull_bands:
            draw_cull(ax, cull_bands(ctraj[ctraj.reef_name == reef], "adults_ha", Q))
        o = obs[obs.reef_name == oname].groupby("year").cotsptow.mean()
        ax.scatter(o.index, o.values, s=18, color=OBS, zorder=5, label="LTMP manta tow (mean per survey year)")
        ax.axhline(OUTBREAK, color=MUTED, linewidth=0.9, linestyle="--", label="outbreak level 0.22 COTS/tow")
        ax.set_title(reef, loc="left")
        style(ax)
        peaks += [{"reef": reef, **r} for r in peak_table(b["mean"], o, windows)]
        cm = cull_bands(ctraj[ctraj.reef_name == reef], "adults_ha", Q)["mean"].loc[YEARS[0]:YEARS[1]] if cull_bands else None
        annotate(ax, add_stats("fig1", reef, "LTMP manta COTS/tow", b["mean"].loc[YEARS[0]:YEARS[1]], o, cm,
                               threshold=OUTBREAK, log_offset=0.01), "upper left")
    for ax in axes[:, 0]:
        ax.set_ylabel("COTS per manta tow")
    h, lab = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, lab, loc="lower center", ncol=3, frameon=False, fontsize=7.5)
    fig.suptitle(f"Calibration reefs: COTS (manta-tow equivalent), {model_label}", x=0.01, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout(rect=(0, 0.09, 1, 1))
    save(fig, "fig1_cots_calibration_reefs")

    # Figure 2: 2020s densities against diver/cull targets.
    reefs2 = list(CALIB) + ["Martin Reef"]
    fig, axes = plt.subplots(1, len(reefs2), figsize=(17, 4.6), sharey=True)
    for ax, reef in zip(axes, reefs2):
        d = traj[traj.reef_name == reef] if reef in CALIB else extra[extra.reef_name == reef]
        b = bands(d, "adults_ha", base, arms).loc[2015:2026]
        ax.fill_between(b.index, b.olo.clip(lower=0.05), b.ohi, color=OUTER, linewidth=0, label=f"{PREFIX['']}±25% parameter envelope (full range)")
        ax.fill_between(b.index, b.lo.clip(lower=0.05), b.hi, color=INNER, linewidth=0, label=f"{PREFIX['']}seed range (5 seeds)")
        ax.plot(b.index, b["mean"].clip(lower=0.05), color=LINE, linewidth=2.0, label=f"{PREFIX['']}model (seed mean)")
        if cull_bands:
            cd = ctraj[ctraj.reef_name == reef] if reef in CALIB else cextra[cextra.reef_name == reef]
            draw_cull(ax, cull_bands(cd, "adults_ha"), lo_clip=0.05, years=(2015, 2026))
        t = targets[(targets.target_reef == reef) & (targets.target_density_ha > 0)]
        for src, mk in (("SALAD", "o"), ("cull_converted", "s")):
            tt = t[t.target_source == src]
            ax.scatter(tt.report_year, tt.target_density_ha, s=24, marker=mk, color=TARGET, zorder=5,
                       label="SALAD diver density" if src == "SALAD" else "cull catch converted to density")
        tv = t.groupby("report_year").target_density_ha.mean()
        cm = (cull_bands(ctraj[ctraj.reef_name == reef] if reef in CALIB else cextra[cextra.reef_name == reef],
                         "adults_ha")["mean"].loc[2015:2026] if cull_bands else None)
        annotate(ax, add_stats("fig2", reef, "diver/cull density (adults/ha)", b["mean"], tv, cm,
                               threshold=OUTBREAK / Q, log_offset=0.5, compact=True), "lower right")
        ax.set_yscale("log")
        ax.set_ylim(0.1, 300)
        ax.set_xlim(2014.5, 2026.5)
        ax.grid(color=GRID, linewidth=0.8)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(colors=MUTED, labelsize=8)
        ax.set_title(reef, loc="left")
    axes[0].set_ylabel("Adult COTS ha$^{-1}$ (log)")
    h, lab = axes[0].get_legend_handles_labels()
    fig.legend(h, lab, loc="lower center", ncol=4, frameon=False, fontsize=7.5)
    fig.suptitle("2020s build-up: model adults against diver and cull density estimates (statistics on the linear scale)",
                 x=0.01, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout(rect=(0, 0.14, 1, 1))
    save(fig, "fig2_cots_2020s_density")

    # Figure 3: additional reefs.
    fig, axes = plt.subplots(2, 3, figsize=(14, 7), sharex=True)
    for ax, (reef, oname) in zip(axes.flat, EXTRA.items()):
        d = extra[extra.reef_name == reef]
        b = bands(d, "adults_ha", base, arms, Q)
        draw(ax, b)
        if cull_bands:
            draw_cull(ax, cull_bands(cextra[cextra.reef_name == reef], "adults_ha", Q))
        if oname is not None:
            o = obs[obs.reef_name == oname].groupby("year").cotsptow.mean()
            ax.scatter(o.index, o.values, s=18, color=OBS, zorder=5, label="LTMP manta tow (mean per survey year)")
            peaks += [{"reef": reef, **r} for r in peak_table(b["mean"], o, windows)]
            cm = cull_bands(cextra[cextra.reef_name == reef], "adults_ha", Q)["mean"].loc[YEARS[0]:YEARS[1]] if cull_bands else None
            annotate(ax, add_stats("fig3", reef, "LTMP manta COTS/tow", b["mean"].loc[YEARS[0]:YEARS[1]], o, cm,
                                   threshold=OUTBREAK, log_offset=0.01), "upper left")
        ax.axhline(OUTBREAK, color=MUTED, linewidth=0.9, linestyle="--", label="outbreak level 0.22 COTS/tow")
        ax.set_title(f"{reef} — {EXTRA_NOTE[reef]}", loc="left")
        style(ax)
    for ax in axes[:, 0]:
        ax.set_ylabel("COTS per manta tow")
    h, lab = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, lab, loc="lower center", ncol=3, frameon=False, fontsize=7.5)
    fig.suptitle(f"Additional reefs (not used in fitting): COTS (manta-tow equivalent), {model_label}", x=0.01, ha="left",
                 fontsize=11, fontweight="bold")
    fig.tight_layout(rect=(0, 0.09, 1, 1))
    save(fig, "fig3_cots_additional_reefs")

    # Figure 4: coral cover at the calibration reefs.
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True, sharey=True)
    for ax, (reef, oname) in zip(axes.flat, CALIB.items()):
        d = traj[traj.reef_name == reef]
        draw(ax, bands(d, "coral_cover", base, arms))
        if cull_bands:
            draw_cull(ax, cull_bands(ctraj[ctraj.reef_name == reef], "coral_cover"))
        m = manta[manta.reef_name == oname].sort_values("report_year")
        ax.errorbar(m.report_year, m["median"], yerr=[m["median"] - m.lower, m.upper - m["median"]], fmt="o", ms=3.5,
                    color=OBS, ecolor=MUTED, elinewidth=0.8, capsize=0, zorder=5, label="LTMP manta hard coral (95% interval)")
        ph = photo[photo.domain_name == oname].groupby("report_year")["median"].mean()
        ax.scatter(ph.index, ph.values, s=22, marker="D", color=TARGET, zorder=6, label="LTMP 9 m photo-transect hard coral")
        cb = bands(d, "coral_cover", base, arms)["mean"].loc[YEARS[0]:YEARS[1]]
        ccb = cull_bands(ctraj[ctraj.reef_name == reef], "coral_cover")["mean"].loc[YEARS[0]:YEARS[1]] if cull_bands else None
        annotate(ax, add_stats("fig4", reef, "LTMP manta hard coral", cb, m.set_index("report_year")["median"], ccb, label="manta")
                 + add_stats("fig4", reef, "LTMP 9 m photo hard coral", cb, ph, ccb, label="photo"), "upper left")
        ax.set_title(reef, loc="left")
        ax.set_ylim(0, 0.8)
        style(ax)
    for ax in axes[:, 0]:
        ax.set_ylabel("Hard coral cover (fraction)")
    h, lab = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, lab, loc="lower center", ncol=3, frameon=False, fontsize=7.5)
    fig.suptitle("Calibration reefs: coral cover (model = fraction of habitable area)", x=0.01, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout(rect=(0, 0.09, 1, 1))
    save(fig, "fig4_coral_calibration_reefs")

    pk = pd.DataFrame(peaks)
    pk.to_csv(out / "peak_timing_summary.csv", index=False)
    files.append(out / "peak_timing_summary.csv")
    pd.DataFrame(stats).to_csv(out / "fit_statistics.csv", index=False)
    files.append(out / "fit_statistics.csv")
    meta = {"run": str(run.relative_to(ROOT)), "base_arms": base, "n_arms": len(arms),
            "cull_run": str(cull_run.relative_to(ROOT)) if cull_run is not None else None,
            "uncertainty": {"inner": "min-max over 5 base seed pairs", "outer": "min-max over base + 12 one-at-a-time +/-25% perturbations x 5 seeds (65 runs)"},
            "observation_model": f"simulated CPUE = {Q} x adults/ha (fixed_q0015)",
            "fit_statistics": "model seed mean vs observations in matching years: n, bias (model - obs), MAE, Pearson r, "
                              "Spearman rho (annotated); fit_statistics.csv adds relative bias %, Nash-Sutcliffe efficiency, "
                              "MAE of log(x + offset) (CPUE 0.01, density 0.5) and the true skill statistic for outbreak state "
                              "(0.22 COTS/tow; 14.7/ha)",
            "data_permissions": "Figure 2 shows values derived from unpublished SALAD and CSIRO cull data; check before publication",
            "inputs": {**{p.name: sha256(p) for p in (OBS_COTS, OBS_MANTA, OBS_PHOTO, TARGETS, run / "trajectories.csv", run / "trajectories_extra.csv")},
                       **({f"cull/{n}": sha256(cull_run / n) for n in ("trajectories.csv", "trajectories_extra.csv")} if cull_run is not None else {})},
            "script_sha256": sha256(Path(__file__)), "outputs": {p.name: sha256(p) for p in files}}
    (out / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(pk.round(2).to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--cull-run", type=Path, default=None)
    parser.add_argument("--model-label", default="W-wave", help="model name in figure titles (e.g. W-wave-lf04)")
    a = parser.parse_args()
    main(a.run.resolve(), a.output.resolve(), a.cull_run.resolve() if a.cull_run else None, a.model_label)
