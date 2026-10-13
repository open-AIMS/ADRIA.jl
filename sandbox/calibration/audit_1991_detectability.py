"""Audit the preregistered manta-detectability observation sensitivity
(crash_mechanism_protocol.md, 2026-10-13): stored run-11 and run-15 trajectories
re-scored under CPUE = q A A / (A + A_half), A_half in {0 (frozen), 1, 3, 10, 30}.
Requires the v2 re-score (raw and smoothed_3y observation treatments)."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from audit_1991_switching_repair import REEFS, TWO_WAVE, sha256  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
BASELINE = ("20261013T_11_cots_rebaseline_1991_v1", "wc_b150_q20")
BETA, Q = (100, 150), (10, 20, 30, 45, 70)
A_HALF = (0.0, 1.0, 3.0, 10.0, 30.0)
RAMP = {0.0: "#52514e", 1.0: "#86b6ef", 3.0: "#3987e5", 10.0: "#1c5cab", 30.0: "#0d366b"}


def main(run: Path) -> None:
    paths = {name: run / name for name in (
        "detectability_audit.csv", "detectability_audit.json", "detectability_trajectories.png")}
    if any(path.exists() for path in paths.values()):
        raise FileExistsError("Refusing to overwrite detectability audit")
    all_scores = pd.read_csv(run / "detectability_scores.csv")
    series = pd.read_csv(run / "detectability_series.csv")
    control = all_scores[all_scores.a_half == 0]
    replay = float((control.loss - control.stored_frozen_loss).abs().max())
    if replay != 0.0:
        raise ValueError(f"Control does not replay stored losses: {replay}")
    cots_obs = pd.read_csv(ROOT / "sandbox/data/reef_cots.csv")
    # Peaks, troughs and loss from the smoothed_3y score (as in every audit); matched
    # peaks also against the raw series, where the gate's 7 observed peaks are defined
    # (smoothing removes North Direction's 1995 peak, leaving 6).
    merge_key = ["run", "treatment", "reef_name", "a_half"]
    raw = all_scores[all_scores.observation_treatment == "raw"].set_index(merge_key)
    scores = all_scores[all_scores.observation_treatment == "smoothed_3y"].reset_index(drop=True)
    scores["matched_peaks_raw"] = raw.matched_peaks.reindex(
        pd.MultiIndex.from_frame(scores[merge_key])).to_numpy()
    scores["observed_peaks_smoothed"] = scores.obs_peak_years.astype(str).str.count(";") + 1
    scores["observed_peaks_raw"] = (raw.obs_peak_years.astype(str).str.count(";") + 1).reindex(
        pd.MultiIndex.from_frame(scores[merge_key])).to_numpy()

    rows = []
    for s in scores.itertuples(index=False):
        x = series[(series.run == s.run) & (series.treatment == s.treatment) &
                   (series.reef_name == s.reef_name) & (series.a_half == s.a_half) &
                   (series.year >= 1992)].sort_values("year").reset_index(drop=True)
        years = [int(v) for v in str(s.sim_peak_years).split(";") if v and v != "nan"]
        heights = [float(v) for v in str(s.sim_peak_heights_cpue).split(";") if v and v != "nan"]
        ratio = float("nan")
        if len(years) >= 2:
            smooth = x.detected_cpue.rolling(3, center=True, min_periods=1).mean()
            ratio = float(smooth[x.year.between(years[0], years[1])].min() / min(heights[:2]))
        cpue = x.set_index("year").detected_cpue
        obs = cots_obs[cots_obs.reef_name == REEFS[s.reef_name]].set_index("year").cotsptow
        rows.append({
            "run": s.run, "treatment": s.treatment, "reef_name": s.reef_name, "a_half": s.a_half,
            "loss": s.loss, "matched_peaks": s.matched_peaks, "matched_peaks_raw": s.matched_peaks_raw,
            "observed_peaks_smoothed": s.observed_peaks_smoothed, "observed_peaks_raw": s.observed_peaks_raw,
            "first_peak_year": years[0] if years else np.nan,
            "first_peak_cpue": heights[0] if heights else np.nan,
            "second_peak_year": years[1] if len(years) > 1 else np.nan,
            "second_peak_cpue": heights[1] if len(heights) > 1 else np.nan,
            "trough_to_smaller_peak": ratio,
            "mean_cpue_2000_2010": float(cpue.loc[2000:2010].mean()),
            "observed_mean_cpue_2000_2010": float(obs.loc[2000:2010].mean()),
            "cpue_1998": float(cpue.loc[1998]), "cpue_1999": float(cpue.loc[1999]),
            "observed_cpue_1998": float(obs.get(1998, np.nan)), "observed_cpue_1999": float(obs.get(1999, np.nan)),
        })
    result = pd.DataFrame(rows)
    result.to_csv(paths["detectability_audit.csv"], index=False)

    key = ["run", "treatment", "a_half"]
    summary = result.groupby(key).agg(
        summed_loss=("loss", "sum"), matched_smoothed=("matched_peaks", "sum"),
        observed_smoothed=("observed_peaks_smoothed", "sum"),
        matched_raw=("matched_peaks_raw", "sum"), observed_raw=("observed_peaks_raw", "sum"),
        reefs_with_first_peak=("first_peak_year", "count"),
        mean_first_peak_cpue=("first_peak_cpue", "mean"),
        mean_cpue_2000_2010=("mean_cpue_2000_2010", "mean"),
        observed_mean_cpue_2000_2010=("observed_mean_cpue_2000_2010", "mean"))
    gate = result[result.reef_name.isin(TWO_WAVE)].groupby(key).agg(
        two_wave_reefs_with_two_peaks=("second_peak_year", "count"),
        worst_trough_ratio=("trough_to_smaller_peak", "max"))
    summary = summary.join(gate)
    two_wave = summary.two_wave_reefs_with_two_peaks == len(TWO_WAVE)
    trough = summary.worst_trough_ratio <= .10
    # As implemented in earlier audits (smoothed count against 7: unattainable, max 6),
    # and as the protocol states it (all observed peaks; 7 on the raw series).
    summary["gate_as_audited_smoothed_vs_7"] = two_wave & (summary.matched_smoothed == 7) & trough
    summary["gate_raw_all_observed_matched"] = two_wave & (summary.matched_raw == summary.observed_raw) & trough
    for reef, short in (("Lizard Island Reef", "liz"), ("MacGillivray Reef", "macg"),
                        ("North Direction Reef", "nd"), ("Eyrie Reef", "eyr")):
        r = result[result.reef_name == reef].set_index(key)
        summary[f"{short}_peaks"] = (r.first_peak_cpue.round(2).astype(str) + " / " +
                                     r.second_peak_cpue.round(2).astype(str))
        summary[f"{short}_trough_ratio"] = r.trough_to_smaller_peak

    # Q re-selection under each A_half with the frozen run-11 rule (per beta).
    q_selection = {}
    run11 = summary.loc[BASELINE[0]]
    for a_half in A_HALF:
        for b in BETA:
            arms = run11.xs(a_half, level="a_half").loc[[f"wc_b{b}_q{q}" for q in Q]]
            eligible = arms[arms.reefs_with_first_peak == len(REEFS)].copy()
            eligible["q"] = [int(i.split("_q")[1]) for i in eligible.index]
            best = eligible.sort_values(["summed_loss", "q"]).index[0] if len(eligible) else None
            q_selection[f"a_half_{a_half:g}_beta_{b}"] = {
                "selected": best, "summed_loss": float(arms.loc[best, "summed_loss"]) if best else None,
                "q_at_grid_boundary": bool(best and int(best.split("_q")[1]) in (min(Q), max(Q)))}

    fig, axes = plt.subplots(4, 2, figsize=(14, 13), sharex=True)
    for index, (reef, observed_name) in enumerate(REEFS.items()):
        lin, log = axes[index]
        for a_half in A_HALF:
            x = series[(series.run == BASELINE[0]) & (series.treatment == BASELINE[1]) &
                       (series.reef_name == reef) & (series.a_half == a_half) & (series.year >= 1992)]
            label = ("frozen linear (A½ = 0)" if a_half == 0 else f"A½ = {a_half:g} /ha") if index == 0 else None
            style = dict(color=RAMP[a_half], linewidth=2, linestyle="--" if a_half == 0 else "-", label=label)
            lin.plot(x.year, x.detected_cpue, **style)
            log.plot(x.year, x.detected_cpue.clip(lower=1e-4), **style)
        obs = cots_obs[(cots_obs.reef_name == observed_name) & cots_obs.year.between(1992, 2024)]
        for axis, values in ((lin, obs.cotsptow), (log, obs.cotsptow.clip(lower=1e-3))):
            axis.scatter(obs.year, values, s=18, color="#0b0b0b", zorder=5,
                         label="LTMP COTS/tow (zeros at 0.001 on log)" if index == 0 and axis is lin else None)
            axis.axhline(.22, color="#52514e", linestyle=":", linewidth=.8)
            axis.grid(color="#e4e3df", linewidth=.8)
            axis.set_xlim(1991, 2024)
        lin.set_ylim(bottom=0)
        log.set_yscale("log")
        log.set_ylim(1e-4, 3)
        lin.set_ylabel(f"{reef}\nCOTS/tow")
    axes[0, 0].set_title(f"Baseline {BASELINE[1]} under each observation model (linear axis)", loc="left", fontsize=10)
    axes[0, 1].set_title("Same, log axis: troughs", loc="left", fontsize=10)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=6, fontsize=8, frameon=False, bbox_to_anchor=(.5, -.01))
    fig.tight_layout(rect=(0, .04, 1, 1))
    fig.savefig(paths["detectability_trajectories.png"], dpi=160, bbox_inches="tight", facecolor="#fcfcfb")
    plt.close(fig)

    flat = summary.reset_index()
    report = {
        "status": "observation_sensitivity_not_promoted",
        "control_replay_max_abs_loss_delta": replay,
        "q_selection_under_each_a_half": q_selection,
        "summary": json.loads(flat.to_json(orient="records")),
        "criteria": "2 peaks on each of Lizard/MacGillivray/North Direction; all observed peaks matched "
                    "(7 on the raw series; the smoothed series has 6); worst trough ratio <=0.10",
        "inputs": {name: sha256(run / name) for name in (
            "metadata.toml", "detectability_scores.csv", "detectability_series.csv")},
        "script_sha256": sha256(Path(__file__)),
        "outputs": {name: sha256(path) for name, path in paths.items() if name != "detectability_audit.json"},
    }
    paths["detectability_audit.json"].write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    pd.set_option("display.width", 260)
    pd.set_option("display.max_columns", 30)
    focus = flat[flat.treatment.isin([BASELINE[1], "wc_b100_q20", "acrpoc_g0"])]
    print(focus.drop(columns=["reefs_with_first_peak"] + [c for c in focus if c.endswith("_peaks")])
          .round(3).to_string(index=False))
    print("Q selection:", json.dumps(q_selection, indent=1))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    main(parser.parse_args().run.resolve())
