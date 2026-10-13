"""Audit the preregistered reef-differentiated 1991 starting state
(crash_mechanism_protocol.md, 2026-10-13): CPUE surface {idw, eb} x inversion {lin,
h10} x Q {20, 30} on the wc_b150 baseline. Needs the run and its detectability re-score
(rescore_detectability.jl over all arms). Each arm is read under the frozen linear
observation model and under A_half = 10; the 'matching' reading pairs lin with A_half 0
and h10 with A_half 10."""

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
from audit_1991_initiation import exact_replay  # noqa: E402
from audit_1991_switching_repair import REEFS, TWO_WAVE, sha256  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
PRIOR = ROOT / "sandbox/calibration/runs/20261013T_11_cots_rebaseline_1991_v1"
ARMS = [f"wc_b150_q{q}_{s}_{i}" for q in (20, 30) for s in ("idw", "eb") for i in ("lin", "h10")]
REFERENCE, PRIOR_REFERENCE = "wc_b150_q20_idw_lin", "wc_b150_q20"
OBSERVED_FIRST = {"Lizard Island Reef": 1.92, "MacGillivray Reef": 0.86,
                  "North Direction Reef": 0.56, "Eyrie Reef": 0.13}
COLORS = {"idw_lin": "#52514e", "eb_lin": "#eb6834", "idw_h10": "#2a78d6", "eb_h10": "#1baf7a"}
N_YEARS = 34


def rank_agreement(model: dict, observed: dict) -> float:
    reefs = list(observed)
    m = pd.Series([model[r] for r in reefs]).rank()
    o = pd.Series([observed[r] for r in reefs]).rank()
    return float(np.corrcoef(m, o)[0, 1])


def main(run: Path, rescore: Path) -> None:
    paths = {name: run / name for name in (
        "initial_state_audit.csv", "initial_state_audit.json", "initial_state_rebaseline.png")}
    if any(path.exists() for path in paths.values()):
        raise FileExistsError("Refusing to overwrite initial-state audit")
    tables = {name: pd.read_csv(run / f"{name}.csv")
              for name in ("trajectories", "flows", "mechanisms", "failures", "initial_reef")}
    if len(tables["failures"]):
        raise ValueError("Failed evaluations in run")
    if set(tables["trajectories"].treatment) != set(ARMS) or len(tables["trajectories"]) != len(ARMS) * 4 * N_YEARS:
        raise ValueError("Missing, extra or incomplete treatments")
    replay = {}
    for name in ("trajectories", "flows", "mechanisms"):
        prior = pd.read_csv(PRIOR / f"{name}.csv")
        columns = [c for c in prior if c not in ("treatment", "reef_name", "year")]
        replay[name] = exact_replay(tables[name][tables[name].treatment == REFERENCE],
                                    prior[prior.treatment == PRIOR_REFERENCE], columns)
    flows = tables["flows"]
    reference = flows[flows.treatment == REFERENCE].sort_values(["reef_name", "year"])
    external = {arm: float(np.max(np.abs(flows[flows.treatment == arm].sort_values(["reef_name", "year"])
                                         .external_immigration.to_numpy() - reference.external_immigration.to_numpy())))
                for arm in ARMS}
    if max(external.values()) > 1e-9:
        raise ValueError("External immigration differs between arms")
    mech = tables["mechanisms"]
    if not np.isfinite(mech[["recruits_ha", "subadults_ha", "adults_ha", "fast_coral"]].to_numpy()).all():
        raise ValueError("Nonfinite state")

    scores = pd.read_csv(rescore / "detectability_scores.csv")
    series = pd.read_csv(rescore / "detectability_series.csv")
    scores = scores[scores.run == run.name]
    series = series[series.run == run.name]
    if set(scores.treatment) != set(ARMS):
        raise ValueError("Re-score does not cover every arm of this run")
    control = scores[scores.a_half == 0]
    if float((control.loss - control.stored_frozen_loss).abs().max()) != 0.0:
        raise ValueError("Re-score control does not reproduce stored losses")
    seeds = tables["initial_reef"].set_index(["treatment", "reef_name"]).seed_adults_ha
    cots_obs = pd.read_csv(ROOT / "sandbox/data/reef_cots.csv")

    rows = []
    smooth = scores[scores.observation_treatment == "smoothed_3y"]
    raw = scores[scores.observation_treatment == "raw"].set_index(["treatment", "reef_name", "a_half"])
    for s in smooth[smooth.a_half.isin([0.0, 10.0])].itertuples(index=False):
        x = series[(series.treatment == s.treatment) & (series.reef_name == s.reef_name) &
                   (series.a_half == s.a_half) & (series.year >= 1992)].sort_values("year").reset_index(drop=True)
        years = [int(v) for v in str(s.sim_peak_years).split(";") if v and v != "nan"]
        heights = [float(v) for v in str(s.sim_peak_heights_cpue).split(";") if v and v != "nan"]
        ratio = float("nan")
        if len(years) >= 2:
            sm = x.detected_cpue.rolling(3, center=True, min_periods=1).mean()
            ratio = float(sm[x.year.between(years[0], years[1])].min() / min(heights[:2]))
        r = raw.loc[(s.treatment, s.reef_name, s.a_half)]
        q, source, inversion = s.treatment.split("_")[2:5]
        rows.append({
            "treatment": s.treatment, "reef_name": s.reef_name, "a_half_observation": s.a_half,
            "q": int(q[1:]), "surface": source, "inversion": inversion,
            "matching_observation": (inversion == "lin") == (s.a_half == 0.0),
            "seed_adults_ha_1991": float(seeds.loc[(s.treatment, s.reef_name)]),
            "loss": s.loss, "matched_smoothed": s.matched_peaks, "matched_raw": int(r.matched_peaks),
            "observed_raw": str(r.obs_peak_years).count(";") + 1,
            "first_peak_year": years[0] if years else np.nan,
            "first_peak_cpue": heights[0] if heights else np.nan,
            "second_peak_cpue": heights[1] if len(heights) > 1 else np.nan,
            "max_cpue_1995_2000": float(x[x.year.between(1995, 2000)].detected_cpue.max()),
            "trough_to_smaller_peak": ratio,
        })
    result = pd.DataFrame(rows)
    result.to_csv(paths["initial_state_audit.csv"], index=False)

    key = ["treatment", "a_half_observation"]
    summary = result.groupby(key).agg(
        matching_observation=("matching_observation", "first"),
        summed_loss=("loss", "sum"), matched_raw=("matched_raw", "sum"), observed_raw=("observed_raw", "sum"),
        matched_smoothed=("matched_smoothed", "sum"), mean_first_peak_year=("first_peak_year", "mean"))
    gate = result[result.reef_name.isin(TWO_WAVE)].groupby(key).agg(
        two_wave_reefs=("second_peak_cpue", "count"), worst_trough_ratio=("trough_to_smaller_peak", "max"))
    summary = summary.join(gate)
    summary["gate_raw_all_observed_matched"] = ((summary.two_wave_reefs == len(TWO_WAVE)) &
                                                (summary.matched_raw == summary.observed_raw) &
                                                (summary.worst_trough_ratio <= .10))
    for reef, short in (("Lizard Island Reef", "liz"), ("MacGillivray Reef", "macg"),
                        ("North Direction Reef", "nd"), ("Eyrie Reef", "eyr")):
        r = result[result.reef_name == reef].set_index(key)
        summary[f"{short}_seed"] = r.seed_adults_ha_1991
        summary[f"{short}_first_1995_2000"] = r.max_cpue_1995_2000
    ranking = {}
    for (arm, a_half), group in result.groupby(key):
        model = dict(zip(group.reef_name, group.max_cpue_1995_2000))
        ranking[f"{arm}|A_half={a_half:g}"] = rank_agreement(model, OBSERVED_FIRST)
    summary["first_wave_rank_corr_vs_observed"] = [ranking[f"{a}|A_half={h:g}"] for a, h in summary.index]

    fig, axes = plt.subplots(4, 1, figsize=(11, 13), sharex=True)
    for index, (reef, observed_name) in enumerate(REEFS.items()):
        axis = axes[index]
        for combo, color in COLORS.items():
            arm = f"wc_b150_q20_{combo}"
            a_half = 0.0 if combo.endswith("lin") else 10.0
            x = series[(series.treatment == arm) & (series.reef_name == reef) &
                       (series.a_half == a_half) & (series.year >= 1992)]
            label = f"{combo.replace('_', ' / ')} (observed with A½ = {a_half:g})" if index == 0 else None
            axis.plot(x.year, x.detected_cpue, color=color, linewidth=2,
                      linestyle="--" if combo == "idw_lin" else "-", label=label)
        obs = cots_obs[(cots_obs.reef_name == observed_name) & cots_obs.year.between(1992, 2024)]
        axis.scatter(obs.year, obs.cotsptow, s=18, color="#0b0b0b", zorder=5,
                     label="LTMP COTS/tow" if index == 0 else None)
        axis.axhline(.22, color="#52514e", linestyle=":", linewidth=.8)
        axis.set_ylabel(f"{reef}\nCOTS/tow")
        axis.set_ylim(bottom=0)
        axis.set_xlim(1991, 2024)
        axis.grid(color="#e4e3df", linewidth=.8)
    axes[0].set_title("1991 starting state (surface / inversion), Q 20, each read with its matching observation model",
                      loc="left", fontsize=10)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, fontsize=8, frameon=False, bbox_to_anchor=(.5, -.01))
    fig.tight_layout(rect=(0, .05, 1, 1))
    fig.savefig(paths["initial_state_rebaseline.png"], dpi=160, bbox_inches="tight", facecolor="#fcfcfb")
    plt.close(fig)

    flat = summary.reset_index()
    report = {
        "status": "bounded_diagnostic_not_promoted",
        "reference_replay_max_rel_delta": replay,
        "external_immigration_max_abs_delta_vs_reference": external,
        "observed_first_wave_cots_tow": OBSERVED_FIRST,
        "summary": json.loads(flat.to_json(orient="records")),
        "inputs": {name: sha256(run / name) for name in (
            "metadata.toml", "trajectories.csv", "flows.csv", "mechanisms.csv", "initial_reef.csv")} |
                  {f"rescore/{name}": sha256(rescore / name) for name in (
                      "metadata.toml", "detectability_scores.csv", "detectability_series.csv")},
        "script_sha256": sha256(Path(__file__)),
        "outputs": {name: sha256(path) for name, path in paths.items() if name != "initial_state_audit.json"},
    }
    paths["initial_state_audit.json"].write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    pd.set_option("display.width", 260)
    pd.set_option("display.max_columns", 30)
    print(flat.round(3).to_string(index=False))
    print("Reference replay max delta:", max(replay.values()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("rescore", type=Path)
    args = parser.parse_args()
    main(args.run.resolve(), args.rescore.resolve())
