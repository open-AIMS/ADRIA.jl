"""Audit the preregistered COTS re-baseline on the working coral configuration
(crash_mechanism_protocol.md, 2026-10-13): run-15 leading candidate with COTS dynamics
on, CoralBlox growth untouched, Mellin + 3%/yr, cyclone kappa 5, bleaching beta {1, 1.5},
predator-pit Q re-checked; ref and single-component arms. COTS metrics as in the
food-grace audit; coral J as in the coral calibration (frozen)."""

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
from audit_1991_coral_calibration import observed_photo  # noqa: E402
from audit_1991_growth_starvation import PHOTO, PHOTO_REEFS, slope  # noqa: E402
from audit_1991_initiation import exact_replay  # noqa: E402
from audit_1991_joint_coral import objective, reef_errors  # noqa: E402
from audit_1991_switching_repair import REEFS, TWO_WAVE, sha256  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
PRIOR = ROOT / "sandbox/calibration/runs/20261012T_15_food_grace_1991_v1"
REFERENCE, PRIOR_REFERENCE = "ref_q30", "acrpoc_g0"
BETA, Q = (100, 150), (10, 20, 30, 45, 70)
WORKING = [f"wc_b{b}_q{q}" for b in BETA for q in Q]
ARMS = [REFERENCE, "k500_q30", "mel3_q30"] + WORKING
OBSERVED_FIRST = {"Lizard Island Reef": 1.92, "MacGillivray Reef": 0.86}
THETA = 0.05
GAP = (2004, 2010)
SECOND_WINDOW = (2010, 2017)
N_YEARS = 34


def main(run: Path) -> None:
    paths = {name: run / name for name in (
        "cots_rebaseline_audit.csv", "cots_rebaseline_audit.json", "cots_rebaseline.png")}
    if any(path.exists() for path in paths.values()):
        raise FileExistsError("Refusing to overwrite re-baseline audit")
    tables = {name: pd.read_csv(run / f"{name}.csv")
              for name in ("trajectories", "flows", "mechanisms", "scores", "failures")}
    if len(tables["failures"]):
        raise ValueError("Failed evaluations in run")
    trajectories, flows, mechanisms = (tables[k] for k in ("trajectories", "flows", "mechanisms"))
    if set(trajectories.treatment) != set(ARMS):
        raise ValueError("Missing or extra treatment")
    if any(len(t) != len(ARMS) * len(REEFS) * N_YEARS for t in (trajectories, flows, mechanisms)):
        raise ValueError("Incomplete reef-year outputs")
    replay = {}
    for name in ("trajectories", "flows", "mechanisms"):
        prior = pd.read_csv(PRIOR / f"{name}.csv")
        columns = [c for c in prior if c not in ("treatment", "reef_name", "year")]
        replay[name] = exact_replay(tables[name][tables[name].treatment == REFERENCE],
                                    prior[prior.treatment == PRIOR_REFERENCE], columns)
    prior_flows = pd.read_csv(PRIOR / "flows.csv")
    reference = prior_flows[prior_flows.treatment == PRIOR_REFERENCE].sort_values(["reef_name", "year"])
    external = {}
    for arm in ARMS:
        current = flows[flows.treatment == arm].sort_values(["reef_name", "year"])
        external[arm] = float(np.max(np.abs(current.external_immigration.to_numpy() -
                                            reference.external_immigration.to_numpy())))
        if external[arm] > 1e-9:
            raise ValueError(f"External immigration differs from Ext x5 in {arm}")
    numeric = mechanisms[["recruits_ha", "subadults_ha", "adults_ha", "fast_coral", "total_coral",
                          "food_survival", "body_condition", "juvenile_predation_ha"]].to_numpy()
    if not np.isfinite(numeric).all() or (numeric < -1e-12).any():
        raise ValueError("Nonfinite or negative state/diagnostic")

    photo = observed_photo()
    cots_obs = pd.read_csv(ROOT / "sandbox/data/reef_cots.csv")
    manta_raw = pd.read_csv(ROOT / "sandbox/data/reef_manta.csv")
    manta_raw = manta_raw[(manta_raw.data_type == "manta") & (manta_raw.domain_category == "reef") &
                          (manta_raw.variable == "HC") & (manta_raw.purpose == "MANTA") &
                          (manta_raw.project_code == "LTMP") & (manta_raw.depth == 9)]
    rows = []
    for score in tables["scores"][tables["scores"].observation_treatment == "smoothed_3y"].itertuples(index=False):
        series = trajectories[(trajectories.treatment == score.treatment) &
                              (trajectories.reef_name == score.reef_name) &
                              (trajectories.year >= 1992)].sort_values("year")
        mech = mechanisms[(mechanisms.treatment == score.treatment) &
                          (mechanisms.reef_name == score.reef_name)].sort_values("year").set_index("year")
        mech = mech.assign(slow_coral=mech.total_coral - mech.fast_coral)
        obs = cots_obs[cots_obs.reef_name == REEFS[score.reef_name]]
        years = [int(v) for v in str(score.sim_peak_years).split(";") if v and v != "nan"]
        heights = [float(v) for v in str(score.sim_peak_heights_cpue).split(";") if v and v != "nan"]
        ratio = float("nan")
        if len(years) >= 2:
            smooth = series.simulated_cpue.rolling(3, center=True, min_periods=1).mean()
            ratio = float(smooth[series.year.between(years[0], years[1])].min() / min(heights[:2]))
        cpue = series[series.year.between(1995, 2006)].simulated_cpue.to_numpy()
        steps = cpue[1:] / np.maximum(cpue[:-1], 1e-12)
        manta = manta_raw[(manta_raw.reef_name == REEFS[score.reef_name]) &
                          manta_raw.report_year.between(1992, 2024)].groupby("report_year")["median"].mean()
        arm = score.treatment
        row = {
            "treatment": arm, "reef_name": score.reef_name,
            "coral": arm.split("_")[0] if not arm.startswith("wc_") else "wc",
            "bleaching_scalar": int(arm.split("_")[1][1:]) / 100 if arm.startswith("wc_") else 1.0,
            "q_cots_ha_yr": int(arm.split("_q")[-1]),
            "first_peak_year": years[0] if years else float("nan"),
            "second_peak_year": years[1] if len(years) > 1 else float("nan"),
            "first_peak_cpue": heights[0] if heights else float("nan"),
            "second_peak_cpue": heights[1] if len(heights) > 1 else float("nan"),
            "max_cpue_1995_2000": float(series[series.year.between(1995, 2000)].simulated_cpue.max()),
            "observed_max_cpue_1995_2000": float(obs[obs.year.between(1995, 2000)].cotsptow.max()),
            "max_cpue_2010_2017": float(series[series.year.between(*SECOND_WINDOW)].simulated_cpue.max()),
            "observed_max_cpue_2010_2017": float(obs[obs.year.between(*SECOND_WINDOW)].cotsptow.max()),
            "trough_to_smaller_peak": ratio,
            "matched_peaks": int(score.matched_peaks),
            "smoothed_loss": float(score.loss),
            "steepest_one_year_drop_ratio_1995_2006": float(steps.min()),
            "min_smoothed_adults_1999_2012": float(mech.loc[1999:2012].adults_ha.rolling(3, center=True, min_periods=1).mean().min()),
            "mean_gap_adults_2004_2010": float(mech.loc[GAP[0]:GAP[1]].adults_ha.mean()),
            "max_subadults_ha": float(mech.loc[1992:2024].subadults_ha.max()),
            "juvenile_predation_ha_1992_2024": float(mech.loc[1992:2024].juvenile_predation_ha.sum()),
            "mean_fast_coral_1992_1996": float(mech.loc[1992:1996].fast_coral.mean()),
            "min_fast_coral_1998_2003": float(mech.loc[1998:2003].fast_coral.min()),
            "mean_fast_coral_2000_2013": float(mech.loc[2000:2013].fast_coral.mean()),
            "years_fast_above_theta_2000_2013": int((mech.loc[2000:2013].fast_coral > THETA).sum()),
            "slow_coral_2003": float(mech.loc[2003, "slow_coral"]),
            "slow_coral_2012": float(mech.loc[2012, "slow_coral"]),
        }
        row |= reef_errors(mech, score.reef_name, photo, manta, (1992, 2024))
        if score.reef_name in PHOTO_REEFS:
            p = photo[photo.model_reef == score.reef_name].set_index("year")
            row |= {"observed_mean_fast_2000_2013": float(p.loc[2000:2013].obs_fast.mean()),
                    "model_fast_slope_2017_2024": slope(range(2017, 2025), mech.fast_coral.loc[2017:2024]),
                    "observed_fast_slope_2017_2024": slope(p.loc[2017:2024].index, p.loc[2017:2024].obs_fast)}
        rows.append(row)
    result = pd.DataFrame(rows)
    result["order"] = result.treatment.map({a: i for i, a in enumerate(ARMS)})
    result = result.sort_values(["order", "reef_name"]).drop(columns="order")
    result.to_csv(paths["cots_rebaseline_audit.csv"], index=False)

    two_wave = result[result.reef_name.isin(TWO_WAVE)]
    summary = result.groupby("treatment").agg(
        bleaching_scalar=("bleaching_scalar", "first"), q_cots_ha_yr=("q_cots_ha_yr", "first"),
        summed_loss=("smoothed_loss", "sum"),
        reefs_with_first_peak=("first_peak_year", "count"),
        reefs_with_two_peaks=("second_peak_year", "count"),
        matched_observed_peaks_of_seven=("matched_peaks", "sum"),
        mean_first_peak_year=("first_peak_year", "mean"),
        mean_first_peak_cpue=("first_peak_cpue", "mean"),
        mean_min_trough_adults_ha=("min_smoothed_adults_1999_2012", "mean"),
        mean_max_cpue_2010_2017=("max_cpue_2010_2017", "mean"),
        mean_fast_2000_2013=("mean_fast_coral_2000_2013", "mean"),
        mean_slow_2012=("slow_coral_2012", "mean"),
        mean_juvenile_predation_ha=("juvenile_predation_ha_1992_2024", "mean"),
    ).reindex(ARMS)
    for reef, short in (("Lizard Island Reef", "liz"), ("MacGillivray Reef", "macg")):
        r = result[result.reef_name == reef].set_index("treatment")
        summary[f"{short}_first_peak_cpue"] = r.first_peak_cpue
        summary[f"{short}_max_2010_2017"] = r.max_cpue_2010_2017
    gate = two_wave.groupby("treatment").agg(
        two_wave_reefs_with_two_peaks=("second_peak_year", "count"),
        worst_available_trough_ratio=("trough_to_smaller_peak", "max"))
    summary = summary.join(gate)
    summary["passes_qualitative_gate"] = (
        (summary.two_wave_reefs_with_two_peaks == len(TWO_WAVE)) &
        (summary.matched_observed_peaks_of_seven == 7) &
        (summary.worst_available_trough_ratio <= .10))
    summary["coral_J"] = objective(result).reindex(ARMS)

    # Frozen selection: Q per beta (min summed loss, first peak on every reef, ties to
    # smaller Q; bounds flagged), then beta by the lower coral J at its selected Q.
    selection = {}
    for b in BETA:
        arms = summary.loc[[f"wc_b{b}_q{q}" for q in Q]]
        eligible = arms[arms.reefs_with_first_peak == len(REEFS)]
        if eligible.empty:
            selection[str(b)] = {"selected": None, "reason": "no Q gives a first peak on every reef"}
            continue
        best = eligible.sort_values(["summed_loss", "q_cots_ha_yr"]).index[0]
        selection[str(b)] = {"selected": best, "summed_loss": float(arms.loc[best, "summed_loss"]),
                             "coral_J": float(arms.loc[best, "coral_J"]),
                             "q_at_grid_boundary": bool(arms.loc[best, "q_cots_ha_yr"] in (min(Q), max(Q)))}
    chosen = [v for v in selection.values() if v["selected"]]
    baseline = min(chosen, key=lambda v: v["coral_J"])["selected"] if chosen else None

    selected = [v["selected"] for v in selection.values() if v["selected"]]
    shown = {REFERENCE: ("#475569", "-", "run-15 candidate (ref)"),
             "k500_q30": ("#d97706", ":", "kappa 5 only, Q 30"),
             "mel3_q30": ("#7c3aed", ":", "Mellin + 3% only, Q 30")}
    for arm, color in zip(selected, ("#0d9488", "#dc2626")):
        shown[arm] = (color, "-", f"working coral {arm} (selected Q)")
    fig, axes = plt.subplots(4, 3, figsize=(17, 13), sharex=True)
    for index, (reef, observed_name) in enumerate(REEFS.items()):
        cots_ax, fast_ax, slow_ax = axes[index]
        for arm, (color, style, label) in shown.items():
            values = trajectories[(trajectories.treatment == arm) &
                                  (trajectories.reef_name == reef)].sort_values("year")
            mech = mechanisms[(mechanisms.treatment == arm) & (mechanisms.reef_name == reef)].sort_values("year")
            kwargs = dict(color=color, linestyle=style, linewidth=1.6, label=label if index == 0 else None)
            cots_ax.plot(values.year.iloc[1:], values.simulated_cpue.iloc[1:], **kwargs)
            fast_ax.plot(mech.year, mech.fast_coral, **kwargs)
            slow_ax.plot(mech.year, mech.total_coral - mech.fast_coral, **kwargs)
        obs_cots = cots_obs[(cots_obs.reef_name == observed_name) & cots_obs.year.between(1992, 2024)]
        cots_ax.scatter(obs_cots.year, obs_cots.cotsptow, s=16, color="#111827", alpha=.7, zorder=5,
                        label="LTMP observations" if index == 0 else None)
        if reef in PHOTO_REEFS:
            p = photo[photo.model_reef == reef]
            fast_ax.scatter(p.year, p.obs_fast, s=16, color="#111827", alpha=.7, zorder=5)
            slow_ax.scatter(p.year, p.obs_slow, s=16, color="#111827", alpha=.7, zorder=5)
        cots_ax.set_ylabel(f"{reef}\nCOTS/tow (q = 0.015)")
        fast_ax.set_ylabel("fast coral (habitable)")
        slow_ax.set_ylabel("slow coral (habitable)")
        cots_ax.axhline(.22, color="#92400e", linewidth=.8, linestyle=":")
        fast_ax.axhline(THETA, color="#92400e", linewidth=.8, linestyle=":")
        for axis in (cots_ax, fast_ax, slow_ax):
            axis.set_xlim(1991, 2024)
            axis.set_ylim(bottom=0)
            axis.grid(alpha=.15)
    axes[0, 0].set_title("COTS (dotted: 0.22 COTS/tow outbreak level)")
    axes[0, 1].set_title(f"Fast coral vs photo transects (dotted: theta {THETA})")
    axes[0, 2].set_title("Slow coral vs photo transects")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, fontsize=9, bbox_to_anchor=(.5, -.02))
    fig.tight_layout(rect=(0, .05, 1, 1))
    fig.savefig(paths["cots_rebaseline.png"], dpi=160, bbox_inches="tight")
    plt.close(fig)

    report = {
        "status": "bounded_diagnostic_not_promoted",
        "reference_replay_max_rel_delta": replay,
        "external_immigration_max_abs_delta_vs_reference": external,
        "q_selection_by_beta": selection,
        "baseline_for_next_screens": baseline,
        "treatments": json.loads(summary.to_json(orient="index")),
        "criteria": "2 peaks on each of Lizard/MacGillivray/North Direction; 7/7 matched; "
                    "worst trough ratio <=0.10; absolute CPUE under q=0.015",
        "observed_first_peaks_cots_tow": OBSERVED_FIRST,
        "inputs": {name: sha256(run / name) for name in (
            "metadata.toml", "scores.csv", "trajectories.csv", "flows.csv", "mechanisms.csv")} |
                  {"reef_photo_transect.csv": sha256(PHOTO)},
        "script_sha256": sha256(Path(__file__)),
        "outputs": {name: sha256(path) for name, path in paths.items()
                    if name != "cots_rebaseline_audit.json"},
    }
    paths["cots_rebaseline_audit.json"].write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 40)
    print(summary.round(3).to_string())
    print("Q selection:", json.dumps(selection, indent=1))
    print("Baseline for next screens:", baseline)
    print("Reference replay max delta:", max(replay.values()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    main(parser.parse_args().run.resolve())
