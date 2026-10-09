"""Audit the four-arm 1991 preferred-prey starvation x adult senescence screen.

Checks the exact archived-control replay, invariant external supply and mechanism
diagnostics, then reports the preregistered peak/trough gate and plots COTS/coral.
See crash_mechanism_protocol.md (2026-10-09 structural comparison).
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

ROOT = Path(__file__).resolve().parents[2]
ARCHIVE = ROOT / "sandbox/calibration/runs/20261007T199111_owen_1991_cohort_history_final"
CONTROL = "cohort_full_control"
ARMS = {
    CONTROL: ("1991 full-cohort control", "#475569", "-"),
    "preferred_starvation": ("Preferred-prey starvation", "#0d9488", "-"),
    "senescence": ("Adult senescence (6+)", "#dc2626", "-"),
    "preferred_starvation_senescence": ("Both", "#7c3aed", "--"),
}
STARVATION_ARMS = {"preferred_starvation", "preferred_starvation_senescence"}
SENESCENCE_ARMS = {"senescence", "preferred_starvation_senescence"}
REEFS = {
    "Lizard Island Reef": "Lizard Isles",
    "MacGillivray Reef": "Macgillivray Reef",
    "North Direction Reef": "North Direction Island",
    "Eyrie Reef": "Eyrie Reef",
}
TWO_WAVE = set(list(REEFS)[:3])
GAP = (2004, 2012)
N_YEARS = 34


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def exact_replay(current: pd.DataFrame, old: pd.DataFrame,
                 columns: list[str]) -> dict[str, float]:
    keys = ["reef_name", "year"]
    joined = current.merge(old, on=keys, suffixes=("_new", "_old"), validate="one_to_one")
    if len(joined) != len(current) or len(joined) != len(old):
        raise ValueError("Archived control does not cover the same reef-years")
    deltas = {column: float(np.max(np.abs(joined[f"{column}_new"] -
                                          joined[f"{column}_old"]))) for column in columns}
    if max(deltas.values()) > 1e-10:
        raise ValueError(f"Archived control did not replay: {deltas}")
    return deltas


def main(run: Path) -> None:
    paths = {name: run / name for name in (
        "structural_crash_audit.csv", "structural_crash_audit.json",
        "cots_coral_structural_crash.png")}
    if any(path.exists() for path in paths.values()):
        raise FileExistsError("Refusing to overwrite structural-crash audit")
    trajectories = pd.read_csv(run / "trajectories.csv")
    flows = pd.read_csv(run / "flows.csv")
    mechanisms = pd.read_csv(run / "mechanisms.csv")
    scores = pd.read_csv(run / "scores.csv")
    failures = pd.read_csv(run / "failures.csv")
    if len(failures):
        raise ValueError("Failed evaluations in run")
    if set(trajectories.treatment) != set(ARMS):
        raise ValueError("Missing or extra treatment")
    if any(len(table) != len(ARMS) * len(REEFS) * N_YEARS
           for table in (trajectories, flows, mechanisms)):
        raise ValueError("Incomplete reef-year outputs")
    if sha256(run / "initial_state_alpha015_cohort_full.csv") != sha256(
            ARCHIVE / "initial_state_alpha015_cohort_full.csv"):
        raise ValueError("Initial state differs from archived cohort anchor")

    old_traj = pd.read_csv(ARCHIVE / "trajectories.csv")
    old_flows = pd.read_csv(ARCHIVE / "flows.csv")
    replay_trajectories = exact_replay(
        trajectories[trajectories.treatment == CONTROL],
        old_traj[old_traj.treatment == "idw_alpha015_cohort_full"],
        ["adults_ha", "simulated_cpue", "coral_cover"])
    flow_columns = [column for column in flows if column not in
                    ("treatment", "reef_name", "year")]
    replay_flows = exact_replay(
        flows[flows.treatment == CONTROL],
        old_flows[old_flows.treatment == "idw_alpha015_cohort_full"], flow_columns)

    control_external = flows[flows.treatment == CONTROL].sort_values(["reef_name", "year"])
    external_deltas = {}
    for arm in ARMS:
        current = flows[flows.treatment == arm].sort_values(["reef_name", "year"])
        delta = float(np.max(np.abs(current.external_immigration.to_numpy() -
                                    control_external.external_immigration.to_numpy())))
        external_deltas[arm] = delta
        if delta > 1e-10:
            raise ValueError(f"External immigration changed in {arm}: {delta}")

    # Mechanism diagnostics must be inert when their switch is off.
    valid = mechanisms[mechanisms.year >= 1992]
    for arm in ARMS:
        rows = valid[valid.treatment == arm]
        if arm in STARVATION_ARMS:
            if not np.isfinite(rows.preferred_food_memory).all():
                raise ValueError(f"Missing preferred-food memory in {arm}")
        elif not rows.preferred_food_memory.isna().all():
            raise ValueError(f"Preferred-food memory active with switch off in {arm}")
        senescent = rows[["senescent_adults_ha", "senescent_deaths_ha"]].to_numpy()
        if arm in SENESCENCE_ARMS:
            if senescent[:, 1].max() <= 0:
                raise ValueError(f"Senescence never removed adults in {arm}")
        elif np.abs(senescent).max() != 0:
            raise ValueError(f"Senescence diagnostics nonzero with switch off in {arm}")
    state_columns = ["recruits_ha", "subadults_ha", "adults_ha", "maturation_ha",
                     "fast_coral", "total_coral", "food_survival",
                     "senescent_adults_ha", "senescent_deaths_ha"]
    state = valid[state_columns].to_numpy()
    if not np.isfinite(state).all() or (state < -1e-12).any():
        raise ValueError("Nonfinite or negative state/diagnostic")
    if (valid.food_survival > 1 + 1e-12).any():
        raise ValueError("Food survival above one")
    if (valid.fast_coral > valid.total_coral + 1e-12).any():
        raise ValueError("Fast coral exceeds total coral")
    merged = valid.merge(trajectories, on=["treatment", "reef_name", "year"],
                         suffixes=("", "_traj"), validate="one_to_one")
    adult_delta = float(np.max(np.abs(merged.adults_ha - merged.adults_ha_traj)))
    if adult_delta > 1e-10:
        raise ValueError(f"Mechanism adults disagree with trajectories: {adult_delta}")

    rows = []
    smoothed = scores[scores.observation_treatment == "smoothed_3y"]
    for score in smoothed.itertuples(index=False):
        series = trajectories[(trajectories.treatment == score.treatment) &
                              (trajectories.reef_name == score.reef_name) &
                              (trajectories.year >= 1992)].sort_values("year")
        mech = mechanisms[(mechanisms.treatment == score.treatment) &
                          (mechanisms.reef_name == score.reef_name) &
                          mechanisms.year.between(*GAP)]
        reef_flows = flows[(flows.treatment == score.treatment) &
                           (flows.reef_name == score.reef_name) &
                           flows.year.between(*GAP)]
        years = [int(value) for value in str(score.sim_peak_years).split(";")
                 if value and value != "nan"]
        heights = [float(value) for value in str(score.sim_peak_heights_cpue).split(";")
                   if value and value != "nan"]
        if len(years) != len(heights) or len(series) != N_YEARS - 1:
            raise ValueError("Peak data or time series incomplete")
        ratio = float("nan")
        low_years = float("nan")
        if len(years) >= 2:
            smooth_cpue = series.simulated_cpue.rolling(3, center=True, min_periods=1).mean()
            between = series.year.between(years[0], years[1])
            ratio = float(smooth_cpue[between].min() / min(heights[:2]))
            # Trough duration: years between the first two peaks below 25% of the smaller.
            low_years = int((smooth_cpue[between] < .25 * min(heights[:2])).sum())
        rows.append({
            "treatment": score.treatment, "reef_name": score.reef_name,
            "first_peak_year": years[0] if years else float("nan"),
            "second_peak_year": years[1] if len(years) > 1 else float("nan"),
            "first_peak_cpue": heights[0] if heights else float("nan"),
            "second_peak_cpue": heights[1] if len(heights) > 1 else float("nan"),
            "trough_to_smaller_peak": ratio,
            "trough_years_below_25pct": low_years,
            "matched_peaks": int(score.matched_peaks),
            "smoothed_loss": float(score.loss),
            "mean_gap_adults_ha": mech.adults_ha.mean(),
            "mean_gap_maturation_ha_per_year": mech.maturation_ha.mean(),
            "mean_gap_food_survival": mech.food_survival.mean(),
            "mean_gap_fast_coral": mech.fast_coral.mean(),
            "mean_gap_senescent_deaths_ha_per_year": mech.senescent_deaths_ha.mean(),
            "mean_gap_internal_immigration_ha_per_year": reef_flows.internal_immigration.mean(),
            "mean_gap_external_immigration_ha_per_year": reef_flows.external_immigration.mean(),
            "model_coral_2012": float(series[series.year == 2012].coral_cover.iloc[0]),
        })
    result = pd.DataFrame(rows)
    control = result[result.treatment == CONTROL].set_index("reef_name")
    for column in ("first_peak_cpue", "second_peak_cpue", "model_coral_2012"):
        result[f"{column}_relative_to_control"] = [
            row[column] / control.loc[row.reef_name, column]
            if pd.notna(row[column]) and control.loc[row.reef_name, column] > 0
            else float("nan") for _, row in result.iterrows()]
    result.to_csv(paths["structural_crash_audit.csv"], index=False)

    two_wave = result[result.reef_name.isin(TWO_WAVE)]
    summary = two_wave.groupby("treatment").agg(
        two_wave_reefs_with_two_peaks=("second_peak_year", "count"),
        worst_available_trough_ratio=("trough_to_smaller_peak", "max"),
        minimum_first_peak_vs_control=("first_peak_cpue_relative_to_control", "min"),
        minimum_second_peak_vs_control=("second_peak_cpue_relative_to_control", "min"),
        minimum_coral_2012_vs_control=("model_coral_2012_relative_to_control", "min"),
    )
    summary["matched_observed_peaks_of_seven"] = result.groupby("treatment").matched_peaks.sum()
    lizard = result[result.reef_name == "Lizard Island Reef"].set_index("treatment")
    summary["lizard_first_peak_year"] = lizard.first_peak_year
    summary["passes_qualitative_gate"] = (
        (summary.two_wave_reefs_with_two_peaks == len(TWO_WAVE)) &
        (summary.matched_observed_peaks_of_seven == 7) &
        (summary.worst_available_trough_ratio <= .10))

    manta = pd.read_csv(ROOT / "sandbox/data/reef_manta.csv")
    observations = pd.read_csv(ROOT / "sandbox/data/reef_cots.csv")
    fig, axes = plt.subplots(4, 2, figsize=(14, 13), sharex=True)
    for index, (reef, observed) in enumerate(REEFS.items()):
        cots_ax, coral_ax = axes[index]
        for arm, (label, color, style) in ARMS.items():
            values = trajectories[(trajectories.treatment == arm) &
                                  (trajectories.reef_name == reef)].sort_values("year")
            fast = mechanisms[(mechanisms.treatment == arm) &
                              (mechanisms.reef_name == reef)].sort_values("year")
            cots_ax.plot(values.year.iloc[1:], values.simulated_cpue.iloc[1:],
                         color=color, linestyle=style, linewidth=1.4,
                         label=label if index == 0 else None)
            coral_ax.plot(values.year, values.coral_cover, color=color,
                          linestyle=style, linewidth=1.4)
            coral_ax.plot(fast.year, fast.fast_coral, color=color,
                          linestyle=":", linewidth=1.0)
        obs_cots = observations[(observations.reef_name == observed) &
                                observations.year.between(1992, 2024)]
        cots_ax.scatter(obs_cots.year, obs_cots.cotsptow, s=16,
                        color="#111827", alpha=.7,
                        label="LTMP COTS/tow" if index == 0 else None)
        obs_coral = manta[(manta.reef_name == observed) &
                          (manta.data_type == "manta") &
                          (manta.domain_category == "reef") &
                          (manta.variable == "HC") &
                          (manta.purpose == "MANTA") &
                          (manta.project_code == "LTMP") &
                          (manta.depth == 9) &
                          manta.report_year.between(1992, 2024)]
        coral_ax.scatter(obs_coral.report_year, obs_coral["median"],
                         s=16, color="#111827", alpha=.6)
        cots_ax.set_ylabel(f"{reef}\nCOTS/tow")
        coral_ax.set_ylabel("Coral fraction (solid total, dotted fast)")
        cots_ax.axhline(.22, color="#92400e", linewidth=.8, linestyle=":")
        for axis in (cots_ax, coral_ax):
            axis.set_xlim(1991, 2024)
            axis.set_ylim(bottom=0)
            axis.grid(alpha=.15)
    axes[0, 0].set_title("COTS / frozen reef observation scale")
    axes[0, 1].set_title("Model coral (habitable fraction) / LTMP 9 m manta")
    axes[-1, 0].set_xlabel("Year")
    axes[-1, 1].set_xlabel("Year")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3,
               fontsize=8, bbox_to_anchor=(.5, -.01))
    fig.tight_layout(rect=(0, .05, 1, 1))
    fig.savefig(paths["cots_coral_structural_crash.png"], dpi=170, bbox_inches="tight")
    plt.close(fig)

    report = {
        "status": "bounded_diagnostic_not_promoted",
        "control_trajectory_replay_max_abs_delta": replay_trajectories,
        "control_flux_replay_max_abs_delta": replay_flows,
        "external_immigration_max_abs_delta": external_deltas,
        "mechanism_adults_vs_trajectory_max_abs_delta": adult_delta,
        "treatments": json.loads(summary.to_json(orient="index")),
        "criteria": "2 peaks on each of Lizard/MacGillivray/North Direction; 7/7 matched; "
                    "worst trough ratio <=0.10; peak and coral guardrails reported, not scored",
        "inputs": {name: sha256(run / name) for name in (
            "metadata.toml", "scores.csv", "trajectories.csv", "flows.csv",
            "mechanisms.csv")},
        "script_sha256": sha256(Path(__file__)),
        "outputs": {name: sha256(path) for name, path in paths.items()
                    if name != "structural_crash_audit.json"},
    }
    paths["structural_crash_audit.json"].write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n")
    pd.set_option("display.width", 200)
    print(summary.to_string())
    print("Control replay max delta:", max(replay_trajectories.values()),
          max(replay_flows.values()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    main(parser.parse_args().run)
