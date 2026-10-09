"""Audit the five-arm 1991 lagged adult-hazard screen and plot COTS/coral."""

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
    "density_h075": ("Density; hazard 0.75", "#0d9488", "-"),
    "density_h150": ("Density; hazard 1.50", "#0f766e", "--"),
    "density_food_h075": ("Density × food; hazard 0.75", "#dc2626", "-"),
    "density_food_h150": ("Density × food; hazard 1.50", "#7c3aed", "--"),
}
REEFS = {
    "Lizard Island Reef": "Lizard Isles",
    "MacGillivray Reef": "Macgillivray Reef",
    "North Direction Reef": "North Direction Island",
    "Eyrie Reef": "Eyrie Reef",
}
TWO_WAVE = set(list(REEFS)[:3])


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
        "lagged_hazard_audit.csv", "lagged_hazard_audit.json",
        "site_hazard_summary.csv", "cots_coral_lagged_hazard.png")}
    if any(path.exists() for path in paths.values()):
        raise FileExistsError("Refusing to overwrite lagged-hazard audit")
    trajectories = pd.read_csv(run / "trajectories.csv")
    flows = pd.read_csv(run / "flows.csv")
    hazards = pd.read_csv(run / "hazards.csv")
    scores = pd.read_csv(run / "scores.csv")
    failures = pd.read_csv(run / "failures.csv")
    site = pd.read_csv(run / "site_hazard.csv")
    if len(failures):
        raise ValueError("Failed evaluations in run")
    if set(trajectories.treatment) != set(ARMS):
        raise ValueError("Missing or extra treatment")
    if any(len(table) != 5 * 4 * 34 for table in (trajectories, flows, hazards)):
        raise ValueError("Incomplete reef-year outputs")
    if len(site) != 5 * 2895 * 34:
        raise ValueError("Incomplete site-year hazard output")
    if site.site_index.min() != 1 or site.site_index.max() != 2895 or set(site.year) != set(range(1991, 2025)):
        raise ValueError("Site/year identity mismatch")
    if site.duplicated(["treatment", "site_index", "year"]).any():
        raise ValueError("Duplicate site-year output")
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
    if site[site.treatment == CONTROL][["burden_ha", "hazard", "adult_deaths_ha"]].abs().to_numpy().max() != 0:
        raise ValueError("Off-by-default control has nonzero hazard diagnostics")
    if (site[["adults_ha", "burden_ha", "hazard", "adult_deaths_ha"]].to_numpy() < 0).any():
        raise ValueError("Negative density/burden/hazard/deaths")
    if not np.isfinite(site[["adults_ha", "burden_ha", "hazard", "adult_deaths_ha"]].to_numpy()).all():
        raise ValueError("Nonfinite site output")
    seed_adults = pd.read_csv(run / "initial_state_alpha015_cohort_full.csv").adult_ha.to_numpy()
    if len(seed_adults) != 2895:
        raise ValueError("Initial adult seed count differs from site output")
    site = site.sort_values(["treatment", "site_index", "year"])
    pre_adults = site.groupby(["treatment", "site_index"]).adults_ha.shift(1)
    is_first_transition = site.year == 1992
    pre_adults.loc[is_first_transition] = seed_adults[site.loc[is_first_transition,
                                                               "site_index"].to_numpy() - 1]
    if (site.loc[site.year >= 1992, "adult_deaths_ha"].to_numpy() >
            pre_adults.loc[site.year >= 1992].to_numpy() + 1e-10).any():
        raise ValueError("Extra adult deaths exceed pre-transition adults")
    if site.loc[site.year == 1992, "burden_ha"].abs().max() != 0:
        raise ValueError("Lagged burden must initialize at zero")
    previous_burden = site.groupby(["treatment", "site_index"]).burden_ha.shift(1)
    previous_pre_adults = pre_adults.groupby([site.treatment, site.site_index]).shift(1)
    treated_later = (site.treatment != CONTROL) & (site.year >= 1993)
    burden_error = np.max(np.abs(site.loc[treated_later, "burden_ha"].to_numpy() -
                                 (.5 * previous_burden.loc[treated_later] +
                                  .5 * previous_pre_adults.loc[treated_later]).to_numpy()))
    if burden_error > 1e-10:
        raise ValueError(f"Lagged burden recurrence failed: {burden_error}")
    if site[site.treatment != CONTROL].adult_deaths_ha.max() <= 0:
        raise ValueError("Hazard did not engage in any treatment")
    control_external = flows[flows.treatment == CONTROL].sort_values(["reef_name", "year"])
    external_deltas = {}
    for arm in ARMS:
        current = flows[flows.treatment == arm].sort_values(["reef_name", "year"])
        delta = float(np.max(np.abs(current.external_immigration.to_numpy() -
                                    control_external.external_immigration.to_numpy())))
        external_deltas[arm] = delta
        if delta > 1e-10:
            raise ValueError(f"External immigration changed in {arm}: {delta}")

    # The first model output year is an initialization placeholder, not a transition.
    valid_site = site[site.year >= 1992]
    site_summary = valid_site.groupby(["treatment", "year"], sort=True).agg(
        n_sites=("site_index", "count"),
        median_adults_ha=("adults_ha", "median"),
        p90_adults_ha=("adults_ha", lambda values: values.quantile(.90)),
        median_burden_ha=("burden_ha", "median"),
        p90_burden_ha=("burden_ha", lambda values: values.quantile(.90)),
        mean_hazard=("hazard", "mean"),
        active_hazard_site_fraction=("hazard", lambda values: (values > 0).mean()),
        mean_extra_adult_deaths_ha=("adult_deaths_ha", "mean"),
    ).reset_index()
    site_summary.to_csv(paths["site_hazard_summary.csv"], index=False)

    rows = []
    smoothed = scores[scores.observation_treatment == "smoothed_3y"]
    for score in smoothed.itertuples(index=False):
        series = trajectories[(trajectories.treatment == score.treatment) &
                              (trajectories.reef_name == score.reef_name) &
                              (trajectories.year >= 1992)].sort_values("year")
        haz = hazards[(hazards.treatment == score.treatment) &
                      (hazards.reef_name == score.reef_name) &
                      hazards.year.between(2004, 2012)]
        reef_flows = flows[(flows.treatment == score.treatment) &
                           (flows.reef_name == score.reef_name) &
                           flows.year.between(2004, 2012)]
        years = [int(value) for value in str(score.sim_peak_years).split(";")
                 if value and value != "nan"]
        heights = [float(value) for value in str(score.sim_peak_heights_cpue).split(";")
                   if value and value != "nan"]
        if len(years) != len(heights) or len(series) != 33:
            raise ValueError("Peak data or time series incomplete")
        ratio = float("nan")
        if len(years) >= 2:
            smooth_cpue = series.simulated_cpue.rolling(3, center=True, min_periods=1).mean()
            ratio = float(smooth_cpue[series.year.between(years[0], years[1])].min() /
                          min(heights[:2]))
        rows.append({
            "treatment": score.treatment, "reef_name": score.reef_name,
            "first_peak_year": years[0] if years else float("nan"),
            "second_peak_year": years[1] if len(years) > 1 else float("nan"),
            "first_peak_cpue": heights[0] if heights else float("nan"),
            "second_peak_cpue": heights[1] if len(heights) > 1 else float("nan"),
            "trough_to_smaller_peak": ratio, "matched_peaks": int(score.matched_peaks),
            "smoothed_loss": float(score.loss),
            "mean_gap_adults_ha": series[series.year.between(2004, 2012)].adults_ha.mean(),
            "model_coral_2012": float(series[series.year == 2012].coral_cover.iloc[0]),
            "mean_gap_burden_ha": haz.burden_ha.mean(),
            "mean_gap_hazard": haz.hazard.mean(),
            "mean_gap_extra_adult_deaths_ha": haz.adult_deaths_ha.mean(),
            "mean_gap_internal_immigration_ha_per_year": reef_flows.internal_immigration.mean(),
            "mean_gap_external_immigration_ha_per_year": reef_flows.external_immigration.mean(),
        })
    result = pd.DataFrame(rows)
    control = result[result.treatment == CONTROL].set_index("reef_name")
    for column in ("first_peak_cpue", "second_peak_cpue", "model_coral_2012"):
        result[f"{column}_relative_to_control"] = [
            row[column] / control.loc[row.reef_name, column]
            if pd.notna(row[column]) and control.loc[row.reef_name, column] > 0
            else float("nan") for _, row in result.iterrows()]
    result.to_csv(paths["lagged_hazard_audit.csv"], index=False)

    manta = pd.read_csv(ROOT / "sandbox/data/reef_manta.csv")
    observations = pd.read_csv(ROOT / "sandbox/data/reef_cots.csv")
    fig, axes = plt.subplots(4, 2, figsize=(14, 13), sharex=True)
    for index, (reef, observed) in enumerate(REEFS.items()):
        cots_ax, coral_ax = axes[index]
        for arm, (label, color, style) in ARMS.items():
            values = trajectories[(trajectories.treatment == arm) &
                                  (trajectories.reef_name == reef)].sort_values("year")
            cots_ax.plot(values.year.iloc[1:], values.simulated_cpue.iloc[1:],
                         color=color, linestyle=style, linewidth=1.4,
                         label=label if index == 0 else None)
            coral_ax.plot(values.year, values.coral_cover, color=color,
                          linestyle=style, linewidth=1.4)
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
        coral_ax.set_ylabel("Hard coral fraction")
        cots_ax.axhline(.22, color="#92400e", linewidth=.8, linestyle=":")
        for axis in (cots_ax, coral_ax):
            axis.set_xlim(1991, 2024)
            axis.set_ylim(bottom=0)
            axis.grid(alpha=.15)
    axes[0, 0].set_title("COTS / frozen reef observation scale")
    axes[0, 1].set_title("Whole-reef model coral / LTMP 9 m manta")
    axes[-1, 0].set_xlabel("Year")
    axes[-1, 1].set_xlabel("Year")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3,
               fontsize=8, bbox_to_anchor=(.5, -.01))
    fig.tight_layout(rect=(0, .05, 1, 1))
    fig.savefig(paths["cots_coral_lagged_hazard.png"], dpi=170, bbox_inches="tight")
    plt.close(fig)

    two_wave = result[result.reef_name.isin(TWO_WAVE)]
    summary = two_wave.groupby("treatment").agg(
        two_wave_reefs_with_two_peaks=("second_peak_year", "count"),
        worst_available_trough_ratio=("trough_to_smaller_peak", "max"),
        minimum_first_peak_vs_control=("first_peak_cpue_relative_to_control", "min"),
        minimum_second_peak_vs_control=("second_peak_cpue_relative_to_control", "min"),
        mean_gap_extra_deaths_ha=("mean_gap_extra_adult_deaths_ha", "mean"),
    )
    summary["matched_observed_peaks_of_seven"] = result.groupby("treatment").matched_peaks.sum()
    report = {
        "status": "bounded_diagnostic_not_promoted",
        "control_trajectory_replay_max_abs_delta": replay_trajectories,
        "control_flux_replay_max_abs_delta": replay_flows,
        "external_immigration_max_abs_delta": external_deltas,
        "site_burden_recurrence_max_abs_delta": float(burden_error),
        "treatments": summary.to_dict(orient="index"),
        "criteria": "2 peaks in each of three two-wave reefs; 7/7 matched; each trough ratio <=0.10; peak and coral guardrails",
        "inputs": {name: sha256(run / name) for name in (
            "metadata.toml", "scores.csv", "trajectories.csv", "flows.csv",
            "hazards.csv", "site_hazard.csv")},
        "script_sha256": sha256(Path(__file__)),
        "outputs": {name: sha256(path) for name, path in paths.items()
                    if name != "lagged_hazard_audit.json"},
    }
    paths["lagged_hazard_audit.json"].write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(summary.to_string())
    print("Control replay max delta:", max(replay_trajectories.values()),
          max(replay_flows.values()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    main(parser.parse_args().run)
