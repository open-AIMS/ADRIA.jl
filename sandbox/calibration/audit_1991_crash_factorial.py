"""Independently audit the frozen 1991 cohort crash factorial and plot COTS/coral."""

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
    CONTROL: ("Full cohort control", "#475569", "-"),
    "cohort_full_stage_gate": ("Stage gate", "#0d9488", "-"),
    "cohort_full_mortality030": ("Adult mortality 0.30", "#dc2626", "--"),
    "cohort_full_stage_gate_mortality030": ("Stage gate + mortality", "#7c3aed", "-"),
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


def audit(run: Path) -> None:
    artifacts = [run / name for name in (
        "crash_factorial_audit.csv", "crash_factorial_audit.json",
        "cots_coral_crash_factorial.png")]
    if any(path.exists() for path in artifacts):
        raise FileExistsError("Refusing to overwrite crash-factorial audit")
    trajectories = pd.read_csv(run / "trajectories.csv")
    flows = pd.read_csv(run / "flows.csv")
    scores = pd.read_csv(run / "scores.csv")
    failures = pd.read_csv(run / "failures.csv")
    if len(failures):
        raise ValueError("Failed evaluations in run")
    if set(trajectories.treatment) != set(ARMS):
        raise ValueError("Missing or extra treatment")
    if len(trajectories) != 4 * 4 * 34 or len(flows) != len(trajectories):
        raise ValueError("Incomplete treatment trajectories or flows")
    if sha256(run / "initial_state_alpha015_cohort_full.csv") != sha256(
        ARCHIVE / "initial_state_alpha015_cohort_full.csv"):
        raise ValueError("Initial cohort state differs from archived anchor")

    reference = pd.read_csv(ARCHIVE / "trajectories.csv")
    reference = reference[reference.treatment == "idw_alpha015_cohort_full"]
    control = trajectories[trajectories.treatment == CONTROL]
    keys = ["reef_name", "year"]
    joined = control.merge(reference, on=keys, suffixes=("_new", "_old"), validate="one_to_one")
    if len(joined) != len(control):
        raise ValueError("Control/reference year rows differ")
    replay = {}
    for column in ("adults_ha", "simulated_cpue", "coral_cover"):
        replay[column] = float(np.max(np.abs(
            joined[f"{column}_new"] - joined[f"{column}_old"])))
    if max(replay.values()) > 1e-10:
        raise ValueError(f"Archived cohort control did not replay: {replay}")
    old_flows = pd.read_csv(ARCHIVE / "flows.csv")
    old_flows = old_flows[old_flows.treatment == "idw_alpha015_cohort_full"]
    flow_join = flows[flows.treatment == CONTROL].merge(
        old_flows, on=keys, suffixes=("_new", "_old"), validate="one_to_one")
    flow_cols = [column for column in flows if column not in ("treatment", *keys)]
    flow_replay = max(float(np.max(np.abs(flow_join[f"{name}_new"] -
                                          flow_join[f"{name}_old"])))
                      for name in flow_cols)
    if flow_replay > 1e-10:
        raise ValueError(f"Archived cohort fluxes did not replay: {flow_replay}")
    expected_external = flows[flows.treatment == CONTROL].sort_values(keys)
    external_max_diff = {}
    for treatment in ARMS:
        current = flows[flows.treatment == treatment].sort_values(keys)
        delta = np.max(np.abs(current.external_immigration.to_numpy() -
                              expected_external.external_immigration.to_numpy()))
        external_max_diff[treatment] = float(delta)
        if delta > 1e-10:
            raise ValueError(f"External supply changed in {treatment}: {delta}")

    rows = []
    smoothed = scores[scores.observation_treatment == "smoothed_3y"]
    for score in smoothed.itertuples(index=False):
        series = trajectories[(trajectories.treatment == score.treatment) &
                              (trajectories.reef_name == score.reef_name) &
                              (trajectories.year >= 1992)].sort_values("year")
        reef_flows = flows[(flows.treatment == score.treatment) &
                           (flows.reef_name == score.reef_name)].sort_values("year")
        if len(series) != 33 or len(reef_flows) != 34:
            raise ValueError("Incomplete reef time series")
        years = [int(value) for value in str(score.sim_peak_years).split(";")
                 if value and value != "nan"]
        heights = [float(value) for value in str(score.sim_peak_heights_cpue).split(";")
                   if value and value != "nan"]
        if len(years) != len(heights):
            raise ValueError("Peak years and heights disagree")
        smooth_cpue = series.simulated_cpue.rolling(3, center=True, min_periods=1).mean()
        ratio = float("nan")
        gap = series[series.year.between(2004, 2012)]
        if len(years) >= 2:
            between = series.year.between(years[0], years[1])
            ratio = float(smooth_cpue[between].min() / min(heights[:2]))
        flow_gap = reef_flows[reef_flows.year.between(2004, 2012)]
        rows.append({
            "treatment": score.treatment, "reef_name": score.reef_name,
            "first_peak_year": years[0] if years else float("nan"),
            "second_peak_year": years[1] if len(years) > 1 else float("nan"),
            "first_peak_cpue": heights[0] if heights else float("nan"),
            "second_peak_cpue": heights[1] if len(heights) > 1 else float("nan"),
            "trough_to_smaller_peak": ratio, "matched_peaks": int(score.matched_peaks),
            "smoothed_loss": float(score.loss),
            "mean_adults_ha_2004_2012": gap.adults_ha.mean(),
            "mean_coral_2004_2012": gap.coral_cover.mean(),
            "model_coral_2012": float(series[series.year == 2012].coral_cover.iloc[0]),
            "minimum_model_coral_1992_2024": series.coral_cover.min(),
            "mean_maturation_2004_2012": flow_gap.maturation.mean(),
            "mean_retained_juveniles_2004_2012": flow_gap.retained_juveniles.mean(),
            "mean_local_fecundity_2004_2012": flow_gap.local_fecundity.mean(),
            "mean_internal_immigration_2004_2012": flow_gap.internal_immigration.mean(),
            "mean_external_immigration_2004_2012": flow_gap.external_immigration.mean(),
        })
    result = pd.DataFrame(rows)
    control_rows = result[result.treatment == CONTROL].set_index("reef_name")
    for column in ("first_peak_cpue", "second_peak_cpue", "model_coral_2012"):
        result[f"{column}_relative_to_control"] = [
            row[column] / control_rows.loc[row.reef_name, column]
            if pd.notna(row[column]) and control_rows.loc[row.reef_name, column] > 0
            else float("nan") for _, row in result.iterrows()]
    result.to_csv(artifacts[0], index=False)

    manta = pd.read_csv(ROOT / "sandbox/data/reef_manta.csv")
    obs = pd.read_csv(ROOT / "sandbox/data/reef_cots.csv")
    fig, axes = plt.subplots(4, 2, figsize=(14, 13), sharex=True)
    for index, (reef, observed) in enumerate(REEFS.items()):
        cots_ax, coral_ax = axes[index]
        for treatment, (label, color, style) in ARMS.items():
            series = trajectories[(trajectories.treatment == treatment) &
                                  (trajectories.reef_name == reef)].sort_values("year")
            cots_ax.plot(series.year.iloc[1:], series.simulated_cpue.iloc[1:],
                         color=color, linestyle=style, linewidth=1.5,
                         label=label if index == 0 else None)
            coral_ax.plot(series.year, series.coral_cover, color=color,
                          linestyle=style, linewidth=1.5)
        observed_cots = obs[(obs.reef_name == observed) & obs.year.between(1992, 2024)]
        cots_ax.scatter(observed_cots.year, observed_cots.cotsptow,
                        color="#111827", alpha=0.7, s=16,
                        label="LTMP COTS/tow" if index == 0 else None)
        coral_obs = manta[(manta.reef_name == observed) &
                          (manta.data_type == "manta") &
                          (manta.domain_category == "reef") &
                          (manta.variable == "HC") &
                          (manta.purpose == "MANTA") &
                          (manta.project_code == "LTMP") &
                          (manta.depth == 9) &
                          manta.report_year.between(1992, 2024)]
        coral_ax.scatter(coral_obs.report_year, coral_obs["median"],
                         color="#111827", alpha=0.6, s=16)
        cots_ax.set_ylabel(f"{reef}\nCOTS/tow")
        coral_ax.set_ylabel("Hard coral fraction")
        cots_ax.axhline(0.22, color="#92400e", linewidth=0.8, linestyle=":")
        for axis in (cots_ax, coral_ax):
            axis.set_xlim(1991, 2024)
            axis.set_ylim(bottom=0)
            axis.grid(alpha=0.15)
    axes[0, 0].set_title("COTS / fixed reef observation scale")
    axes[0, 1].set_title("Whole-reef model coral / LTMP 9 m manta")
    axes[-1, 0].set_xlabel("Year")
    axes[-1, 1].set_xlabel("Year")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3,
               fontsize=8, bbox_to_anchor=(0.5, -0.01))
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(artifacts[2], dpi=170, bbox_inches="tight")
    plt.close(fig)

    two_wave = result[result.reef_name.isin(TWO_WAVE)]
    summary = two_wave.groupby("treatment").agg(
        max_trough_ratio=("trough_to_smaller_peak", "max"),
        two_wave_matched_peaks=("matched_peaks", "sum"),
        min_first_peak_vs_control=("first_peak_cpue_relative_to_control", "min"),
        min_second_peak_vs_control=("second_peak_cpue_relative_to_control", "min"),
        min_coral_2012_vs_control=("model_coral_2012_relative_to_control", "min"),
    )
    summary["all_reef_matched_peaks"] = result.groupby("treatment").matched_peaks.sum()
    summary["two_wave_reefs_with_two_peaks"] = two_wave.groupby("treatment").second_peak_year.count()
    report = {
        "status": "diagnostic_only_no_promotion",
        "replay_max_abs_delta": replay,
        "flow_replay_max_abs_delta": flow_replay,
        "external_immigration_max_abs_delta": external_max_diff,
        "screen_criterion": "all 7 observed peaks matched and each two-wave trough ratio <=0.10 without substantial peak/coral degradation",
        "treatments": summary.to_dict(orient="index"),
        "inputs": {name: sha256(run / name) for name in
                   ("trajectories.csv", "flows.csv", "scores.csv", "metadata.toml")},
        "script_sha256": sha256(Path(__file__)),
        "outputs": {path.name: sha256(path) for path in (artifacts[0], artifacts[2])},
    }
    artifacts[1].write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(summary.to_string())
    print("Control replay max absolute delta:", max(replay.values()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    audit(parser.parse_args().run)
