"""Audit the 1991 q=0.015 unit migration and consumption factorial.

Checks that the compatibility arm reproduces the archived control times k, that
external supply and switch diagnostics behave as specified, then reports the
preregistered peak/trough gate in absolute CPUE (q=0.015 on every reef).
See crash_mechanism_protocol.md (2026-10-09 observation-unit migration).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import tomllib
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
ARCHIVE = ROOT / "sandbox/calibration/runs/20261007T199111_owen_1991_cohort_history_final"
COMPAT = "compat_rescaled_exact"
CONTROL = "migrated_control"
ARMS = {
    COMPAT: ("Compatibility (archived ×k)", "#94a3b8", ":"),
    CONTROL: ("Migrated control", "#475569", "-"),
    "migrated_structural": ("Migrated + starvation + senescence", "#0d9488", "-"),
    "per_capita_legacy": ("Per-capita grazing", "#dc2626", "-"),
    "per_capita_structural": ("Per-capita + starvation + senescence", "#7c3aed", "--"),
}
STRUCTURAL_ARMS = {"migrated_structural", "per_capita_structural"}
REEFS = {
    "Lizard Island Reef": "Lizard Isles",
    "MacGillivray Reef": "Macgillivray Reef",
    "North Direction Reef": "North Direction Island",
    "Eyrie Reef": "Eyrie Reef",
}
TWO_WAVE = set(list(REEFS)[:3])
GAP = (2004, 2012)
N_YEARS = 34
Q = 0.015
SCALE_INVARIANT_FLOWS = {"settlement_gate"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def scaled_replay(current: pd.DataFrame, old: pd.DataFrame, columns: list[str],
                  factors: dict[str, float]) -> dict[str, float]:
    """Maximum relative difference between current/factor and the archive."""
    joined = current.merge(old, on=["reef_name", "year"], suffixes=("_new", "_old"),
                           validate="one_to_one")
    if len(joined) != len(current) or len(joined) != len(old):
        raise ValueError("Archived control does not cover the same reef-years")
    deltas = {}
    for column in columns:
        new = joined[f"{column}_new"].to_numpy() / factors.get(column, 1.0)
        ref = joined[f"{column}_old"].to_numpy()
        deltas[column] = float(np.max(np.abs(new - ref) / np.maximum(1.0, np.abs(ref))))
    if max(deltas.values()) > 1e-9:
        raise ValueError(f"Compatibility arm did not replay archive x k: {deltas}")
    return deltas


def main(run: Path) -> None:
    paths = {name: run / name for name in (
        "unit_migration_audit.csv", "unit_migration_audit.json",
        "cots_coral_unit_migration.png")}
    if any(path.exists() for path in paths.values()):
        raise FileExistsError("Refusing to overwrite unit-migration audit")
    metadata = tomllib.loads((run / "metadata.toml").read_text(encoding="utf-8"))
    k = float(metadata["unit_migration"]["k_model_unit_to_cots_ha"])
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
        raise ValueError("Real-unit initial state differs from archived cohort anchor")
    scaled_initial = pd.read_csv(run / "initial_state_alpha015_cohort_full_times_k.csv")
    archived_initial = pd.read_csv(ARCHIVE / "initial_state_alpha015_cohort_full.csv")
    initial_delta = float(np.max(np.abs(
        scaled_initial[["juvenile_ha", "subadult_ha", "adult_ha"]].to_numpy() / k -
        archived_initial[["juvenile_ha", "subadult_ha", "adult_ha"]].to_numpy())))
    if initial_delta > 1e-9:
        raise ValueError("Compatibility initial state is not archived x k")

    old_traj = pd.read_csv(ARCHIVE / "trajectories.csv")
    old_traj = old_traj[old_traj.treatment == "idw_alpha015_cohort_full"]
    old_flows = pd.read_csv(ARCHIVE / "flows.csv")
    old_flows = old_flows[old_flows.treatment == "idw_alpha015_cohort_full"]
    compat_traj = trajectories[trajectories.treatment == COMPAT]
    replay_traj = scaled_replay(compat_traj, old_traj, ["adults_ha", "coral_cover"],
                                {"adults_ha": k})
    flow_columns = [column for column in flows if column not in
                    ("treatment", "reef_name", "year")]
    replay_flows = scaled_replay(
        flows[flows.treatment == COMPAT], old_flows, flow_columns,
        {column: k for column in flow_columns if column not in SCALE_INVARIANT_FLOWS})
    lizard_new = compat_traj[compat_traj.reef_name == "Lizard Island Reef"].sort_values("year")
    lizard_old = old_traj[old_traj.reef_name == "Lizard Island Reef"].sort_values("year")
    lizard_cpue_delta = float(np.max(np.abs(lizard_new.simulated_cpue.to_numpy() -
                                            lizard_old.simulated_cpue.to_numpy())))
    if lizard_cpue_delta > 1e-9:
        raise ValueError(f"Lizard CPUE changed under compatible migration: {lizard_cpue_delta}")

    reference = flows[flows.treatment == CONTROL].sort_values(["reef_name", "year"])
    external_deltas = {}
    for arm in ARMS:
        current = flows[flows.treatment == arm].sort_values(["reef_name", "year"])
        delta = float(np.max(np.abs(current.external_immigration.to_numpy() -
                                    reference.external_immigration.to_numpy())))
        external_deltas[arm] = delta
        if delta > 1e-9:
            raise ValueError(f"External immigration changed in {arm}: {delta}")

    valid = mechanisms[mechanisms.year >= 1992]
    for arm in ARMS:
        rows = valid[valid.treatment == arm]
        if arm in STRUCTURAL_ARMS:
            if not np.isfinite(rows.preferred_food_memory).all():
                raise ValueError(f"Missing preferred-food memory in {arm}")
            if rows.senescent_deaths_ha.max() <= 0:
                raise ValueError(f"Senescence never removed adults in {arm}")
        else:
            if not rows.preferred_food_memory.isna().all():
                raise ValueError(f"Preferred-food memory active with switch off in {arm}")
            if rows[["senescent_adults_ha", "senescent_deaths_ha"]].abs().to_numpy().max() != 0:
                raise ValueError(f"Senescence diagnostics nonzero with switch off in {arm}")
    state_columns = ["recruits_ha", "subadults_ha", "adults_ha", "maturation_ha",
                     "fast_coral", "total_coral", "food_survival",
                     "senescent_adults_ha", "senescent_deaths_ha"]
    state = valid[state_columns].to_numpy()
    if not np.isfinite(state).all() or (state < -1e-12).any():
        raise ValueError("Nonfinite or negative state/diagnostic")
    scale_check = trajectories[trajectories.year >= 1992]
    scale_delta = float(np.max(np.abs(scale_check.simulated_cpue - Q * scale_check.adults_ha)))
    if scale_delta > 1e-12:
        raise ValueError("Simulated CPUE is not q times adults")

    observations = pd.read_csv(ROOT / "sandbox/data/reef_cots.csv")
    rows = []
    smoothed = scores[scores.observation_treatment == "smoothed_3y"]
    for score in smoothed.itertuples(index=False):
        series = trajectories[(trajectories.treatment == score.treatment) &
                              (trajectories.reef_name == score.reef_name) &
                              (trajectories.year >= 1992)].sort_values("year")
        mech = mechanisms[(mechanisms.treatment == score.treatment) &
                          (mechanisms.reef_name == score.reef_name) &
                          mechanisms.year.between(*GAP)]
        obs = observations[(observations.reef_name == REEFS[score.reef_name]) &
                           observations.year.between(1992, 2024)]
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
            low_years = int((smooth_cpue[between] < .25 * min(heights[:2])).sum())
        rows.append({
            "treatment": score.treatment, "reef_name": score.reef_name,
            "first_peak_year": years[0] if years else float("nan"),
            "second_peak_year": years[1] if len(years) > 1 else float("nan"),
            "first_peak_cpue": heights[0] if heights else float("nan"),
            "second_peak_cpue": heights[1] if len(heights) > 1 else float("nan"),
            "observed_max_cpue": float(obs.cotsptow.max()),
            "trough_to_smaller_peak": ratio,
            "trough_years_below_25pct": low_years,
            "matched_peaks": int(score.matched_peaks),
            "smoothed_loss": float(score.loss),
            "max_adults_ha": float(series.adults_ha.max()),
            "mean_gap_adults_ha": mech.adults_ha.mean(),
            "mean_gap_food_survival": mech.food_survival.mean(),
            "mean_gap_fast_coral": mech.fast_coral.mean(),
            "mean_gap_total_coral": mech.total_coral.mean(),
            "mean_gap_senescent_deaths_ha_per_year": mech.senescent_deaths_ha.mean(),
            "model_coral_2012": float(series[series.year == 2012].coral_cover.iloc[0]),
        })
    result = pd.DataFrame(rows)
    result.to_csv(paths["unit_migration_audit.csv"], index=False)

    two_wave = result[result.reef_name.isin(TWO_WAVE)]
    summary = two_wave.groupby("treatment").agg(
        two_wave_reefs_with_two_peaks=("second_peak_year", "count"),
        worst_available_trough_ratio=("trough_to_smaller_peak", "max"),
        minimum_coral_2012=("model_coral_2012", "min"),
    )
    summary["matched_observed_peaks_of_seven"] = result.groupby("treatment").matched_peaks.sum()
    lizard = result[result.reef_name == "Lizard Island Reef"].set_index("treatment")
    summary["lizard_first_peak_year"] = lizard.first_peak_year
    summary["lizard_first_peak_cpue"] = lizard.first_peak_cpue
    summary["lizard_max_adults_ha"] = lizard.max_adults_ha
    summary["passes_qualitative_gate"] = (
        (summary.two_wave_reefs_with_two_peaks == len(TWO_WAVE)) &
        (summary.matched_observed_peaks_of_seven == 7) &
        (summary.worst_available_trough_ratio <= .10))

    manta = pd.read_csv(ROOT / "sandbox/data/reef_manta.csv")
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
            if arm != COMPAT:
                coral_ax.plot(fast.year, fast.fast_coral, color=color,
                              linestyle=":", linewidth=.9)
        obs_cots = observations[(observations.reef_name == observed) &
                                observations.year.between(1992, 2024)]
        cots_ax.scatter(obs_cots.year, obs_cots.cotsptow, s=16, color="#111827",
                        alpha=.7, label="LTMP COTS/tow" if index == 0 else None)
        obs_coral = manta[(manta.reef_name == observed) &
                          (manta.data_type == "manta") &
                          (manta.domain_category == "reef") &
                          (manta.variable == "HC") &
                          (manta.purpose == "MANTA") &
                          (manta.project_code == "LTMP") &
                          (manta.depth == 9) &
                          manta.report_year.between(1992, 2024)]
        coral_ax.scatter(obs_coral.report_year, obs_coral["median"], s=16,
                         color="#111827", alpha=.6)
        cots_ax.set_ylabel(f"{reef}\nCOTS/tow")
        coral_ax.set_ylabel("Coral fraction (solid total, dotted fast)")
        cots_ax.axhline(.22, color="#92400e", linewidth=.8, linestyle=":")
        for axis in (cots_ax, coral_ax):
            axis.set_xlim(1991, 2024)
            axis.set_ylim(bottom=0)
            axis.grid(alpha=.15)
    axes[0, 0].set_title("COTS at fixed q = 0.015 COTS/tow per COTS/ha")
    axes[0, 1].set_title("Model coral (habitable fraction) / LTMP 9 m manta")
    axes[-1, 0].set_xlabel("Year")
    axes[-1, 1].set_xlabel("Year")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, fontsize=8,
               bbox_to_anchor=(.5, -.01))
    fig.tight_layout(rect=(0, .05, 1, 1))
    fig.savefig(paths["cots_coral_unit_migration.png"], dpi=170, bbox_inches="tight")
    plt.close(fig)

    report = {
        "status": "bounded_diagnostic_not_promoted",
        "observation_model": "fixed_q0015",
        "k_model_unit_to_cots_ha": k,
        "compat_initial_state_max_abs_delta": initial_delta,
        "compat_trajectory_replay_max_rel_delta": replay_traj,
        "compat_flux_replay_max_rel_delta": replay_flows,
        "compat_lizard_cpue_max_abs_delta": lizard_cpue_delta,
        "external_immigration_max_abs_delta_vs_migrated_control": external_deltas,
        "simulated_cpue_equals_q_times_adults_max_abs_delta": scale_delta,
        "treatments": json.loads(summary.to_json(orient="index")),
        "criteria": "2 peaks on each of Lizard/MacGillivray/North Direction; 7/7 matched; "
                    "worst trough ratio <=0.10; absolute CPUE under q=0.015",
        "inputs": {name: sha256(run / name) for name in (
            "metadata.toml", "scores.csv", "trajectories.csv", "flows.csv",
            "mechanisms.csv")},
        "script_sha256": sha256(Path(__file__)),
        "outputs": {name: sha256(path) for name, path in paths.items()
                    if name != "unit_migration_audit.json"},
    }
    paths["unit_migration_audit.json"].write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n")
    pd.set_option("display.width", 220)
    print(summary.to_string())
    print("Compat replay max rel delta:", max(replay_traj.values()), max(replay_flows.values()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    main(parser.parse_args().run)
