"""Audit paired 1991 cycle-capacity runs without changing the frozen score."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
GATE = ROOT / "sandbox/calibration/runs/20261001T113401_lizard_domain_gate/scores.csv"
TWO_WAVE = {"Lizard Island Reef", "MacGillivray Reef", "North Direction Reef"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def numbers(value: object, kind=float) -> list:
    if pd.isna(value) or str(value).strip() == "":
        return []
    return [kind(item) for item in str(value).split(";") if item]


def main(run: Path, allow_incomplete: bool = False) -> None:
    output = run / "capacity_audit.csv"
    treatment_output = run / "capacity_treatment_summary.csv"
    summary_path = run / "capacity_audit.json"
    if output.exists() or treatment_output.exists() or summary_path.exists():
        raise FileExistsError("Refusing to overwrite existing capacity audit")
    design = pd.read_csv(run / "design.csv")
    scores = pd.read_csv(run / "scores.csv")
    trajectories = pd.read_csv(run / "trajectories.csv")
    flows = pd.read_csv(run / "flows.csv")
    failures = pd.read_csv(run / "failures.csv")
    expected_treatments = {
        f"pair_{int(pair):03d}_theta_{theta}"
        for pair in design.pair for theta in (3, 1)
    }
    actual_treatments = set(scores.treatment)
    missing_treatments = sorted(expected_treatments - actual_treatments)
    if missing_treatments and not allow_incomplete:
        raise ValueError(f"Incomplete screen: {missing_treatments}")
    for table, columns in ((scores, ("loss", "sim_peak_count", "matched_peaks")),
                           (trajectories, ("adults_ha", "simulated_cpue", "coral_cover")),
                           (flows, ("local_retention", "internal_immigration",
                                    "external_immigration"))):
        if not np.isfinite(table[list(columns)].to_numpy(dtype=float)).all():
            raise ValueError(f"Nonfinite values in {columns}")
    if trajectories.duplicated(["treatment", "reef_name", "year"]).any():
        raise ValueError("Duplicate treatment/reef/year trajectory")
    gate = pd.read_csv(GATE)
    observed = gate[(gate.seed == 20260930) &
                    (gate.treatment == "v1_legacy_closed") &
                    (gate.observation_treatment == "smoothed_3y")]
    obs_peaks = {row.reef_name: (numbers(row.obs_peak_years, int),
                                 numbers(row.obs_peak_heights_cpue))
                 for row in observed.itertuples(index=False)}
    if len(obs_peaks) != 4:
        raise ValueError("Incomplete frozen observed peaks")
    rows = []
    for score in scores.itertuples(index=False):
        series = trajectories[(trajectories.treatment == score.treatment) &
                              (trajectories.reef_name == score.reef_name) &
                              (trajectories.year >= 1992)].sort_values("year")
        flux = flows[(flows.treatment == score.treatment) &
                     (flows.reef_name == score.reef_name)]
        if len(series) != 33 or len(flux) != 34:
            raise ValueError(f"Incomplete treatment: {score.treatment}/{score.reef_name}")
        years = numbers(score.sim_peak_years, int)
        heights = numbers(score.sim_peak_heights_cpue)
        if len(years) != len(heights):
            raise ValueError("Peak year/height length mismatch")
        expected_years, expected_heights = obs_peaks[score.reef_name]
        smooth_cpue = series.simulated_cpue.rolling(3, center=True,
                                                     min_periods=1).mean()
        if len(years) >= 2 and min(heights[:2]) > 0:
            between = series.year.between(years[0], years[1])
            trough = float(smooth_cpue[between].min())
            ratio = trough / min(heights[:2])
        else:
            trough = np.nan
            ratio = np.nan
        if len(years) == len(expected_years):
            max_year_error = max(abs(a - b) for a, b in zip(years, expected_years))
            height_ratios = [a / b for a, b in zip(heights, expected_heights)]
            height_in_half_to_double = all(0.5 <= value <= 2 for value in height_ratios)
        else:
            max_year_error = np.nan
            height_ratios = []
            height_in_half_to_double = False
        coral_2015 = series.loc[series.year == 2015, "coral_cover"]
        if len(coral_2015) != 1:
            raise ValueError("Missing 2015 coral")
        quiet = flux[flux.year.between(2004, 2009)]
        rows.append({
            "treatment": score.treatment,
            "pair": score.pair,
            "theta_ha": score.theta_ha,
            "reef_name": score.reef_name,
            "q_tow_per_adult_ha": float(design.loc[design.pair == score.pair,
                                                     "q_tow_per_adult_ha"].iloc[0]),
            "sim_peak_count": score.sim_peak_count,
            "expected_peak_count": len(expected_years),
            "matched_peaks": score.matched_peaks,
            "sim_peak_years": score.sim_peak_years,
            "observed_peak_years": ";".join(map(str, expected_years)),
            "sim_peak_heights_cpue": score.sim_peak_heights_cpue,
            "observed_peak_heights_cpue": ";".join(map(str, expected_heights)),
            "max_peak_year_error": max_year_error,
            "height_ratios": ";".join(f"{value:.5f}" for value in height_ratios),
            "all_heights_half_to_double": height_in_half_to_double,
            "smoothed_interpeak_min_cpue": trough,
            "trough_to_smaller_peak": ratio,
            "frozen_loss": score.loss,
            "minimum_coral_1992_2024": float(series.coral_cover.min()),
            "coral_2015": float(coral_2015.iloc[0]),
            "mean_external_settlement_2004_2009_ha_year":
                float(quiet.external_immigration.mean()),
            "mean_internal_settlement_2004_2009_ha_year":
                float(quiet.internal_immigration.mean()),
            "mean_local_retention_2004_2009_ha_year":
                float(quiet.local_retention.mean()),
        })
    result = pd.DataFrame(rows)
    result.to_csv(output, index=False)
    comparisons = []
    for treatment, group in result.groupby("treatment"):
        if len(group) != 4:
            raise ValueError(f"Incomplete score group: {treatment}")
        two = group[group.reef_name.isin(TWO_WAVE)]
        exact_peaks = bool((group.sim_peak_count == group.expected_peak_count).all())
        seven_matched = int(group.matched_peaks.sum()) == 7
        deep = bool((two.trough_to_smaller_peak <= 0.10).all())
        timing = bool((group.max_peak_year_error <= 4).all())
        amplitude = bool(group.all_heights_half_to_double.all())
        comparisons.append({
            "treatment": treatment,
            "pair": int(group.pair.iloc[0]),
            "theta_ha": float(group.theta_ha.iloc[0]),
            "exact_peak_counts": exact_peaks,
            "seven_matched_peaks": seven_matched,
            "all_three_troughs_below_0p10": deep,
            "max_trough_ratio": float(two.trough_to_smaller_peak.max())
                if two.trough_to_smaller_peak.notna().any() else None,
            "all_peak_years_within_4y": timing,
            "all_peak_heights_within_factor_2": amplitude,
            "qualitative_cycle_pass": exact_peaks and seven_matched and deep and
                timing and amplitude,
            "mean_frozen_loss": float(group.frozen_loss.mean()),
            "minimum_coral": float(group.minimum_coral_1992_2024.min()),
        })
    comparisons.sort(key=lambda item: (not item["qualitative_cycle_pass"],
                                       not item["exact_peak_counts"],
                                       not item["seven_matched_peaks"],
                                       item["max_trough_ratio"] if
                                       item["max_trough_ratio"] is not None else 10,
                                       item["mean_frozen_loss"]))
    pd.DataFrame(comparisons).to_csv(treatment_output, index=False)
    summary = {
        "status": "incomplete_diagnostic_screen_no_promotion" if missing_treatments
            else "diagnostic_screen_no_promotion",
        "missing_treatments": missing_treatments,
        "criteria": "Exact 2/2/2/1 peak counts; seven matched peaks; three trough ratios <=0.10; all peak years within 4 years; all heights within 0.5-2x observed. Coral reported separately because whole-reef versus 9m manta support differs.",
        "completed_treatments": len(comparisons),
        "failed_treatments": int(len(failures)),
        "qualitative_passes_theta3": sum(item["qualitative_cycle_pass"] and
                                           item["theta_ha"] == 3 for item in comparisons),
        "qualitative_passes_theta1_counterfactual": sum(
            item["qualitative_cycle_pass"] and item["theta_ha"] == 1
            for item in comparisons),
        "treatments": comparisons,
        "inputs_sha256": {path.name: sha256(path) for path in
                          (run / "design.csv", run / "scores.csv",
                           run / "trajectories.csv", run / "flows.csv",
                           run / "failures.csv", GATE)},
        "script_sha256": sha256(Path(__file__)),
        "output_sha256": sha256(output),
        "treatment_output_sha256": sha256(treatment_output),
    }
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(pd.DataFrame(comparisons).head(12).to_string(index=False))
    print("Qualitative passes: theta=3:", summary["qualitative_passes_theta3"],
          "theta=1 counterfactual:", summary["qualitative_passes_theta1_counterfactual"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("--allow-incomplete", action="store_true")
    args = parser.parse_args()
    main(args.run, args.allow_incomplete)
