"""Audit the bounded opt-in size-weighted 1991 screen against frozen peaks."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from audit_1991_cycle_capacity import GATE, TWO_WAVE, numbers, sha256


def main(run: Path) -> None:
    output = run / "size_audit.csv"
    summary_path = run / "size_audit.json"
    if output.exists() or summary_path.exists():
        raise FileExistsError("Refusing to overwrite existing size audit")
    design = pd.read_csv(run / "design.csv")
    scores = pd.read_csv(run / "scores.csv")
    trajectories = pd.read_csv(run / "trajectories.csv")
    flows = pd.read_csv(run / "flows.csv")
    failures = pd.read_csv(run / "failures.csv")
    expected = set(design.treatment)
    completed = set(scores.treatment)
    missing = sorted(expected - completed)
    if missing:
        raise ValueError(f"Incomplete size screen: {missing}")
    if len(failures):
        raise ValueError(f"Failed evaluations: {failures.treatment.tolist()}")
    if trajectories.duplicated(["treatment", "reef_name", "year"]).any():
        raise ValueError("Duplicate treatment/reef/year trajectory")
    if flows.duplicated(["treatment", "reef_name", "year"]).any():
        raise ValueError("Duplicate treatment/reef/year flow")
    for table, fields in ((scores, ("loss", "sim_peak_count", "matched_peaks")),
                          (trajectories, ("adults_ha", "simulated_cpue", "coral_cover")),
                          (flows, ("small_adults", "large_adults", "effective_breeders",
                                   "small_adult_growth", "local_fecundity",
                                   "internal_immigration", "external_immigration"))):
        if not np.isfinite(table[list(fields)].to_numpy(dtype=float)).all():
            raise ValueError(f"Nonfinite {fields}")
    if ((flows.small_adults < -1e-10) | (flows.large_adults < -1e-10) |
            (flows.effective_breeders < -1e-10)).any():
        raise ValueError("Negative size density or effective breeders")
    merged = trajectories.merge(flows, on=["treatment", "reef_name", "year"],
                                validate="one_to_one")
    if len(merged) != len(trajectories):
        raise ValueError("Missing size diagnostics")
    if not np.allclose(merged.small_adults + merged.large_adults,
                       merged.adults_ha, rtol=1e-8, atol=1e-8):
        raise ValueError("Size density does not sum to adult density")

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
        all_series = merged[(merged.treatment == score.treatment) &
                            (merged.reef_name == score.reef_name)].sort_values("year")
        series = all_series[all_series.year >= 1992]
        if len(all_series) != 34 or len(series) != 33:
            raise ValueError(f"Incomplete {score.treatment}/{score.reef_name}")
        years = numbers(score.sim_peak_years, int)
        heights = numbers(score.sim_peak_heights_cpue)
        expected_years, expected_heights = obs_peaks[score.reef_name]
        smooth_cpue = series.simulated_cpue.rolling(3, center=True,
                                                     min_periods=1).mean()
        trough = np.nan
        ratio = np.nan
        if len(years) >= 2 and min(heights[:2]) > 0:
            between = series.year.between(years[0], years[1])
            trough = float(smooth_cpue[between].min())
            ratio = trough / min(heights[:2])
        timing = np.nan
        amplitude = False
        if len(years) == len(expected_years):
            timing = max(abs(a - b) for a, b in zip(years, expected_years))
            amplitude = all(0.5 <= a / b <= 2.0 for a, b in
                            zip(heights, expected_heights))
        gap = all_series[all_series.year.between(2004, 2009)]
        rows.append({
            "treatment": score.treatment,
            "reef_name": score.reef_name,
            "sim_peak_count": score.sim_peak_count,
            "expected_peak_count": len(expected_years),
            "matched_peaks": score.matched_peaks,
            "sim_peak_years": score.sim_peak_years,
            "observed_peak_years": ";".join(map(str, expected_years)),
            "trough_cpue": trough,
            "trough_to_smaller_peak": ratio,
            "max_peak_year_error": timing,
            "all_heights_half_to_double": amplitude,
            "frozen_loss": score.loss,
            "minimum_coral_1992_2024": float(series.coral_cover.min()),
            "coral_2015": float(series.loc[series.year == 2015, "coral_cover"].iloc[0]),
            "gap_effective_breeders_ha": float(gap.effective_breeders.mean()),
            "gap_large_adults_ha": float(gap.large_adults.mean()),
            "gap_local_production_ha_year": float(gap.local_fecundity.mean()),
            "gap_internal_settlement_ha_year": float(gap.internal_immigration.mean()),
            "gap_external_settlement_ha_year": float(gap.external_immigration.mean()),
        })
    result = pd.DataFrame(rows)
    result.to_csv(output, index=False)
    comparison = []
    for treatment, group in result.groupby("treatment", sort=False):
        if len(group) != 4:
            raise ValueError(f"Incomplete four-reef group: {treatment}")
        two = group[group.reef_name.isin(TWO_WAVE)]
        exact = bool((group.sim_peak_count == group.expected_peak_count).all())
        matched = int(group.matched_peaks.sum()) == 7
        deep = bool((two.trough_to_smaller_peak <= 0.10).all())
        timing = bool((group.max_peak_year_error <= 4).all())
        amplitude = bool(group.all_heights_half_to_double.all())
        comparison.append({
            "treatment": treatment,
            "exact_peak_counts": exact,
            "seven_matched_peaks": matched,
            "all_three_troughs_below_0p10": deep,
            "max_trough_ratio": float(two.trough_to_smaller_peak.max())
                if two.trough_to_smaller_peak.notna().any() else None,
            "all_peak_years_within_4y": timing,
            "all_peak_heights_within_factor_2": amplitude,
            "qualitative_cycle_pass": exact and matched and deep and timing and amplitude,
            "mean_frozen_loss": float(group.frozen_loss.mean()),
            "minimum_coral": float(group.minimum_coral_1992_2024.min()),
            "matched_peak_total": int(group.matched_peaks.sum()),
        })
    summary = {
        "status": "bounded_qualitative_screen_not_promoted",
        "criteria": "Frozen 2/2/2/1 peaks, seven matched, all three trough ratios <=0.10, peak years within 4 years, heights within 0.5-2x. Coral separately reported because observation support differs.",
        "treatments": comparison,
        "inputs_sha256": {path.name: sha256(path) for path in
                          (run / "design.csv", run / "scores.csv",
                           run / "trajectories.csv", run / "flows.csv",
                           run / "failures.csv", GATE)},
        "script_sha256": sha256(Path(__file__)),
        "audit_sha256": sha256(output),
    }
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(pd.DataFrame(comparison).to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    main(parser.parse_args().run)
