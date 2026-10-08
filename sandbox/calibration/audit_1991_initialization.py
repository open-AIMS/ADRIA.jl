"""Audit per-reef 1991 initialization peaks and inter-wave troughs."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

TWO_WAVE_REEFS = {"Lizard Island Reef", "MacGillivray Reef", "North Direction Reef"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def audit(run: Path) -> None:
    output = run / "peak_trough_audit.csv"
    summary_path = run / "peak_trough_audit.json"
    if output.exists() or summary_path.exists():
        raise FileExistsError("Refusing to overwrite existing 1991 audit")
    scores = pd.read_csv(run / "scores.csv")
    trajectories = pd.read_csv(run / "trajectories.csv")
    initial = pd.read_csv(run / "initial_reef.csv")
    scores = scores[scores.observation_treatment == "smoothed_3y"]
    rows = []
    for score in scores.itertuples(index=False):
        series = trajectories[(trajectories.treatment == score.treatment) &
                              (trajectories.reef_name == score.reef_name) &
                              (trajectories.year >= 1992)].sort_values("year")
        if len(series) != 33:
            raise ValueError(f"Incomplete trajectory for {score.treatment} / {score.reef_name}")
        years = [int(y) for y in str(score.sim_peak_years).split(";") if y and y != "nan"]
        heights = [float(v) for v in str(score.sim_peak_heights_cpue).split(";")
                   if v and v != "nan"]
        if len(years) != len(heights):
            raise ValueError("Peak-year/height counts differ")
        smoothed = series.simulated_cpue.rolling(3, center=True, min_periods=1).mean()
        trough = np.nan
        ratio = np.nan
        if len(years) >= 2:
            between = series.year.between(years[0], years[1])
            trough = float(smoothed[between].min())
            ratio = trough / min(heights[:2]) if min(heights[:2]) > 0 else np.nan
        seed = initial[(initial.treatment == score.treatment) &
                       (initial.reef_name == score.reef_name)]
        if len(seed) != 1:
            raise ValueError("Missing reef initial-state summary")
        rows.append({"treatment": score.treatment, "reef_name": score.reef_name,
                     "first_peak_year": years[0] if years else np.nan,
                     "second_peak_year": years[1] if len(years) > 1 else np.nan,
                     "first_peak_cpue": heights[0] if heights else np.nan,
                     "second_peak_cpue": heights[1] if len(heights) > 1 else np.nan,
                     "smoothed_interpeak_min_cpue": trough,
                     "trough_to_smaller_peak": ratio,
                     "smoothed_loss": float(score.loss),
                     "matched_peaks": int(score.matched_peaks),
                     "initial_adults_ha": float(seed.seed_adults_ha.iloc[0]),
                     "initial_coral_fraction": float(seed.seed_corals_fraction.iloc[0]),
                     "minimum_model_coral_fraction_1992_2024": float(series.coral_cover.min())})
    frame = pd.DataFrame(rows)
    control = frame[frame.treatment == "inherited_1991"].set_index("reef_name")
    frame["first_peak_shift_vs_1991_control_years"] = [
        row.first_peak_year - control.loc[row.reef_name, "first_peak_year"]
        for row in frame.itertuples(index=False)]
    frame.to_csv(output, index=False)
    two_wave = frame[frame.reef_name.isin(TWO_WAVE_REEFS)]
    treatment_summary = two_wave.groupby("treatment").agg(
        max_trough_ratio=("trough_to_smaller_peak", "max"),
        mean_first_peak_shift_years=("first_peak_shift_vs_1991_control_years", "mean"),
        matched_peaks=("matched_peaks", "sum"))
    summary = {
        "status": "no_1991_treatment_passed_qualitative_gate",
        "criteria": "All three two-wave reefs need trough/smaller-adjacent-peak <= 0.1 without losing or delaying peaks or degrading coral",
        "treatments": treatment_summary.to_dict(orient="index"),
        "inputs": {name: sha256(run / name) for name in
                   ("scores.csv", "trajectories.csv", "initial_reef.csv")},
        "script_sha256": sha256(Path(__file__)),
        "output_sha256": sha256(output),
    }
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(treatment_summary.to_string())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    args = parser.parse_args()
    audit(args.run)
