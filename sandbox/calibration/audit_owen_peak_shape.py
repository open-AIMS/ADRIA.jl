"""Describe inter-peak trough depth without changing the frozen COTS score."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
GATE = ROOT / "sandbox/calibration/runs/20261001T113401_lizard_domain_gate"


def digest(path: Path) -> str:
    sha = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            sha.update(block)
    return sha.hexdigest()


def peak_years(value: str | float) -> list[int]:
    if pd.isna(value) or str(value).strip() == "":
        return []
    return [int(item) for item in str(value).split(";")]


def shape(series: pd.Series, first: int, second: int) -> dict[str, float | int]:
    between = series[(series.index > first) & (series.index < second)]
    if between.empty:
        raise ValueError("No measured year between peaks")
    low_year = int(between.idxmin())
    low = float(between.min())
    smaller_peak = min(float(series.loc[first]), float(series.loc[second]))
    return {
        "first_peak_year": first,
        "second_peak_year": second,
        "first_peak_cpue": float(series.loc[first]),
        "second_peak_cpue": float(series.loc[second]),
        "trough_year": low_year,
        "trough_cpue": low,
        "trough_to_smaller_peak": low / smaller_peak if smaller_peak > 0 else np.nan,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-id", required=True)
    args = parser.parse_args()
    if not args.run_id.replace("_", "").replace("-", "").isalnum():
        raise ValueError("Invalid run ID")
    run = ROOT / "sandbox/calibration/runs" / args.run_id
    out = run / "shape_audit"
    out.mkdir(exist_ok=False)
    sources = {
        "scores": run / "scores.csv",
        "trajectories": run / "trajectories.csv",
        "gate_scores": GATE / "scores.csv",
        "gate_trajectories": GATE / "trajectories.csv",
        "observations": GATE / "observations.csv",
    }
    scores = pd.read_csv(sources["scores"])
    scores = scores[scores.observation_treatment == "smoothed_3y"]
    trajectories = pd.read_csv(sources["trajectories"])
    gate_scores = pd.read_csv(sources["gate_scores"])
    gate_scores = gate_scores[(gate_scores.seed == 20260930) &
                              (gate_scores.treatment == "v2_owen_boundary") &
                              (gate_scores.observation_treatment == "smoothed_3y")]
    gate_trajectories = pd.read_csv(sources["gate_trajectories"])
    gate_trajectories = gate_trajectories[(gate_trajectories.seed == 20260930) &
                                          (gate_trajectories.treatment == "v2_owen_boundary")]
    observations = pd.read_csv(sources["observations"])

    rows: list[dict[str, object]] = []
    for treatment, reef, score, data in (
        [("observed", str(row.reef_name), row, observations[
            observations.reef_name == row.reef_name]) for row in gate_scores.itertuples()]
        + [("v2_owen_boundary", str(row.reef_name), row, gate_trajectories[
            gate_trajectories.reef_name == row.reef_name]) for row in gate_scores.itertuples()]
        + [(str(row.treatment), str(row.reef_name), row, trajectories[
            (trajectories.treatment == row.treatment) &
            (trajectories.reef_name == row.reef_name)]) for row in scores.itertuples()]
    ):
        years = peak_years(score.obs_peak_years if treatment == "observed"
                           else score.sim_peak_years)
        if len(years) < 2:
            continue
        if treatment == "observed":
            series = data.set_index("year").smoothed_3y_cpue.sort_index()
        else:
            series = data.set_index("year").simulated_cpue.sort_index()
            series = series.rolling(window=3, center=True, min_periods=1).mean()
        if len(years) != 2:
            raise ValueError(f"Expected at most two peaks: {treatment}, {reef}")
        if not all(year in series.index for year in years):
            raise ValueError(f"Missing detected peak year: {treatment}, {reef}")
        rows.append({"treatment": treatment, "reef_name": reef,
                     **shape(series, years[0], years[1])})
    summary = pd.DataFrame(rows).sort_values(["reef_name", "treatment"])
    summary.to_csv(out / "interpeak_troughs.csv", index=False)
    metadata = {
        "status": "diagnostic_only_not_part_of_frozen_score",
        "definition": "minimum 3-year-smoothed COTS/tow strictly between first two detected peaks, divided by smaller peak",
        "observed_note": "LTMP smoothed survey years only; model uses complete annual sequence",
        "source_paths": {key: str(path) for key, path in sources.items()},
        "source_sha256": {key: digest(path) for key, path in sources.items()},
        "script_sha256": digest(Path(__file__)),
        "output_sha256": digest(out / "interpeak_troughs.csv"),
    }
    (out / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n",
                                       encoding="utf-8")
    print(summary[["treatment", "reef_name", "first_peak_year", "second_peak_year",
                   "trough_cpue", "trough_to_smaller_peak"]].to_string(index=False))


if __name__ == "__main__":
    main()
