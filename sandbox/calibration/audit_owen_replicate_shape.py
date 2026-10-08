"""Check replicated Owen peak troughs with parameters held fixed across seeds."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from audit_owen_peak_shape import peak_years, shape

ROOT = Path(__file__).resolve().parents[2]


def digest(path: Path) -> str:
    sha = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            sha.update(block)
    return sha.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-id", required=True)
    args = parser.parse_args()
    if not args.run_id.replace("_", "").replace("-", "").isalnum():
        raise ValueError("Invalid run ID")
    run = ROOT / "sandbox/calibration/runs" / args.run_id
    out = run / "shape_audit"
    out.mkdir(exist_ok=False)
    sources = {"scores": run / "scores.csv", "trajectories": run / "trajectories.csv"}
    scores = pd.read_csv(sources["scores"])
    trajectories = pd.read_csv(sources["trajectories"])
    rows = []
    for score in scores[scores.observation_treatment == "smoothed_3y"].itertuples():
        years = peak_years(score.sim_peak_years)
        if len(years) != 2:
            raise ValueError(f"Expected two peaks: {score.seed}, {score.treatment}, {score.reef_name}")
        series = trajectories[(trajectories.seed == score.seed) &
                              (trajectories.treatment == score.treatment) &
                              (trajectories.reef_name == score.reef_name)]
        if len(series) != 40:
            raise ValueError("Incomplete annual trajectory")
        smoothed = series.set_index("year").simulated_cpue.sort_index().rolling(
            3, center=True, min_periods=1).mean()
        rows.append({"seed": score.seed, "treatment": score.treatment,
                     "reef_name": score.reef_name, **shape(smoothed, *years)})
    result = pd.DataFrame(rows).sort_values(["seed", "treatment", "reef_name"])
    result.to_csv(out / "interpeak_troughs.csv", index=False)
    metadata = {
        "status": "diagnostic_only_not_part_of_frozen_score",
        "definition": "three-year-smoothed minimum between two detected peaks / smaller peak",
        "source_paths": {key: str(path) for key, path in sources.items()},
        "source_sha256": {key: digest(path) for key, path in sources.items()},
        "code_sha256": {str(path): digest(path) for path in
                        (Path(__file__), Path(__file__).with_name("audit_owen_peak_shape.py"))},
        "output_sha256": digest(out / "interpeak_troughs.csv"),
    }
    (out / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n",
                                       encoding="utf-8")
    print(result.groupby("treatment").trough_to_smaller_peak.agg(
        ["min", "median", "max"]).to_string())


if __name__ == "__main__":
    main()
