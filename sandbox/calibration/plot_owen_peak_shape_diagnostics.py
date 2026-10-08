"""Plot observed wave shape against fixed-seed parameter and boundary diagnostics."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
GATE = ROOT / "sandbox/calibration/runs/20261001T113401_lizard_domain_gate"
BASE = ROOT / "sandbox/calibration/runs/20261001T163000_owen_embedded_replicates"
SHUTDOWN = ROOT / "sandbox/calibration/runs/20261002T090000_owen_boundary_shutdown"
REEFS = ["Lizard Island Reef", "MacGillivray Reef",
         "North Direction Reef", "Eyrie Reef"]


def digest(path: Path) -> str:
    sha = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            sha.update(block)
    return sha.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-id", required=True,
                        help="Crash-parameter screen run ID; plot is written there")
    args = parser.parse_args()
    if not args.run_id.replace("_", "").replace("-", "").isalnum():
        raise ValueError("Invalid run ID")
    run = ROOT / "sandbox/calibration/runs" / args.run_id
    paths = {
        "parameter_trajectories": run / "trajectories.csv",
        "baseline_trajectories": BASE / "trajectories.csv",
        "shutdown_trajectories": SHUTDOWN / "trajectories.csv",
        "observations": GATE / "observations.csv",
    }
    parameter = pd.read_csv(paths["parameter_trajectories"])
    baseline = pd.read_csv(paths["baseline_trajectories"])
    baseline = baseline[(baseline.seed == 20260930) &
                        (baseline.treatment == "embedded_joint_4x")]
    shutdown = pd.read_csv(paths["shutdown_trajectories"])
    observations = pd.read_csv(paths["observations"])
    out = run / "plots"
    out.mkdir(exist_ok=False)

    lines = [
        ("Joint 4× baseline", baseline, "#7e22ce", "-"),
        ("No external supply, 2000–10", shutdown, "#d97706", "--"),
        ("Adult mortality upper", parameter[parameter.treatment == "adult_mortality_upper"],
         "#dc2626", "-"),
        ("Ricker dependence upper", parameter[parameter.treatment == "density_dependence_upper"],
         "#2563eb", "-"),
        ("Both upper", parameter[parameter.treatment == "both_upper"],
         "#15803d", "-"),
    ]
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.spines.top": False, "axes.spines.right": False})
    figure, axes = plt.subplots(4, 1, figsize=(12.5, 12.8), sharex=True,
                               constrained_layout=True)
    for idx, (axis, reef) in enumerate(zip(axes, REEFS)):
        obs = observations[observations.reef_name == reef].sort_values("year")
        axis.plot(obs.year, obs.smoothed_3y_cpue, color="#111827", linewidth=2,
                  label="LTMP 3-year mean" if idx == 0 else None)
        axis.scatter(obs.year, obs.raw_cpue, color="#111827", s=12, alpha=0.32,
                     label="LTMP annual tow" if idx == 0 else None)
        for label, data, color, style in lines:
            subset = data[data.reef_name == reef].sort_values("year")
            if len(subset) != 40:
                raise ValueError(f"Incomplete trajectory for {label}, {reef}")
            axis.plot(subset.year, subset.simulated_cpue, color=color,
                      linestyle=style, linewidth=1.35,
                      label=label if idx == 0 else None)
        axis.axhline(0.22, color="#92400e", linestyle=":", linewidth=1,
                     label="0.22 COTS/tow" if idx == 0 else None)
        axis.set_ylabel(f"{reef}\nCOTS/tow")
        axis.set_ylim(bottom=0)
        axis.grid(axis="y", alpha=0.25)
    axes[0].legend(loc="upper left", ncol=4, fontsize=8, frameon=False)
    axes[0].set_title("Do existing parameters reproduce the outbreak crash? (seed 20260930)")
    axes[-1].set_xlabel("Report/model year")
    filename = out / "cots_peak_shape_diagnostics.png"
    figure.savefig(filename, dpi=175)
    plt.close(figure)
    metadata = {
        "scope": "diagnostic only; boundary shutdown is an oracle counterfactual",
        "source_paths": {key: str(path) for key, path in paths.items()},
        "source_sha256": {key: digest(path) for key, path in paths.items()},
        "script_sha256": digest(Path(__file__)),
        "output_sha256": digest(filename),
    }
    (out / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n",
                                       encoding="utf-8")
    print(filename)


if __name__ == "__main__":
    main()
