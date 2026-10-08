"""Plot fixed-scale COTS peaks and coral guardrails for the Owen assumption screen."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from plot_lizard_ltmp_coral import REEF_ALIASES, ROOT, load_ltmp_coral

GATE = ROOT / "sandbox/calibration/runs/20261001T113401_lizard_domain_gate"
SERIES = [
    ("v2_owen_boundary", "RME-on control", "#64748b"),
    ("embedded_candidate", "Embedded: candidate production", "#0e7490"),
    ("embedded_mean_matched", "Embedded: mean-matched", "#2563eb"),
    ("embedded_local_2x", "Embedded: local 2×", "#d97706"),
    ("embedded_boundary_2x", "Embedded: boundary 2×", "#15803d"),
    ("embedded_joint_2x", "Embedded: both 2×", "#dc2626"),
    ("embedded_joint_4x", "Embedded: both 4×", "#7e22ce"),
]


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
    paths = {
        "screen": run / "trajectories.csv",
        "gate": GATE / "trajectories.csv",
        "observations": GATE / "observations.csv",
        "manta": ROOT / "sandbox/data/reef_manta.csv",
        "photo": ROOT / "sandbox/data/reef_photo_transect.csv",
    }
    screen = pd.read_csv(paths["screen"])
    gate = pd.read_csv(paths["gate"])
    gate = gate[(gate.seed == 20260930) & (gate.treatment == "v2_owen_boundary")]
    gate = gate[["treatment", "reef_name", "year", "simulated_cpue", "coral_cover"]]
    trajectories = pd.concat([screen, gate], ignore_index=True)
    if trajectories.groupby(["treatment", "reef_name", "year"]).size().ne(1).any():
        raise ValueError("Duplicate treatment/reef/year trajectory")
    observations = pd.read_csv(paths["observations"])
    coral = load_ltmp_coral(paths["manta"], paths["photo"])
    out = run / "plots"
    out.mkdir(exist_ok=False)

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.spines.top": False, "axes.spines.right": False})
    figure, axes = plt.subplots(4, 1, figsize=(12.5, 13), sharex=True,
                               constrained_layout=True)
    for axis, reef in zip(axes, REEF_ALIASES):
        obs = observations[observations.reef_name == reef].sort_values("year")
        axis.scatter(obs.year, obs.raw_cpue, s=13, color="#111827", alpha=0.35,
                     label="LTMP COTS/tow" if reef == next(iter(REEF_ALIASES)) else None)
        axis.plot(obs.year, obs.smoothed_3y_cpue, color="#111827", linewidth=1.4,
                  label="LTMP 3-year mean" if reef == next(iter(REEF_ALIASES)) else None)
        for name, label, color in SERIES:
            data = trajectories[(trajectories.reef_name == reef) &
                                (trajectories.treatment == name)].sort_values("year")
            if len(data) != 40:
                raise ValueError(f"Incomplete trajectory: {reef}, {name}")
            axis.plot(data.year, data.simulated_cpue, color=color,
                      linewidth=1.25 if name != "v2_owen_boundary" else 1.75,
                      label=label if reef == next(iter(REEF_ALIASES)) else None)
        axis.axhline(0.22, color="#92400e", linestyle="--", linewidth=0.9,
                     label="0.22 COTS/tow" if reef == next(iter(REEF_ALIASES)) else None)
        axis.set_ylabel(f"{reef}\nCOTS/tow")
        axis.set_ylim(bottom=0)
        axis.grid(axis="y", alpha=0.25)
    axes[-1].set_xlabel("Report/model year")
    axes[0].legend(loc="upper left", ncol=4, fontsize=8, frameon=False)
    axes[0].set_title("Owen embedded-mortality assumption: peak capacity (seed 20260930)")
    figure.savefig(out / "cots_peak_capacity.png", dpi=175)
    plt.close(figure)

    coral_series = [SERIES[i] for i in (0, 2, 5, 6)]
    figure, axes = plt.subplots(4, 1, figsize=(12.5, 12), sharex=True,
                               constrained_layout=True)
    for axis, reef in zip(axes, REEF_ALIASES):
        for name, label, color in coral_series:
            data = trajectories[(trajectories.reef_name == reef) &
                                (trajectories.treatment == name)].sort_values("year")
            axis.plot(data.year, data.coral_cover, color=color, linewidth=1.6,
                      label=label if reef == next(iter(REEF_ALIASES)) else None)
        for method, marker, color, label in (
            ("manta", "o", "#111827", "LTMP manta, 9 m"),
            ("photo_transect", "D", "#a21caf", "LTMP photo, 9 m"),
        ):
            obs = coral[(coral.model_reef_name == reef) & (coral.method == method)]
            if not obs.empty:
                axis.scatter(obs.year, obs["median"], marker=marker, s=19,
                             color=color, alpha=0.75,
                             label=label if reef == next(iter(REEF_ALIASES)) else None)
        axis.set_ylabel(f"{reef}\nHard coral cover")
        axis.set_ylim(0, 0.7)
        axis.grid(axis="y", alpha=0.25)
    axes[-1].set_xlabel("Report/model year")
    axes[0].legend(loc="upper left", ncol=3, fontsize=8, frameon=False)
    axes[0].set_title("Coral guardrail: LTMP 9 m surveys versus whole-reef model")
    figure.savefig(out / "coral_guardrail.png", dpi=175)
    plt.close(figure)

    metadata = {
        "scope": "diagnostic plots; coral survey and model spatial supports differ",
        "source_paths": {key: str(path) for key, path in paths.items()},
        "source_sha256": {key: digest(path) for key, path in paths.items()},
        "plot_script_sha256": digest(Path(__file__)),
        "output_sha256": {name: digest(out / name) for name in
                          ("cots_peak_capacity.png", "coral_guardrail.png")},
    }
    (out / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(out)


if __name__ == "__main__":
    main()
