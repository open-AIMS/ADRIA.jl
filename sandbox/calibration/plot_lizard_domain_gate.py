"""Plot observed and replicated Lizard COTS peaks from a frozen gate run."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
REEFS = ["Lizard Island Reef", "MacGillivray Reef", "North Direction Reef", "Eyrie Reef"]
LEFT = [
    ("v1_legacy_closed", "V1 legacy / closed", "#64748b"),
    ("v2_legacy_closed", "V2 legacy / closed", "#2563eb"),
]
RIGHT = [
    ("v2_legacy_closed", "Legacy / closed", "#64748b"),
    ("v2_legacy_boundary", "Legacy / boundary", "#2563eb"),
    ("v2_owen_closed", "Owen / closed", "#d97706"),
    ("v2_owen_boundary", "Owen / boundary", "#15803d"),
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parts(value: object, cast):
    if pd.isna(value) or str(value) == "":
        return []
    return [cast(item) for item in str(value).split(";")]


def add_observed_segments(ax, observed):
    """Do not draw a continuous observation line across gaps over three years."""
    years = observed.year.to_numpy()
    values = observed.smoothed_3y_cpue.to_numpy()
    starts = np.r_[0, np.flatnonzero(np.diff(years) > 3) + 1]
    ends = np.r_[starts[1:], len(years)]
    for start, end in zip(starts, ends):
        ax.plot(years[start:end], values[start:end], color="#111827",
                linewidth=2.1, marker=".", markersize=5, zorder=7)


def add_ensemble(ax, trajectories, scores, reef, treatments):
    for treatment, label, color in treatments:
        selected = trajectories[(trajectories.reef_name == reef) &
                                (trajectories.treatment == treatment)]
        if selected.empty:
            raise ValueError(f"No trajectories for {reef}: {treatment}")
        summary = selected.groupby("year").simulated_cpue.agg(
            median="median", low=lambda x: x.quantile(0.1), high=lambda x: x.quantile(0.9)
        ).reset_index()
        ax.plot(summary.year, summary["median"], color=color, linewidth=1.8,
                label=label, zorder=3)
        ax.fill_between(summary.year.to_numpy(), summary.low.to_numpy(),
                        summary.high.to_numpy(), color=color, alpha=0.09, linewidth=0)
        peaks = scores[(scores.reef_name == reef) & (scores.treatment == treatment)]
        for row in peaks.itertuples():
            years = parts(row.sim_peak_years, int)
            heights = parts(row.sim_peak_heights_cpue, float)
            if len(years) != len(heights):
                raise ValueError(f"Peak years/heights mismatch for {reef}: {treatment}")
            ax.scatter(years, heights, marker="x", s=24, linewidths=0.9,
                       color=color, alpha=0.7, zorder=5)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-id", required=True)
    args = parser.parse_args()
    if not args.run_id.replace("_", "").replace("-", "").isalnum():
        raise ValueError("Invalid run ID")
    run_dir = REPO_ROOT / "sandbox" / "calibration" / "runs" / args.run_id
    paths = {name: run_dir / f"{name}.csv" for name in
             ("scores", "trajectories", "observations")}
    if any(not path.is_file() for path in paths.values()):
        raise FileNotFoundError("Comparison outputs are incomplete")
    scores = pd.read_csv(paths["scores"])
    trajectories = pd.read_csv(paths["trajectories"])
    observations = pd.read_csv(paths["observations"])
    scores = scores[scores.observation_treatment == "smoothed_3y"].copy()
    seeds = sorted(scores.seed.unique().tolist())
    if len(seeds) < 3:
        raise ValueError("At least three replicate seeds are required for this figure")
    if scores.groupby(["seed", "treatment", "reef_name"]).size().ne(1).any():
        raise ValueError("Duplicate score rows")
    if observations.groupby(["reef_name", "year"]).size().ne(1).any():
        raise ValueError("Duplicate observed reef-year rows")

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(4, 2, figsize=(15, 14), sharex=True,
                             constrained_layout=True)
    for row_idx, reef in enumerate(REEFS):
        obs = observations[observations.reef_name == reef].sort_values("year")
        if obs.empty:
            raise ValueError(f"No observations for {reef}")
        reef_scores = scores[scores.reef_name == reef]
        first = reef_scores.iloc[0]
        obs_peak_years = parts(first.obs_peak_years, int)
        obs_peak_heights = parts(first.obs_peak_heights_cpue, float)
        for col_idx, treatments in enumerate((LEFT, RIGHT)):
            ax = axes[row_idx, col_idx]
            ax.axhline(0.22, color="#b45309", linestyle="--", linewidth=1,
                       alpha=0.75, zorder=0)
            add_observed_segments(ax, obs)
            ax.scatter(obs.year, obs.raw_cpue, color="#374151", marker="o",
                       s=19, alpha=0.48, zorder=8)
            ax.scatter(obs_peak_years, obs_peak_heights, marker="*", s=85,
                       color="#111827", edgecolor="white", linewidth=0.45, zorder=9)
            add_ensemble(ax, trajectories, reef_scores, reef, treatments)
            ax.set_xlim(1985, 2024)
            ax.set_ylim(bottom=0)
            ax.grid(axis="y", color="#cbd5e1", linewidth=0.5, alpha=0.65)
            if col_idx == 0:
                ax.set_ylabel(f"{reef}\nCOTS/tow")
            if row_idx == 0:
                ax.set_title("Domain change: closed legacy control" if col_idx == 0
                             else "V2 transport and upstream boundary", fontsize=11)
                ax.legend(loc="upper right", frameon=True, framealpha=0.9,
                          fontsize=8, title="Model treatments", title_fontsize=8)
            if row_idx == 3:
                ax.set_xlabel("Calendar year")
        ymax = max(axes[row_idx, 0].get_ylim()[1], axes[row_idx, 1].get_ylim()[1])
        for ax in axes[row_idx, :]:
            ax.set_ylim(0, ymax)

    legend = [
        Line2D([0], [0], color="#111827", linewidth=2, label="Observed 3-year mean"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#374151",
               alpha=0.5, label="Raw survey CPUE"),
        Line2D([0], [0], marker="*", color="none", markerfacecolor="#111827",
               markersize=10, label="Detected observed peak"),
        Line2D([0], [0], marker="x", color="#2563eb", linestyle="none",
               label="Detected model peak (each seed)"),
        Line2D([0], [0], color="#b45309", linestyle="--", label="Outbreak: 0.22 COTS/tow"),
    ]
    fig.legend(handles=legend, loc="upper center", ncol=5,
               bbox_to_anchor=(0.5, 1.025), frameon=False, fontsize=8)
    fig.suptitle("Lizard COTS peak timing and amplitude | three matched seeds",
                 y=1.055, fontsize=14)
    fig.text(0.5, -0.012,
             "Lines: median model CPUE; shaded: 10–90% across seeds. Model peak markers use the scoring 3-year mean. "
             "The same V1-control CPUE scale is held for every treatment.",
             ha="center", va="top", fontsize=8, color="#475569")

    plot_dir = run_dir / "plots"
    plot_dir.mkdir(exist_ok=True)
    output = plot_dir / "cots_peak_comparison_gapaware.png"
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite existing plot: {output}")
    fig.savefig(output, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    metadata = {
        "plot": output.name,
        "script_sha256": sha256(Path(__file__)),
        "inputs": {name: sha256(path) for name, path in paths.items()},
        "output_sha256": sha256(output),
        "seeds": seeds,
        "observation_treatment": "smoothed_3y",
        "cpue_outbreak_threshold": 0.22,
    }
    metadata_path = plot_dir / "plot_metadata_gapaware.json"
    if metadata_path.exists():
        raise FileExistsError(f"Refusing to overwrite existing metadata: {metadata_path}")
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"Saved {output}")


if __name__ == "__main__":
    main()
