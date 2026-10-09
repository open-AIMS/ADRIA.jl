"""Compare opt-in size-weighted Owen trajectories with COTS and coral observations."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from plot_1991_cycle_capacity import ROOT, TARGETS, calendar_mean, sha256


def main(run: Path) -> None:
    output = run / "cots_coral_size_weighted_pair020.png"
    metadata_path = run / "plot_size_weighted_metadata.json"
    if output.exists() or metadata_path.exists():
        raise FileExistsError("Refusing to overwrite size-weighted plot")
    trajectories = pd.read_csv(run / "trajectories.csv")
    names = ("control", "size_fixed_growth", "size_food_growth",
             "size_food_growth_large75", "size_food_growth_d350")
    for name in names:
        if len(trajectories[trajectories.treatment == name]) != 4 * 34:
            raise ValueError(f"Incomplete trajectory: {name}")
    cots_path = ROOT / "sandbox/data/reef_cots.csv"
    coral_path = ROOT / "sandbox/data/reef_manta.csv"
    cots = pd.read_csv(cots_path)
    manta = pd.read_csv(coral_path)
    manta = manta[(manta.data_type == "manta") &
                  (manta.domain_category == "reef") &
                  (manta.variable == "HC") & (manta.purpose == "MANTA") &
                  (manta.project_code == "LTMP") & (manta.depth == 9)]
    styles = {
        "control": ("Unweighted control", "#64748b", "--"),
        "size_fixed_growth": ("Size, fixed growth", "#1d4ed8", "-"),
        "size_food_growth": ("Size, food-mediated growth", "#15803d", "-"),
        "size_food_growth_large75": ("Food growth, 75% large initially", "#b45309", ":"),
        "size_food_growth_d350": ("Food growth, large = 350 mm", "#7e22ce", ":"),
    }
    fig, axes = plt.subplots(4, 2, figsize=(14, 13), sharex=True)
    for row, (model_reef, observed_reef) in enumerate(TARGETS.items()):
        cots_ax, coral_ax = axes[row]
        obs_cots = cots[(cots.reef_name == observed_reef) &
                        cots.year.between(1992, 2024)].sort_values("year")
        obs_coral = manta[(manta.reef_name == observed_reef) &
                          manta.report_year.between(1992, 2024)]
        cots_ax.scatter(obs_cots.year, obs_cots.cotsptow,
                        color="#111827", s=16, alpha=0.5, zorder=5,
                        label="Observed COTS/tow" if row == 0 else None)
        cots_ax.plot(obs_cots.year, calendar_mean(obs_cots),
                     color="#111827", linewidth=1.8,
                     label="Observed 3-year mean" if row == 0 else None)
        coral_ax.scatter(obs_coral.report_year, obs_coral["median"],
                         color="#111827", s=16, alpha=0.55, zorder=5,
                         label="LTMP 9 m manta coral" if row == 0 else None)
        for name, (label, color, linestyle) in styles.items():
            series = trajectories[(trajectories.treatment == name) &
                                  (trajectories.reef_name == model_reef)].sort_values("year")
            cots_ax.plot(series.year.iloc[1:], series.simulated_cpue.iloc[1:],
                         color=color, linestyle=linestyle, linewidth=1.35,
                         label=label if row == 0 else None)
            coral_ax.plot(series.year, series.coral_cover,
                          color=color, linestyle=linestyle, linewidth=1.35)
        cots_ax.axhline(0.22, color="#9a3412", linewidth=0.8, linestyle=":")
        cots_ax.set_ylabel(f"{model_reef}\nCOTS/tow")
        coral_ax.set_ylabel("Hard-coral fraction")
        cots_ax.set_xlim(1991, 2024)
        coral_ax.set_xlim(1991, 2024)
        cots_ax.set_ylim(bottom=0)
        coral_ax.set_ylim(0, 1)
        cots_ax.grid(alpha=0.16)
        coral_ax.grid(alpha=0.16)
    axes[0, 0].set_title("Adult size and COTS peaks/troughs")
    axes[0, 1].set_title("Model whole-reef coral vs LTMP 9 m manta")
    axes[-1, 0].set_xlabel("Year")
    axes[-1, 1].set_xlabel("Year")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    handles2, labels2 = axes[0, 1].get_legend_handles_labels()
    fig.legend(handles + handles2, labels + labels2,
               loc="lower center", ncol=3, fontsize=9, bbox_to_anchor=(0.5, -0.005))
    fig.tight_layout(rect=(0, 0.055, 1, 1))
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)
    metadata = {
        "status": "one_seed_size_weighted_diagnostic_not_promoted",
        "series": list(names),
        "observation_note": "COTS black line is calendar-aware 3-year mean; coral black points are 9 m manta versus whole-reef model",
        "input_sha256": {name: sha256(path) for name, path in {
            "trajectories": run / "trajectories.csv",
            "cots_observations": cots_path,
            "coral_observations": coral_path,
            "script": Path(__file__),
        }.items()},
        "plot_sha256": sha256(output),
    }
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    print(output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    main(parser.parse_args().run)
