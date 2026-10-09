"""Plot one paired capacity treatment with observed COTS and coral trajectories."""

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
TARGETS = {
    "Lizard Island Reef": "Lizard Isles",
    "MacGillivray Reef": "Macgillivray Reef",
    "North Direction Reef": "North Direction Island",
    "Eyrie Reef": "Eyrie Reef",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def calendar_mean(frame: pd.DataFrame) -> list[float]:
    years = frame.year.to_numpy()
    cpue = frame.cotsptow.to_numpy(dtype=float)
    return [float(cpue[abs(years - year) <= 1].mean()) for year in years]


def main(run: Path, pair: int) -> None:
    output = run / f"cots_coral_capacity_pair_{pair:03d}.png"
    metadata_path = run / f"plot_capacity_pair_{pair:03d}_metadata.json"
    if output.exists() or metadata_path.exists():
        raise FileExistsError("Refusing to overwrite capacity plot")
    trajectories = pd.read_csv(run / "trajectories.csv")
    names = ("pair_000_theta_3", f"pair_{pair:03d}_theta_3",
             f"pair_{pair:03d}_theta_1")
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
        names[0]: ("1991 control, 3/ha", "#64748b", "--"),
        names[1]: (f"Pair {pair}, 3/ha", "#1d4ed8", "-"),
        names[2]: (f"Pair {pair}, 1/ha counterfactual", "#d97706", "-"),
    }
    fig, axes = plt.subplots(4, 2, figsize=(13, 13), sharex=True)
    for row, (model_reef, observed_reef) in enumerate(TARGETS.items()):
        cots_ax, coral_ax = axes[row]
        obs_cots = cots[(cots.reef_name == observed_reef) &
                        cots.year.between(1992, 2024)].sort_values("year")
        obs_coral = manta[(manta.reef_name == observed_reef) &
                          manta.report_year.between(1992, 2024)]
        cots_ax.scatter(obs_cots.year, obs_cots.cotsptow,
                        color="#111827", s=18, alpha=0.5, zorder=5,
                        label="Observed COTS/tow" if row == 0 else None)
        cots_ax.plot(obs_cots.year, calendar_mean(obs_cots),
                     color="#111827", linewidth=1.8,
                     label="Observed 3-year mean" if row == 0 else None)
        coral_ax.scatter(obs_coral.report_year, obs_coral["median"],
                         color="#111827", s=18, alpha=0.55, zorder=5,
                         label="LTMP 9 m manta coral" if row == 0 else None)
        for name, (label, color, linestyle) in styles.items():
            series = trajectories[(trajectories.treatment == name) &
                                  (trajectories.reef_name == model_reef)].sort_values("year")
            cots_ax.plot(series.year.iloc[1:], series.simulated_cpue.iloc[1:],
                         color=color, linestyle=linestyle, linewidth=1.5,
                         label=label if row == 0 else None)
            coral_ax.plot(series.year, series.coral_cover,
                          color=color, linestyle=linestyle, linewidth=1.5)
        cots_ax.axhline(0.22, color="#9a3412", linewidth=0.8, linestyle=":")
        cots_ax.set_ylabel(f"{model_reef}\nCOTS/tow")
        coral_ax.set_ylabel("Hard-coral fraction")
        cots_ax.set_xlim(1991, 2024)
        coral_ax.set_xlim(1991, 2024)
        cots_ax.set_ylim(bottom=0)
        coral_ax.set_ylim(0, 1)
        cots_ax.grid(alpha=0.16)
        coral_ax.grid(alpha=0.16)
    axes[0, 0].set_title("COTS peaks and inter-wave troughs")
    axes[0, 1].set_title("Model whole-reef coral vs LTMP 9 m manta")
    axes[-1, 0].set_xlabel("Year")
    axes[-1, 1].set_xlabel("Year")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    handles2, labels2 = axes[0, 1].get_legend_handles_labels()
    fig.legend(handles + handles2, labels + labels2,
               loc="lower center", ncol=3, fontsize=9, bbox_to_anchor=(0.5, -0.005))
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)
    metadata = {
        "status": "incomplete_screen_illustration_not_promoted",
        "pair": pair,
        "series": list(names),
        "observation_note": "COTS black line is calendar-aware 3-year mean; coral black points are 9m manta versus whole-reef model",
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
    parser.add_argument("--pair", type=int, default=20)
    args = parser.parse_args()
    main(args.run, args.pair)
