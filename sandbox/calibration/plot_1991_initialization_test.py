"""Compare COTS peaks and coral cover in the bounded 1991 initial-state screen."""

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
STYLES = {
    "inherited_1991": ("1991 inherited state", "#475569", "--"),
    "idw_alpha015_old_coral": ("IDW COTS, low conversion, old coral", "#d97706", "-"),
    "idw_alpha015_idw_coral": ("IDW COTS + coral, low conversion", "#dc2626", "-"),
    "idw_alpha150_old_coral": ("IDW COTS, high conversion, old coral", "#2563eb", "-"),
    "idw_alpha150_idw_coral": ("IDW COTS + coral, high conversion", "#15803d", "-"),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main(run: Path) -> None:
    if not run.is_dir():
        raise ValueError(f"Missing run directory: {run}")
    output = run / "cots_coral_1991_initialization.png"
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite {output}")
    trajectories = pd.read_csv(run / "trajectories.csv")
    if len(trajectories) != 5 * 4 * 34:
        raise ValueError("Incomplete 1991 trajectories")
    cots = pd.read_csv(ROOT / "sandbox" / "data" / "reef_cots.csv")
    manta = pd.read_csv(ROOT / "sandbox" / "data" / "reef_manta.csv")
    manta = manta[(manta.data_type == "manta") & (manta.domain_category == "reef") &
                  (manta.variable == "HC") & (manta.purpose == "MANTA") &
                  (manta.project_code == "LTMP") & (manta.depth == 9)]
    fig, axes = plt.subplots(4, 2, figsize=(13, 13), sharex=True,
                             gridspec_kw={"width_ratios": [1, 1]})
    for row, (model_reef, observed_reef) in enumerate(TARGETS.items()):
        cots_ax, coral_ax = axes[row]
        obs_cots = cots[(cots.reef_name == observed_reef) & cots.year.between(1992, 2024)]
        obs_coral = manta[(manta.reef_name == observed_reef) &
                          manta.report_year.between(1992, 2024)]
        cots_ax.scatter(obs_cots.year, obs_cots.cots / obs_cots.tows,
                        color="#111827", s=19, alpha=0.75, zorder=5,
                        label="LTMP COTS/tow" if row == 0 else None)
        coral_ax.scatter(obs_coral.report_year, obs_coral["median"],
                         color="#111827", s=19, alpha=0.65, zorder=5,
                         label="LTMP 9 m manta coral" if row == 0 else None)
        for name, (label, color, linestyle) in STYLES.items():
            series = trajectories[(trajectories.reef_name == model_reef) &
                                  (trajectories.treatment == name)].sort_values("year")
            cots_ax.plot(series.year.iloc[1:], series.simulated_cpue.iloc[1:],
                         color=color, linestyle=linestyle, linewidth=1.35,
                         label=label if row == 0 else None)
            coral_ax.plot(series.year, series.coral_cover, color=color,
                          linestyle=linestyle, linewidth=1.35)
        cots_ax.axhline(0.22, color="#92400e", linewidth=0.8, linestyle=":")
        cots_ax.set_ylabel(f"{model_reef}\nCOTS/tow")
        coral_ax.set_ylabel("Hard-coral fraction")
        cots_ax.grid(alpha=0.15)
        coral_ax.grid(alpha=0.15)
        cots_ax.set_xlim(1991, 2024)
        coral_ax.set_xlim(1991, 2024)
        cots_ax.set_ylim(bottom=0)
        coral_ax.set_ylim(0, 1)
    axes[0, 0].set_title("COTS: frozen observation scale; 1991 seed year excluded")
    axes[0, 1].set_title("Coral: model whole reef vs LTMP 9 m manta")
    axes[-1, 0].set_xlabel("Year")
    axes[-1, 1].set_xlabel("Year")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    handles2, labels2 = axes[0, 1].get_legend_handles_labels()
    fig.legend(handles + handles2, labels + labels2, loc="lower center", ncol=3,
               fontsize=8, bbox_to_anchor=(0.5, -0.01))
    fig.tight_layout(rect=(0, 0.035, 1, 1))
    fig.savefig(output, dpi=170, bbox_inches="tight")
    plt.close(fig)
    metadata = {"plot": output.name, "sha256": sha256(output),
                "inputs": {name: sha256(path) for name, path in {
                    "trajectories": run / "trajectories.csv",
                    "cots_observations": ROOT / "sandbox" / "data" / "reef_cots.csv",
                    "coral_observations": ROOT / "sandbox" / "data" / "reef_manta.csv",
                    "script": Path(__file__),
                }.items()}}
    (run / "plot_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    args = parser.parse_args()
    main(args.run)
