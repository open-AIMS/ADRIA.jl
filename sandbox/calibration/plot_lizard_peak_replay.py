"""Plot fixed-scale COTS and coral cover for the bounded V1 peak replay."""

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
REEFS = ["Lizard Island Reef", "MacGillivray Reef", "North Direction Reef", "Eyrie Reef"]
TREATMENTS = [
    ("old_coral_allee1_counterfactual", "Old candidate / coral / Allee 1 (counterfactual)", "#dc2626"),
    ("old_coral_allee3", "Old candidate / coral / Allee 3", "#d97706"),
    ("old_cots_allee3", "Old candidate / COTS / Allee 3", "#2563eb"),
    ("current_cots_legacy_allee3", "Current candidate / COTS / Allee 3", "#475569"),
    ("current_cots_separated_survival_allee3", "Current separated + survival / Allee 3", "#15803d"),
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-id", required=True)
    args = parser.parse_args()
    if not args.run_id.replace("_", "").replace("-", "").isalnum():
        raise ValueError("Invalid run ID")
    run_dir = ROOT / "sandbox/calibration/runs" / args.run_id
    paths = {
        "replay": run_dir / "trajectories.csv",
        "observations": ROOT / "sandbox/calibration/runs/20261001T113401_lizard_domain_gate/observations.csv",
        "empirical_coral": ROOT / "sandbox/data/calibration_results_emp_coral.csv",
        "archive": ROOT / "sandbox/data/best_bbo_simulated_trajectories.csv",
    }
    for path in paths.values():
        if not path.is_file():
            raise FileNotFoundError(path)
    replay = pd.read_csv(paths["replay"])
    observed = pd.read_csv(paths["observations"])
    empirical = pd.read_csv(paths["empirical_coral"])
    archive = pd.read_csv(paths["archive"])
    empirical = empirical[(empirical.variable == "HC") & (empirical.data_type == "manta") &
                          (empirical.reef_name.isin(REEFS)) &
                          (empirical.report_year.between(1985, 2024))]
    if replay.groupby(["treatment", "reef_name", "year"]).size().ne(1).any():
        raise ValueError("Duplicate replay trajectories")
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(4, 2, figsize=(15, 13), sharex=True,
                             constrained_layout=True)
    for row, reef in enumerate(REEFS):
        cots_ax, coral_ax = axes[row]
        obs = observed[observed.reef_name == reef].sort_values("year")
        cots_ax.scatter(obs.year, obs.raw_cpue, s=18, alpha=0.4, color="#111827",
                        label="Survey COTS/tow" if row == 0 else None, zorder=7)
        cots_ax.plot(obs.year, obs.smoothed_3y_cpue, linestyle="none", marker=".",
                     markersize=4, color="#111827", label="Survey 3-year mean" if row == 0
                     else None, zorder=8)
        cots_ax.axhline(0.22, color="#92400e", linestyle="--", linewidth=1,
                        label="0.22 outbreak threshold" if row == 0 else None)
        for treatment, label, color in TREATMENTS:
            subset = replay[(replay.reef_name == reef) & (replay.treatment == treatment)]
            if len(subset) != 40:
                raise ValueError(f"Missing or incomplete {reef}: {treatment}")
            subset = subset.sort_values("year")
            cots_ax.plot(subset.year, subset.simulated_cpue, color=color, linewidth=1.7,
                         label=label if row == 0 else None)
            coral_ax.plot(subset.year, subset.coral_cover, color=color, linewidth=1.7,
                          label=label if row == 0 else None)
        old = archive[archive.reef_name == reef].sort_values("year")
        coral_ax.plot(old.year, old.sim_coral_cover, color="#7c3aed", linestyle=":",
                      linewidth=1.4, label="Archived cover (unweighted)" if row == 0
                      else None)
        e = empirical[empirical.reef_name == reef].sort_values("report_year")
        if not e.empty:
            coral_ax.errorbar(e.report_year, e["median"],
                              yerr=[e["median"] - e.lower, e.upper - e["median"]],
                              fmt="o", markersize=4, color="#111827", capsize=2,
                              linewidth=0.8, label="Eyrie LTMP manta HC (9 m)", zorder=8)
        else:
            coral_ax.text(0.02, 0.95, "No exact-match coral survey series",
                          transform=coral_ax.transAxes, va="top", fontsize=8,
                          color="#6b7280")
        cots_ax.set_ylabel(f"{reef}\nCOTS/tow")
        coral_ax.set_ylabel("Hard coral cover (fraction)")
        cots_ax.set_ylim(bottom=0)
        coral_ax.set_ylim(0, 0.65)
        for ax in (cots_ax, coral_ax):
            ax.set_xlim(1985, 2024)
            ax.grid(axis="y", color="#cbd5e1", linewidth=0.5, alpha=0.6)
        if row == 0:
            cots_ax.set_title("Old-wave replay: fixed-scale COTS/tow")
            coral_ax.set_title("Matched coral-cover trajectories")
        if row == 3:
            cots_ax.set_xlabel("Calendar year")
            coral_ax.set_xlabel("Calendar year")
    handles, labels = [], []
    for ax in (axes[0, 0], axes[0, 1], axes[3, 1]):
        h, l = ax.get_legend_handles_labels()
        handles.extend(h)
        labels.extend(l)
    unique = dict(zip(labels, handles))
    fig.legend(unique.values(), unique.keys(), loc="upper center", ncol=3,
               bbox_to_anchor=(0.5, 1.06), frameon=False, fontsize=8)
    fig.suptitle("What restores the old Lizard waves?", y=1.10, fontsize=14)
    fig.text(0.5, -0.014,
             "All replay model lines are V1 polygon-area-weighted reef means with a frozen CPUE scale. "
             "Only the red line lowers Allee to 1 COTS/ha; it is a historical counterfactual, "
             "not a proposed biological threshold. Archived cover is unweighted. Eyrie 9 m survey "
             "cover is not a whole-reef estimate.",
             ha="center", va="top", fontsize=8, color="#475569", wrap=True)
    plot_dir = run_dir / "plots"
    plot_dir.mkdir(exist_ok=True)
    output = plot_dir / "cots_coral_peak_replay.png"
    metadata = plot_dir / "cots_coral_peak_replay.json"
    if output.exists() or metadata.exists():
        raise FileExistsError("Replay plot or metadata already exists")
    fig.savefig(output, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    metadata.write_text(json.dumps({
        "script_sha256": sha256(Path(__file__)),
        "inputs": {name: sha256(path) for name, path in paths.items()},
        "output_sha256": sha256(output),
        "observation_mapping": "frozen V1 polygon-area-weighted reef density and CPUE scale",
        "counterfactual": "1 COTS/ha Allee; all other replay lines 3 COTS/ha",
    }, indent=2) + "\n", encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()
