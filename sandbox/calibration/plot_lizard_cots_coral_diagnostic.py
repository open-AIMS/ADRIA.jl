"""Compare COTS peaks and coral cover from one immutable Lizard gate run.

Coral observations are plotted only where the archived LTMP manta series has an
exact calibration-reef name. They are not treated as an area-weighted reef truth.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
REEFS = ["Lizard Island Reef", "MacGillivray Reef", "North Direction Reef", "Eyrie Reef"]
TREATMENTS = [
    ("v1_legacy_closed", "V1 legacy closed", "#64748b"),
    ("v2_legacy_closed", "V2 legacy closed", "#2563eb"),
    ("v2_owen_closed", "V2 Owen closed", "#d97706"),
    ("v2_owen_boundary", "V2 Owen boundary", "#15803d"),
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
    inputs = {
        "trajectories": run_dir / "trajectories.csv",
        "observations": run_dir / "observations.csv",
        "empirical_coral": ROOT / "sandbox/data/calibration_results_emp_coral.csv",
        "archived_trajectory": ROOT / "sandbox/data/best_bbo_simulated_trajectories.csv",
    }
    for path in inputs.values():
        if not path.is_file():
            raise FileNotFoundError(path)
    trajectories = pd.read_csv(inputs["trajectories"])
    observed = pd.read_csv(inputs["observations"])
    coral_obs = pd.read_csv(inputs["empirical_coral"])
    archive = pd.read_csv(inputs["archived_trajectory"])
    required = {"reef_name", "year", "coral_cover", "simulated_cpue", "treatment", "seed"}
    if not required.issubset(trajectories):
        raise ValueError(f"Missing trajectory columns: {required - set(trajectories)}")
    if trajectories.groupby(["reef_name", "year", "treatment", "seed"]).size().ne(1).any():
        raise ValueError("Duplicate model trajectory rows")
    coral_obs = coral_obs[(coral_obs.variable == "HC") &
                            (coral_obs.data_type == "manta") &
                            (coral_obs.reef_name.isin(REEFS)) &
                            (coral_obs.report_year.between(1985, 2024))].copy()
    if coral_obs.groupby(["reef_name", "report_year"]).size().gt(1).any():
        raise ValueError("Ambiguous empirical coral reef-year rows")

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(4, 2, figsize=(15, 13), sharex=True,
                             constrained_layout=True)
    for row, reef in enumerate(REEFS):
        cots_ax, coral_ax = axes[row]
        reef_obs = observed[observed.reef_name == reef].sort_values("year")
        cots_ax.scatter(reef_obs.year, reef_obs.raw_cpue, s=17, color="#111827",
                        alpha=0.45, label="Observed COTS/tow" if row == 0 else None, zorder=6)
        cots_ax.plot(reef_obs.year, reef_obs.smoothed_3y_cpue, color="#111827",
                     linewidth=1.5, linestyle="none", marker=".", markersize=4,
                     label="Observed 3-year mean" if row == 0 else None, zorder=7)
        cots_ax.axhline(0.22, color="#b45309", linestyle="--", linewidth=1,
                        label="Outbreak threshold" if row == 0 else None)
        for treatment, label, color in TREATMENTS:
            subset = trajectories[(trajectories.reef_name == reef) &
                                  (trajectories.treatment == treatment)]
            if subset.empty:
                raise ValueError(f"Missing {reef}, {treatment}")
            for ax, column in ((cots_ax, "simulated_cpue"), (coral_ax, "coral_cover")):
                summary = subset.groupby("year")[column].agg(
                    median="median", low=lambda x: x.quantile(0.1),
                    high=lambda x: x.quantile(0.9)).reset_index()
                ax.plot(summary.year, summary["median"], color=color, linewidth=1.65,
                        label=label if row == 0 else None)
                ax.fill_between(summary.year.to_numpy(), summary.low.to_numpy(),
                                summary.high.to_numpy(), color=color, alpha=0.09,
                                linewidth=0)
        old = archive[archive.reef_name == reef].sort_values("year")
        if not old.empty:
            coral_ax.plot(old.year, old.sim_coral_cover, color="#7c3aed",
                          linestyle=":", linewidth=1.5,
                          label="Archived old model (unweighted)" if row == 0 else None)
        empirical = coral_obs[coral_obs.reef_name == reef].sort_values("report_year")
        if not empirical.empty:
            coral_ax.errorbar(empirical.report_year, empirical["median"],
                              yerr=[empirical["median"] - empirical.lower,
                                    empirical.upper - empirical["median"]],
                              fmt="o", markersize=4, color="#111827", capsize=2,
                              linewidth=0.8, label="LTMP manta HC (9 m)" if row == 0
                              else "LTMP manta HC (9 m)", zorder=8)
        else:
            coral_ax.text(0.02, 0.95, "No exact-match coral survey series",
                          transform=coral_ax.transAxes, va="top", fontsize=8,
                          color="#6b7280")
        cots_ax.set_ylabel(f"{reef}\nCOTS/tow")
        coral_ax.set_ylabel("Hard coral cover (fraction)")
        cots_ax.set_ylim(bottom=0)
        coral_ax.set_ylim(0, max(0.55, coral_ax.get_ylim()[1]))
        for ax in (cots_ax, coral_ax):
            ax.set_xlim(1985, 2024)
            ax.grid(axis="y", color="#cbd5e1", alpha=0.6, linewidth=0.5)
        if row == 0:
            cots_ax.set_title("COTS peaks: survey CPUE and model")
            coral_ax.set_title("Coral food / response: model and available survey")
        if row == 3:
            cots_ax.set_xlabel("Calendar year")
            coral_ax.set_xlabel("Calendar year")

    handles, labels = [], []
    for ax in (axes[0, 0], axes[3, 1]):
        h, l = ax.get_legend_handles_labels()
        handles.extend(h)
        labels.extend(l)
    unique = dict(zip(labels, handles))
    fig.legend(unique.values(), unique.keys(), loc="upper center", ncol=4,
               bbox_to_anchor=(0.5, 1.055), frameon=False, fontsize=8)
    fig.suptitle("Lizard COTS peaks alongside coral cover", y=1.085, fontsize=14)
    fig.text(0.5, -0.014,
             "Model lines: area-weighted reef median across three transport seeds; shade: 10–90%. "
             "Old model line: unweighted archived site mean, not a like-for-like control. "
             "Eyrie LTMP manta hard coral (9 m) is a survey proxy, not an area-weighted reef estimate.",
             ha="center", va="top", fontsize=8, color="#475569", wrap=True)
    out_dir = run_dir / "plots"
    out_dir.mkdir(exist_ok=True)
    output = out_dir / "cots_coral_cover_diagnostic.png"
    metadata_path = out_dir / "cots_coral_cover_diagnostic.json"
    if output.exists() or metadata_path.exists():
        raise FileExistsError("Diagnostic plot or metadata already exists; use a new run")
    fig.savefig(output, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    metadata_path.write_text(json.dumps({
        "script_sha256": sha256(Path(__file__)),
        "inputs": {name: sha256(path) for name, path in inputs.items()},
        "output_sha256": sha256(output),
        "coral_observation_reefs": sorted(coral_obs.reef_name.unique().tolist()),
        "coral_observation_limit": "LTMP manta HC at 9 m; not area-weighted reef cover",
        "archived_line_limit": "archived model cover is an unweighted site mean",
    }, indent=2).encode("utf-8").decode("utf-8") + "\n", encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()
