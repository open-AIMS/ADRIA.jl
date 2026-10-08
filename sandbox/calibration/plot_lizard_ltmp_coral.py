"""Plot Lizard COTS and model coral with mapped LTMP manta/photo observations.

This is an observation revision of the earlier coral diagnostics. It does not
change model trajectories, the CPUE conversion, or the frozen COTS score.
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
GATE_RUN = "20261001T113401_lizard_domain_gate"
REEF_ALIASES = {
    "Lizard Island Reef": "Lizard Isles",
    "MacGillivray Reef": "Macgillivray Reef",
    "North Direction Reef": "North Direction Island",
    "Eyrie Reef": "Eyrie Reef",
}
TREATMENTS = {
    "gate": [
        ("v1_legacy_closed", "V1 legacy closed", "#64748b"),
        ("v2_legacy_closed", "V2 legacy closed", "#2563eb"),
        ("v2_owen_closed", "V2 Owen closed", "#d97706"),
        ("v2_owen_boundary", "V2 Owen boundary", "#15803d"),
    ],
    "replay": [
        ("old_coral_allee1_counterfactual", "Old/coral/Allee 1 (counterfactual)", "#dc2626"),
        ("old_coral_allee3", "Old/coral/Allee 3", "#d97706"),
        ("old_cots_allee3", "Old/COTS/Allee 3", "#2563eb"),
        ("current_cots_legacy_allee3", "Current/COTS/Allee 3", "#475569"),
        ("current_cots_separated_survival_allee3", "Separated+survival/Allee 3", "#15803d"),
    ],
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_ltmp_coral(manta_path: Path, photo_path: Path) -> pd.DataFrame:
    selections = []
    for method, path, variable, purpose in (
        ("manta", manta_path, "HC", "MANTA"),
        ("photo_transect", photo_path, "HARD CORAL", "GROUP_LEVEL"),
    ):
        data = pd.read_csv(path)
        required = {"data_type", "domain_category", "reef_name", "report_year",
                    "date", "variable", "purpose", "project_code", "depth",
                    "median", "lower", "upper", "id"}
        if not required.issubset(data.columns):
            raise ValueError(f"Missing LTMP columns in {path}: {required - set(data.columns)}")
        selected = data[(data.data_type == ("manta" if method == "manta" else "photo-transect")) &
                        (data.domain_category == "reef") &
                        (data.reef_name.isin(REEF_ALIASES.values())) &
                        (data.variable == variable) & (data.purpose == purpose) &
                        (data.project_code == "LTMP") & (data.depth == 9) &
                        (data.report_year.between(1985, 2024))].copy()
        selected["method"] = method
        selected["model_reef_name"] = selected.reef_name.map(
            {alias: model for model, alias in REEF_ALIASES.items()})
        selections.append(selected)
    coral = pd.concat(selections, ignore_index=True)
    if coral.groupby(["model_reef_name", "method", "report_year"]).size().ne(1).any():
        raise ValueError("Duplicate LTMP reef/method/report-year coral values")
    if coral[coral.method == "manta"].model_reef_name.nunique() != len(REEF_ALIASES):
        raise ValueError("Manta coral coverage is missing a calibration reef")
    for column in ("median", "lower", "upper"):
        if not np.isfinite(coral[column]).all():
            raise ValueError(f"Nonfinite LTMP {column}")
    if ((coral.lower < 0) | (coral.upper > 1) |
            (coral.lower > coral["median"]) | (coral["median"] > coral.upper)).any():
        raise ValueError("Invalid LTMP cover fraction or source-reported interval")
    coral = coral.rename(columns={"reef_name": "source_reef_name",
                                  "report_year": "year", "date": "survey_date",
                                  "id": "source_id"})
    columns = ["model_reef_name", "source_reef_name", "method", "year",
               "survey_date", "median", "lower", "upper", "source_id"]
    return coral[columns].sort_values(["model_reef_name", "method", "year"])


def observed_cots(ax, observations: pd.DataFrame, reef: str, label: bool) -> None:
    subset = observations[observations.reef_name == reef].sort_values("year")
    if subset.empty:
        raise ValueError(f"Missing COTS observations for {reef}")
    ax.scatter(subset.year, subset.raw_cpue, s=16, color="#111827", alpha=0.4,
               label="Observed COTS/tow" if label else None, zorder=7)
    years = subset.year.to_numpy()
    values = subset.smoothed_3y_cpue.to_numpy()
    starts = np.r_[0, np.flatnonzero(np.diff(years) > 3) + 1]
    ends = np.r_[starts[1:], len(years)]
    for idx, (start, end) in enumerate(zip(starts, ends)):
        ax.plot(years[start:end], values[start:end], color="#111827",
                linewidth=1.35, marker=".", markersize=3.5,
                label="Observed COTS 3-year mean" if label and idx == 0 else None,
                zorder=8)
    ax.axhline(0.22, color="#92400e", linestyle="--", linewidth=1,
               label="0.22 COTS/tow outbreak" if label else None)


def model_lines(ax, trajectories: pd.DataFrame, reef: str, kind: str,
                column: str, label: bool) -> None:
    for treatment, name, color in TREATMENTS[kind]:
        subset = trajectories[(trajectories.reef_name == reef) &
                              (trajectories.treatment == treatment)]
        if subset.empty:
            raise ValueError(f"Missing model trajectory for {reef}: {treatment}")
        summary = subset.groupby("year")[column].agg(
            median="median", low=lambda x: x.quantile(0.1),
            high=lambda x: x.quantile(0.9)).reset_index()
        ax.plot(summary.year, summary["median"], color=color, linewidth=1.7,
                label=name if label else None)
        if subset.seed.nunique() > 1:
            ax.fill_between(summary.year.to_numpy(), summary.low.to_numpy(),
                            summary.high.to_numpy(), color=color, alpha=0.08,
                            linewidth=0)


def observed_coral(ax, coral: pd.DataFrame, reef: str, label: bool) -> None:
    for method, marker, color, name in (
        ("manta", "o", "#111827", "LTMP manta hard coral, 9 m"),
        ("photo_transect", "D", "#a21caf", "LTMP photo-transect hard coral, 9 m"),
    ):
        subset = coral[(coral.model_reef_name == reef) & (coral.method == method)]
        if subset.empty:
            continue
        ax.errorbar(subset.year.to_numpy(), subset["median"].to_numpy(),
                    yerr=[(subset["median"] - subset.lower).to_numpy(),
                          (subset.upper - subset["median"]).to_numpy()],
                    fmt=marker, markersize=4.1, color=color, ecolor=color,
                    alpha=0.9, capsize=1.5, elinewidth=0.75,
                    label=name if label else None, zorder=9)
    if reef == "Eyrie Reef":
        ax.text(0.02, 0.95, "Photo-transect series unavailable",
                transform=ax.transAxes, va="top", fontsize=8, color="#6b7280")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--kind", required=True, choices=TREATMENTS)
    args = parser.parse_args()
    if not args.run_id.replace("_", "").replace("-", "").isalnum():
        raise ValueError("Invalid run ID")
    run_dir = ROOT / "sandbox/calibration/runs" / args.run_id
    paths = {
        "trajectories": run_dir / "trajectories.csv",
        "cots_observations": ROOT / "sandbox/calibration/runs" / GATE_RUN / "observations.csv",
        "reef_manta": ROOT / "sandbox/data/reef_manta.csv",
        "reef_photo_transect": ROOT / "sandbox/data/reef_photo_transect.csv",
        "archived_model": ROOT / "sandbox/data/best_bbo_simulated_trajectories.csv",
    }
    for path in paths.values():
        if not path.is_file():
            raise FileNotFoundError(path)
    trajectories = pd.read_csv(paths["trajectories"])
    if args.kind == "replay":
        trajectories["seed"] = 20260930
    if trajectories.groupby(["reef_name", "treatment", "year", "seed"]).size().ne(1).any():
        raise ValueError("Duplicate model trajectory rows")
    observations = pd.read_csv(paths["cots_observations"])
    coral = load_ltmp_coral(paths["reef_manta"], paths["reef_photo_transect"])
    archive = pd.read_csv(paths["archived_model"])

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(4, 2, figsize=(15, 13), sharex=True,
                             constrained_layout=True)
    for row, reef in enumerate(REEF_ALIASES):
        cots_ax, coral_ax = axes[row]
        observed_cots(cots_ax, observations, reef, row == 0)
        model_lines(cots_ax, trajectories, reef, args.kind, "simulated_cpue", row == 0)
        model_lines(coral_ax, trajectories, reef, args.kind, "coral_cover", row == 0)
        prior = archive[archive.reef_name == reef].sort_values("year")
        coral_ax.plot(prior.year, prior.sim_coral_cover, color="#64748b",
                      linestyle=":", linewidth=1.5,
                      label="Archived model (unweighted)" if row == 0 else None)
        observed_coral(coral_ax, coral, reef, row == 0)
        cots_ax.set_ylabel(f"{reef}\nCOTS/tow")
        coral_ax.set_ylabel("Hard coral cover (fraction)")
        cots_ax.set_ylim(bottom=0)
        coral_ax.set_ylim(0, max(0.6, coral_ax.get_ylim()[1]))
        for ax in (cots_ax, coral_ax):
            ax.set_xlim(1985, 2024)
            ax.grid(axis="y", color="#cbd5e1", linewidth=0.5, alpha=0.6)
        if row == 0:
            cots_ax.set_title("COTS survey and fixed-scale model")
            coral_ax.set_title("LTMP coral surveys and area-weighted model")
        if row == 3:
            cots_ax.set_xlabel("Report/model year")
            coral_ax.set_xlabel("Report/model year")

    handles, labels = [], []
    for ax in (axes[0, 0], axes[0, 1]):
        h, l = ax.get_legend_handles_labels()
        handles.extend(h)
        labels.extend(l)
    unique = dict(zip(labels, handles))
    fig.legend(unique.values(), unique.keys(), loc="upper center", ncol=3,
               bbox_to_anchor=(0.5, 1.06), frameon=False, fontsize=8)
    fig.suptitle("Lizard COTS peaks and LTMP coral cover",
                 y=1.10, fontsize=14)
    fig.text(0.5, -0.015,
             "LTMP points use report year and source-reported intervals; manta and photo-transect "
             "are distinct 9 m methods. Model cover is a polygon-area-weighted whole-reef mean; "
             "archived cover is an unweighted site mean. Coral is not part of the frozen COTS score.",
             ha="center", va="top", fontsize=8, color="#475569", wrap=True)

    plot_dir = run_dir / "plots"
    plot_dir.mkdir(exist_ok=True)
    suffix = "ltmp_manta_photo"
    output = plot_dir / f"cots_coral_{args.kind}_{suffix}.png"
    coral_out = plot_dir / f"coral_observations_{suffix}.csv"
    metadata_out = plot_dir / f"cots_coral_{args.kind}_{suffix}.json"
    if any(path.exists() for path in (output, coral_out, metadata_out)):
        raise FileExistsError("Refusing to overwrite a previous observation revision")
    coral.to_csv(coral_out, index=False)
    fig.savefig(output, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    coverage = {
        reef: {method: int(((coral.model_reef_name == reef) &
                            (coral.method == method)).sum())
               for method in ("manta", "photo_transect")}
        for reef in REEF_ALIASES
    }
    metadata_out.write_text(json.dumps({
        "kind": args.kind,
        "status": "observation_plot_revision_no_score_change",
        "script_sha256": sha256(Path(__file__)),
        "inputs": {name: sha256(path) for name, path in paths.items()},
        "outputs": {output.name: sha256(output), coral_out.name: sha256(coral_out)},
        "reef_aliases": REEF_ALIASES,
        "coverage_1985_2024": coverage,
        "method_filters": {
            "manta": "LTMP reef HC, purpose MANTA, depth 9 m",
            "photo_transect": "LTMP reef HARD CORAL, purpose GROUP_LEVEL, depth 9 m",
        },
        "year_convention": "LTMP report_year aligned to model calendar year; fractional survey date retained in CSV",
        "comparability": "model area-weighted whole-reef cover versus 9 m field methods; no coral scoring",
    }, indent=2) + "\n", encoding="utf-8")
    print(f"Saved {output}")
    print(f"LTMP 1985-2024 coverage: {coverage}")


if __name__ == "__main__":
    main()
