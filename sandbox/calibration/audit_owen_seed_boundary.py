"""Score the preregistered Owen spatial-seed x outside-supply diagnostic."""

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

from audit_owen_peak_shape import peak_years, shape
from plot_lizard_ltmp_coral import load_ltmp_coral

ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / "sandbox/calibration/runs"
REEFS = ("Lizard Island Reef", "MacGillivray Reef", "North Direction Reef", "Eyrie Reef")
TWO_WAVE = set(REEFS[:3])
TREATMENTS = ("reference_full", "reference_10pct", "reference_zero",
              "quarter_full", "quarter_10pct", "quarter_zero")
BOUNDARY_COLORS = {"full": "#335A81", "10pct": "#BD7824", "zero": "#478877"}


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
    run = RUNS / args.run_id
    out = run / "seed_boundary_audit"
    out.mkdir(exist_ok=False)
    metadata_path = run / "metadata.json"
    if not metadata_path.exists():
        metadata_path = run / "metadata.toml"
    paths = {
        "scores": run / "scores.csv",
        "trajectories": run / "trajectories.csv",
        "flows": run / "flows.csv",
        "initialization": run / "initialization.csv",
        "failures": run / "failures.csv",
        "observations": RUNS / "20261001T113401_lizard_domain_gate/observations.csv",
        "manta": ROOT / "sandbox/data/reef_manta.csv",
        "photo": ROOT / "sandbox/data/reef_photo_transect.csv",
        "run_metadata": metadata_path,
    }
    scores = pd.read_csv(paths["scores"])
    scores = scores[scores.observation_treatment == "smoothed_3y"]
    trajectories = pd.read_csv(paths["trajectories"])
    flows = pd.read_csv(paths["flows"])
    initial = pd.read_csv(paths["initialization"])
    failures = pd.read_csv(paths["failures"])
    observations = pd.read_csv(paths["observations"])
    coral_obs = load_ltmp_coral(paths["manta"], paths["photo"])
    if len(failures) or len(scores) != len(TREATMENTS) * len(REEFS):
        raise ValueError("Failed or incomplete screen")
    if len(trajectories) != len(TREATMENTS) * len(REEFS) * 40:
        raise ValueError("Incomplete trajectories")
    if len(flows) != len(trajectories) or len(initial) != len(TREATMENTS) * len(REEFS):
        raise ValueError("Incomplete flows or initialization")

    rows = []
    for score in scores.itertuples():
        data = trajectories[(trajectories.treatment == score.treatment) &
                            (trajectories.reef_name == score.reef_name)].sort_values("year")
        series = data.set_index("year").simulated_cpue.rolling(
            3, center=True, min_periods=1).mean()
        peaks = peak_years(score.sim_peak_years)
        trough = shape(series, *peaks) if len(peaks) == 2 else None
        seed = initial[(initial.treatment == score.treatment) &
                       (initial.reef_name == score.reef_name)]
        if len(seed) != 1:
            raise ValueError("Missing initial state")
        flows_reef = flows[(flows.treatment == score.treatment) &
                           (flows.reef_name == score.reef_name)]
        if score.treatment.endswith("zero") and flows_reef.external_pelagic.abs().max() > 1e-12:
            raise ValueError("Zero boundary leaked outside supply")
        gap = flows_reef[flows_reef.year.between(2000, 2010)]
        rows.append({
            "treatment": score.treatment, "reef_name": score.reef_name,
            "initial_adults_ha": float(seed.adult_ha.iloc[0]),
            "initial_juveniles_ha": float(seed.juvenile_ha.iloc[0]),
            "1989_model_cpue": float(data.loc[data.year == 1989, "simulated_cpue"].iloc[0]),
            "sim_peak_count": int(score.sim_peak_count),
            "sim_peak_years": score.sim_peak_years,
            "matched_peaks": int(score.matched_peaks),
            "frozen_loss": float(score.loss),
            "first_peak_cpue": trough["first_peak_cpue"] if trough else np.nan,
            "second_peak_cpue": trough["second_peak_cpue"] if trough else np.nan,
            "trough_year": trough["trough_year"] if trough else np.nan,
            "trough_cpue": trough["trough_cpue"] if trough else np.nan,
            "trough_to_smaller_peak": trough["trough_to_smaller_peak"] if trough else np.nan,
            "trough_screen_pass": bool(score.reef_name in TWO_WAVE and trough and
                                       trough["trough_to_smaller_peak"] <= 0.10),
            "gap_local_fecundity": float(gap.local_fecundity.mean()),
            "gap_internal_immigration": float(gap.internal_immigration.mean()),
            "gap_external_pelagic": float(gap.external_pelagic.mean()),
            "model_coral_2015": float(data.loc[data.year == 2015, "coral_cover"].iloc[0]),
        })
    summary = pd.DataFrame(rows).sort_values(["treatment", "reef_name"])
    summary.to_csv(out / "peak_shape_flux_and_coral.csv", index=False)

    fig, axes = plt.subplots(4, 2, figsize=(14, 14), sharex="col")
    for row, reef in enumerate(REEFS):
        cots_ax, coral_ax = axes[row]
        obs = observations[observations.reef_name == reef].sort_values("year")
        cots_ax.scatter(obs.year, obs.raw_cpue, s=10, color="#999999", alpha=0.5,
                        label="LTMP tow points" if row == 0 else None)
        cots_ax.plot(obs.year, obs.smoothed_3y_cpue, color="black", linewidth=2,
                     label="LTMP 3-year mean" if row == 0 else None)
        for treatment in TREATMENTS:
            data = trajectories[(trajectories.reef_name == reef) &
                                (trajectories.treatment == treatment)].sort_values("year")
            boundary = treatment.split("_", 1)[1]
            style = "-" if treatment.startswith("reference") else "--"
            label = treatment.replace("_", " ") if row == 0 else None
            color = BOUNDARY_COLORS[boundary]
            cots_ax.plot(data.year, data.simulated_cpue.rolling(
                3, center=True, min_periods=1).mean(), color=color,
                linestyle=style, linewidth=1.55, label=label)
            coral_ax.plot(data.year, data.coral_cover, color=color,
                          linestyle=style, linewidth=1.55)
        manta = coral_obs[(coral_obs.model_reef_name == reef) &
                          (coral_obs.method == "manta")]
        photo = coral_obs[(coral_obs.model_reef_name == reef) &
                          (coral_obs.method == "photo_transect")]
        coral_ax.scatter(manta.year, manta["median"], s=17, color="black",
                         label="LTMP manta 9 m" if row == 0 else None)
        if len(photo):
            coral_ax.scatter(photo.year, photo["median"], s=15, marker="D",
                             color="#8C43A1", label="LTMP photo 9 m" if row == 0 else None)
        cots_ax.set_title(reef)
        cots_ax.set_ylabel("COTS / tow")
        coral_ax.set_ylabel("Hard coral fraction")
        for ax in (cots_ax, coral_ax):
            ax.set_xlim(1989, 2024)
            ax.set_ylim(bottom=0)
            ax.grid(axis="y", alpha=0.16)
    axes[-1, 0].set_xlabel("Calendar year")
    axes[-1, 1].set_xlabel("Calendar year")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    coral_handles, coral_labels = axes[0, 1].get_legend_handles_labels()
    fig.legend(handles + coral_handles, labels + coral_labels,
               loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 0.0))
    fig.suptitle("Owen seed and outside supply: COTS waves and coral", fontsize=14)
    fig.tight_layout(rect=(0, 0.085, 1, 0.975))
    plot = out / "cots_coral_seed_boundary.png"
    fig.savefig(plot, dpi=170, bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(4, 1, figsize=(10, 11), sharex=True)
    for ax, reef in zip(axes, REEFS):
        obs = observations[observations.reef_name == reef].sort_values("year")
        ax.scatter(obs.year, obs.raw_cpue, s=11, color="#999999", alpha=0.5)
        ax.plot(obs.year, obs.smoothed_3y_cpue, color="black", linewidth=2,
                label="LTMP 3-year mean")
        for treatment in TREATMENTS:
            data = trajectories[(trajectories.reef_name == reef) &
                                (trajectories.treatment == treatment)].sort_values("year")
            ax.plot(data.year, data.simulated_cpue.rolling(3, center=True,
                    min_periods=1).mean(), color=BOUNDARY_COLORS[treatment.split("_", 1)[1]],
                    linestyle="-" if treatment.startswith("reference") else "--",
                    linewidth=1.55, label=treatment.replace("_", " "))
        ax.set_title(reef)
        ax.set_ylabel("COTS / tow")
        ax.set_xlim(1989, 2002)
        ax.set_ylim(bottom=0)
        ax.grid(axis="y", alpha=0.16)
    axes[-1].set_xlabel("Calendar year")
    axes[0].legend(loc="upper left", ncol=2, frameon=False, fontsize=8)
    fig.suptitle("First-wave zoom: 1989–2002", fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    zoom = out / "first_wave_zoom.png"
    fig.savefig(zoom, dpi=170, bbox_inches="tight")
    plt.close(fig)

    meta = {
        "status": "diagnostic_only_not_promotion",
        "trough_definition": "3-year-smoothed minimum / smaller adjacent peak",
        "screen_cutoff": 0.10,
        "coral_support_caution": "model area-weighted whole reef; LTMP manta/photo at 9 m",
        "source_sha256": {name: digest(path) for name, path in paths.items()},
        "script_sha256": digest(Path(__file__)),
        "output_sha256": {path.name: digest(path) for path in
                          (out / "peak_shape_flux_and_coral.csv", plot, zoom)},
    }
    (out / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    for treatment in TREATMENTS:
        subset = summary[summary.treatment == treatment]
        two_wave = subset[subset.reef_name.isin(TWO_WAVE)]
        print(treatment, "matches", int(subset.matched_peaks.sum()),
              "trough_pass", bool(two_wave.trough_screen_pass.all()),
              "ratios", dict(zip(two_wave.reef_name,
                                  two_wave.trough_to_smaller_peak.round(3))))
    print("Saved", out)


if __name__ == "__main__":
    main()
