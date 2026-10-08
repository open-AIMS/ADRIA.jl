"""Audit paired external condition-proxy COTS runs and plot wave/coral checks."""

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
TREATMENTS = ("joint4_reference", "external_condition_proxy")
COLORS = ("#253B53", "#C47A29")


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
    out = run / "source_condition_audit"
    out.mkdir(exist_ok=False)
    sources = {
        "scores": run / "scores.csv",
        "trajectories": run / "trajectories.csv",
        "flows": run / "flows.csv",
        "failures": run / "failures.csv",
        "observations": RUNS / "20261001T113401_lizard_domain_gate/observations.csv",
        "manta": ROOT / "sandbox/data/reef_manta.csv",
        "photo": ROOT / "sandbox/data/reef_photo_transect.csv",
        "run_metadata": run / "metadata.toml",
    }
    scores = pd.read_csv(sources["scores"])
    scores = scores[scores.observation_treatment == "smoothed_3y"]
    trajectories = pd.read_csv(sources["trajectories"])
    flows = pd.read_csv(sources["flows"])
    failures = pd.read_csv(sources["failures"])
    observations = pd.read_csv(sources["observations"])
    coral_obs = load_ltmp_coral(sources["manta"], sources["photo"])
    if len(failures) or len(scores) != len(TREATMENTS) * len(REEFS):
        raise ValueError("Failed or incomplete paired screen")
    if len(trajectories) != len(TREATMENTS) * len(REEFS) * 40:
        raise ValueError("Incomplete annual trajectories")
    if len(flows) != len(trajectories):
        raise ValueError("Incomplete annual flows")

    rows = []
    for score in scores.itertuples():
        data = trajectories[(trajectories.treatment == score.treatment) &
                            (trajectories.reef_name == score.reef_name)].sort_values("year")
        series = data.set_index("year").simulated_cpue.rolling(
            3, center=True, min_periods=1).mean()
        peaks = peak_years(score.sim_peak_years)
        trough = shape(series, *peaks) if len(peaks) == 2 else None
        year_2015 = data[data.year == 2015]
        if len(year_2015) != 1:
            raise ValueError("Missing 2015 model coral")
        manta_2015 = coral_obs[(coral_obs.model_reef_name == score.reef_name) &
                               (coral_obs.method == "manta") & (coral_obs.year == 2015)]
        rows.append({
            "treatment": score.treatment,
            "reef_name": score.reef_name,
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
            "model_coral_2015": float(year_2015.coral_cover.iloc[0]),
            "manta_9m_median_coral_2015": (
                float(manta_2015["median"].iloc[0]) if len(manta_2015) == 1 else np.nan),
        })
    result = pd.DataFrame(rows).sort_values(["treatment", "reef_name"])
    result.to_csv(out / "peak_shape_and_coral.csv", index=False)
    gap = flows[flows.year.between(2000, 2010)]
    gap_summary = gap.groupby(["treatment", "reef_name"], as_index=False).agg(
        mean_local_fecundity=("local_fecundity", "mean"),
        mean_local_retention=("local_retention", "mean"),
        mean_internal_immigration=("internal_immigration", "mean"),
        mean_external_pelagic=("external_pelagic", "mean"),
        mean_maturation=("maturation", "mean"),
        mean_settled_recruits=("settled_recruits", "mean"),
    )
    gap_summary.to_csv(out / "gap_flows.csv", index=False)

    fig, axes = plt.subplots(4, 2, figsize=(13, 13), sharex="col")
    for row, reef in enumerate(REEFS):
        cots_ax, coral_ax = axes[row]
        obs = observations[observations.reef_name == reef].sort_values("year")
        cots_ax.scatter(obs.year, obs.raw_cpue, s=10, color="#969696", alpha=0.5,
                        label="LTMP tow points" if row == 0 else None)
        cots_ax.plot(obs.year, obs.smoothed_3y_cpue, color="black", linewidth=2,
                     label="LTMP 3-year mean" if row == 0 else None)
        for treatment, color in zip(TREATMENTS, COLORS):
            data = trajectories[(trajectories.reef_name == reef) &
                                (trajectories.treatment == treatment)].sort_values("year")
            cots_ax.plot(data.year, data.simulated_cpue.rolling(
                3, center=True, min_periods=1).mean(), color=color, linewidth=1.6,
                label=treatment.replace("_", " ") if row == 0 else None)
            coral_ax.plot(data.year, data.coral_cover, color=color, linewidth=1.6)
        manta = coral_obs[(coral_obs.model_reef_name == reef) &
                          (coral_obs.method == "manta")]
        photo = coral_obs[(coral_obs.model_reef_name == reef) &
                          (coral_obs.method == "photo_transect")]
        coral_ax.scatter(manta.year, manta["median"], s=16, color="black",
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
               loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.0))
    fig.suptitle("Owen source-condition proxy: COTS waves and coral", fontsize=14)
    fig.tight_layout(rect=(0, 0.06, 1, 0.975))
    plot = out / "cots_coral_source_condition_proxy.png"
    fig.savefig(plot, dpi=170, bbox_inches="tight")
    plt.close(fig)
    metadata = {
        "status": "diagnostic_only_not_frozen_score_or_promotion",
        "trough_definition": "3-year-smoothed minimum strictly between two detected peaks / smaller peak",
        "screening_cutoff": 0.10,
        "coral_support_caution": "model area-weighted whole reef; LTMP manta/photo at 9 m",
        "source_sha256": {key: digest(path) for key, path in sources.items()},
        "script_sha256": digest(Path(__file__)),
        "output_sha256": {path.name: digest(path) for path in
                          (out / "peak_shape_and_coral.csv", out / "gap_flows.csv", plot)},
    }
    (out / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n",
                                       encoding="utf-8")
    for treatment in TREATMENTS:
        subset = result[(result.treatment == treatment) &
                        (result.reef_name.isin(TWO_WAVE))]
        matches = int(result[result.treatment == treatment].matched_peaks.sum())
        passed = bool(len(subset) == 3 and subset.trough_screen_pass.all() and matches >= 7)
        print(f"{treatment}: matched={matches}/7, trough_gate={passed}, "
              f"ratios={dict(zip(subset.reef_name, subset.trough_to_smaller_peak.round(3)))}")
    print("Saved", out)


if __name__ == "__main__":
    main()
