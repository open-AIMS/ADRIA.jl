"""Audit the predeclared inter-wave crash screen and plot the oracle comparisons."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from audit_owen_peak_shape import peak_years, shape

ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / "sandbox/calibration/runs"
GATE = RUNS / "20261001T113401_lizard_domain_gate"
OBSERVED_SHAPE = (
    RUNS / "20261001T162000_owen_embedded_mortality_screen"
    / "shape_audit/interpeak_troughs.csv"
)
REEFS = ["Lizard Island Reef", "MacGillivray Reef", "North Direction Reef", "Eyrie Reef"]
TWO_WAVE_REEFS = set(REEFS[:3])
TREATMENTS = [
    "joint4_reference", "external_off", "internal_off", "both_off", "both_off_m3_upper"
]
COLORS = ["#253B53", "#D88B2E", "#429B8A", "#B25279", "#7E65B1"]


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
    out = run / "crash_audit"
    out.mkdir(exist_ok=False)
    sources = {
        "scores": run / "scores.csv",
        "trajectories": run / "trajectories.csv",
        "flows": run / "flows.csv",
        "failures": run / "failures.csv",
        "observations": GATE / "observations.csv",
        "observed_troughs": OBSERVED_SHAPE,
    }
    scores = pd.read_csv(sources["scores"])
    scores = scores[scores.observation_treatment == "smoothed_3y"]
    trajectories = pd.read_csv(sources["trajectories"])
    flows = pd.read_csv(sources["flows"])
    failures = pd.read_csv(sources["failures"])
    observations = pd.read_csv(sources["observations"])
    observed_shape = pd.read_csv(sources["observed_troughs"])
    observed_shape = observed_shape[observed_shape.treatment == "observed"]
    if len(failures) or len(scores) != len(TREATMENTS) * len(REEFS):
        raise ValueError("Failures or incomplete treatment scores")
    if len(trajectories) != len(TREATMENTS) * len(REEFS) * 40:
        raise ValueError("Incomplete annual trajectory output")
    if len(flows) != len(trajectories):
        raise ValueError("Incomplete flow output")

    rows: list[dict[str, object]] = []
    for score in scores.itertuples():
        data = trajectories[(trajectories.treatment == score.treatment) &
                            (trajectories.reef_name == score.reef_name)]
        series = data.set_index("year").simulated_cpue.sort_index()
        smooth = series.rolling(3, center=True, min_periods=1).mean()
        peaks = peak_years(score.sim_peak_years)
        trough = shape(smooth, *peaks) if len(peaks) == 2 else None
        rows.append({
            "treatment": score.treatment,
            "reef_name": score.reef_name,
            "sim_peak_count": int(score.sim_peak_count),
            "matched_peaks": int(score.matched_peaks),
            "frozen_loss": float(score.loss),
            "first_peak_year": trough["first_peak_year"] if trough else np.nan,
            "second_peak_year": trough["second_peak_year"] if trough else np.nan,
            "first_peak_cpue": trough["first_peak_cpue"] if trough else np.nan,
            "second_peak_cpue": trough["second_peak_cpue"] if trough else np.nan,
            "trough_year": trough["trough_year"] if trough else np.nan,
            "trough_cpue": trough["trough_cpue"] if trough else np.nan,
            "trough_to_smaller_peak": trough["trough_to_smaller_peak"] if trough else np.nan,
            "trough_screen_pass": bool(
                score.reef_name in TWO_WAVE_REEFS and trough and
                trough["trough_to_smaller_peak"] <= 0.10
            ),
        })
    result = pd.DataFrame(rows).sort_values(["treatment", "reef_name"])
    result.to_csv(out / "peak_shape.csv", index=False)

    gap = flows[flows.year.between(2000, 2010)]
    summary = gap.groupby(["treatment", "reef_name"], as_index=False).agg(
        mean_local_fecundity=("local_fecundity", "mean"),
        mean_local_retention=("local_retention", "mean"),
        mean_internal_immigration=("internal_immigration", "mean"),
        mean_external_pelagic=("external_pelagic", "mean"),
        mean_maturation=("maturation", "mean"),
        mean_settled_recruits=("settled_recruits", "mean"),
    )
    summary.to_csv(out / "gap_flows.csv", index=False)
    for treatment in TREATMENTS:
        subset = result[(result.treatment == treatment) &
                        (result.reef_name.isin(TWO_WAVE_REEFS))]
        matches = result[result.treatment == treatment].matched_peaks.sum()
        passed = bool(len(subset) == 3 and subset.trough_screen_pass.all() and matches >= 7)
        print(f"{treatment}: matched={matches}/7, trough_gate={passed}, "
              f"ratios={dict(zip(subset.reef_name, subset.trough_to_smaller_peak.round(3)))}")

    fig, axes = plt.subplots(2, 2, figsize=(13.5, 8.4), sharex=True)
    for ax, reef in zip(axes.flat, REEFS):
        obs = observations[observations.reef_name == reef]
        ax.scatter(obs.year, obs.raw_cpue, color="#B0AAA0", s=10, alpha=0.55, zorder=2)
        ax.plot(obs.year, obs.smoothed_3y_cpue, color="#111111", linewidth=2.0,
                label="LTMP 3-year mean", zorder=5)
        for treatment, color in zip(TREATMENTS, COLORS):
            data = trajectories[(trajectories.treatment == treatment) &
                                (trajectories.reef_name == reef)].sort_values("year")
            smoothed_model = data.simulated_cpue.rolling(3, center=True,
                                                        min_periods=1).mean()
            ax.plot(data.year, smoothed_model, color=color, linewidth=1.4,
                    label=treatment.replace("_", " "))
        ax.axvspan(2000, 2010, color="#999999", alpha=0.10, linewidth=0)
        ax.set_title(reef)
        ax.set_xlim(1989, 2024)
        ax.set_ylim(bottom=0)
        ax.grid(axis="y", alpha=0.18)
        ax.set_ylabel("COTS per manta tow")
        observed = observed_shape[observed_shape.reef_name == reef]
        if len(observed):
            ax.text(0.02, 0.96, f"Observed trough / peak: "
                    f"{observed.iloc[0].trough_to_smaller_peak:.3f}",
                    transform=ax.transAxes, va="top", fontsize=9,
                    bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.75})
    for ax in axes[-1]:
        ax.set_xlabel("Calendar year")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, -0.015))
    fig.suptitle("Inter-wave COTS crash: source-blackout oracle screen", fontsize=14)
    fig.tight_layout(rect=(0, 0.075, 1, 0.96))
    plot = out / "cots_recruitment_blackout.png"
    fig.savefig(plot, dpi=180, bbox_inches="tight")
    plt.close(fig)

    metadata = {
        "status": "diagnostic_only_not_frozen_score_or_historical_reconstruction",
        "trough_definition": "3-year-smoothed minimum strictly between two detected peaks / smaller peak",
        "screening_cutoff": 0.10,
        "two_wave_reefs": sorted(TWO_WAVE_REEFS),
        "blackout_calendar_years_inclusive": [2000, 2010],
        "source_paths": {key: str(path) for key, path in sources.items()},
        "source_sha256": {key: digest(path) for key, path in sources.items()},
        "script_sha256": digest(Path(__file__)),
        "output_sha256": {path.name: digest(path) for path in
                          (out / "peak_shape.csv", out / "gap_flows.csv", plot)},
    }
    (out / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n",
                                       encoding="utf-8")


if __name__ == "__main__":
    main()
