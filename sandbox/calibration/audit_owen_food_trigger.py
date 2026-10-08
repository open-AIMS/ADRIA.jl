"""Screen three pre-registered low-food mortality triggers on archived site states.

This is a no-feedback diagnostic. It does not change the COTS transition or
claim that a hypothetical mortality effect is calibrated or evidence-based.
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
RUNS = ROOT / "sandbox/calibration/runs"
REEFS = ("Lizard Island Reef", "MacGillivray Reef", "North Direction Reef", "Eyrie Reef")
TWO_WAVE = REEFS[:3]
WINDOWS = {"first_wave": (1994, 1997), "trough": (2004, 2009),
           "second_wave": (2012, 2015)}
THRESHOLDS = {"prey_cover": 0.10, "condition": 0.50,
              "prey_per_adult": 0.08}
HYPOTHETICAL_PENALTY = 0.50
COLORS = {"prey_cover": "#429B8A", "condition": "#A85179",
          "prey_per_adult": "#C47A29"}


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
    out = run / "trigger_audit"
    out.mkdir(exist_ok=False)
    sources = {
        "site_diagnostics": run / "site_food_diagnostics.csv",
        "diagnostic_metadata": run / "metadata.toml",
        "candidate": RUNS / "expanded_cotsconn_pilot_seed20260930/best_summary.csv",
    }
    sites = pd.read_csv(sources["site_diagnostics"])
    best = pd.read_csv(sources["candidate"]).sort_values("loss").iloc[0]
    if set(sites.reef_name) != set(REEFS) or sites.year.nunique() != 39:
        raise ValueError("Unexpected diagnostic site or year coverage")
    threshold_legacy = 0.15 * float(best.C_max)
    pre_cover = sites.prey_cover_before_feeding_est.to_numpy(float)
    saturated = ~np.isfinite(pre_cover)
    adult = sites.adult_before_ha.to_numpy(float)
    juvenile = sites.juvenile_before_ha.to_numpy(float)
    adult_floor = np.maximum(adult, 0.01)
    # When condition is saturated, prey cover is at least C_max/2. Verify that
    # this lower bound excludes both cover and per-adult stress in these rows.
    minimum_cover = 0.5 * float(best.C_max)
    if minimum_cover < THRESHOLDS["prey_cover"] or np.any(
        minimum_cover / adult_floor[saturated] < THRESHOLDS["prey_per_adult"]
    ):
        raise ValueError("Saturated condition leaves stress response unidentified")
    legacy_food_survival = np.ones(len(sites))
    low = ~saturated & (pre_cover <= threshold_legacy)
    legacy_food_survival[low] = (1 - float(best.p_tilde)) + float(best.p_tilde) * (
        pre_cover[low] / threshold_legacy
    ) ** 3
    observed_food_survival = sites.food_survival.to_numpy(float)
    checked = np.isfinite(observed_food_survival)
    replay_error = float(np.max(np.abs(
        legacy_food_survival[checked] - observed_food_survival[checked]
    )))
    if replay_error > 1e-8:
        raise ValueError(f"Legacy food-survival reconstruction failed: {replay_error}")

    stresses = {
        "prey_cover": np.clip(1 - np.where(saturated, minimum_cover, pre_cover) /
                              THRESHOLDS["prey_cover"], 0, 1),
        "condition": np.clip(1 - sites.condition_before.to_numpy(float) /
                             THRESHOLDS["condition"], 0, 1),
        "prey_per_adult": np.clip(1 - np.where(saturated, minimum_cover, pre_cover) /
                                  adult_floor / THRESHOLDS["prey_per_adult"], 0, 1),
    }
    annual_rows = []
    for trigger, stress in stresses.items():
        data = sites[["reef_name", "year", "area_weight"]].copy()
        data["adult_mass"] = sites.area_weight.to_numpy(float) * adult
        data["weighted_stress"] = data.adult_mass.to_numpy(float) * stress
        data["hypothetical_extra_adult_deaths_ha"] = (
            sites.area_weight.to_numpy(float) * adult *
            (1 - float(best.m3)) * legacy_food_survival *
            HYPOTHETICAL_PENALTY * stress
        )
        data["hypothetical_extra_juvenile_deaths_ha"] = (
            sites.area_weight.to_numpy(float) * juvenile *
            (1 - float(best.m2)) * legacy_food_survival *
            HYPOTHETICAL_PENALTY * stress
        )
        yearly = data.groupby(["reef_name", "year"], as_index=False).agg(
            adult_mass=("adult_mass", "sum"),
            weighted_stress=("weighted_stress", "sum"),
            hypothetical_extra_adult_deaths_ha=("hypothetical_extra_adult_deaths_ha", "sum"),
            hypothetical_extra_juvenile_deaths_ha=("hypothetical_extra_juvenile_deaths_ha", "sum"),
        )
        yearly["adult_weighted_stress"] = yearly.weighted_stress / yearly.adult_mass
        yearly["trigger"] = trigger
        annual_rows.append(yearly.drop(columns=["adult_mass", "weighted_stress"]))
    annual = pd.concat(annual_rows, ignore_index=True)
    annual.to_csv(out / "trigger_reef_year.csv", index=False)

    window_rows = []
    for trigger in THRESHOLDS:
        for reef in REEFS:
            series = annual[(annual.trigger == trigger) &
                            (annual.reef_name == reef)]
            windows = {name: series[series.year.between(*span)].adult_weighted_stress.mean()
                       for name, span in WINDOWS.items()}
            ratio = windows["trough"] / max(windows["first_wave"],
                                             windows["second_wave"], 1e-12)
            window_rows.append({
                "trigger": trigger,
                "reef_name": reef,
                "first_wave_stress": windows["first_wave"],
                "trough_stress": windows["trough"],
                "second_wave_stress": windows["second_wave"],
                "trough_to_larger_wave_stress": ratio,
                "selectivity_pass": bool(reef in TWO_WAVE and
                                         windows["trough"] >= 0.10 and ratio >= 2.0),
            })
    windows = pd.DataFrame(window_rows)
    windows.to_csv(out / "trigger_window_summary.csv", index=False)

    fig, axes = plt.subplots(4, 1, figsize=(11, 10), sharex=True)
    for ax, reef in zip(axes, REEFS):
        for trigger in THRESHOLDS:
            series = annual[(annual.trigger == trigger) &
                            (annual.reef_name == reef)].sort_values("year")
            ax.plot(series.year, series.adult_weighted_stress, color=COLORS[trigger],
                    label=trigger.replace("_", " "))
        for label, (start, end) in WINDOWS.items():
            ax.axvspan(start, end, color="#253B53", alpha=0.06)
        ax.set_title(reef)
        ax.set_ylabel("Adult-weighted stress")
        ax.set_ylim(0, 1)
        ax.grid(axis="y", alpha=0.18)
    axes[-1].set_xlabel("Calendar year")
    axes[0].legend(ncol=3, frameon=False)
    fig.suptitle("Diagnostic low-food triggers on the unchanged Owen trajectory")
    fig.tight_layout()
    plot = out / "food_trigger_timing.png"
    fig.savefig(plot, dpi=170, bbox_inches="tight")
    plt.close(fig)

    metadata = {
        "status": "no_feedback_precore_screen_not_promoted",
        "stress_definition": "max(0, 1 - input/threshold)",
        "thresholds_diagnostic_not_evidence_bounded": THRESHOLDS,
        "hypothetical_extra_survival_penalty_strength": HYPOTHETICAL_PENALTY,
        "adult_density_floor_for_ratio_cots_ha": 0.01,
        "window_years": {k: list(v) for k, v in WINDOWS.items()},
        "gate": "trough stress >= 0.10 and >= 2x both wave stresses on all three two-wave reefs",
        "food_survival_replay_error": replay_error,
        "source_sha256": {key: digest(path) for key, path in sources.items()},
        "script_sha256": digest(Path(__file__)),
        "output_sha256": {path.name: digest(path) for path in
                          (out / "trigger_reef_year.csv", out / "trigger_window_summary.csv", plot)},
    }
    (out / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n",
                                       encoding="utf-8")
    for trigger in THRESHOLDS:
        subset = windows[(windows.trigger == trigger) &
                         (windows.reef_name.isin(TWO_WAVE))]
        print(trigger, "gate:", bool(subset.selectivity_pass.all()),
              "trough/wave ratios:",
              dict(zip(subset.reef_name, subset.trough_to_larger_wave_stress.round(3))))
    print("Saved", out)


if __name__ == "__main__":
    main()
