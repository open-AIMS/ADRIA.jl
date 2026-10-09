"""Plot the exact-control 1991 cohort readout without changing model dynamics."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=Path)
    args = parser.parse_args()
    run_dir = args.run_dir.resolve(strict=True)
    input_files = [run_dir / name for name in ("cohorts.csv", "flows.csv", "trajectories.csv")]
    cohorts, flows, trajectories = (pd.read_csv(path) for path in input_files)
    treatment = "cohort_full_control"
    reef = "Lizard Island Reef"
    cohorts = cohorts.query("treatment == @treatment and reef_name == @reef and 1996 <= year <= 2015")
    flows = flows.query("treatment == @treatment and reef_name == @reef and 1996 <= year <= 2015")
    trajectories = trajectories.query("treatment == @treatment and reef_name == @reef and 1996 <= year <= 2015")
    if not (len(cohorts) == len(flows) == len(trajectories) == 20):
        raise ValueError("Unexpected cohort audit coverage")
    if not (cohorts.year.to_numpy() == flows.year.to_numpy()).all():
        raise ValueError("Cohort and flow calendars differ")
    if not (cohorts.year.to_numpy() == trajectories.year.to_numpy()).all():
        raise ValueError("Cohort and trajectory calendars differ")

    scale = float(trajectories.simulated_cpue.iloc[0] / trajectories.adults_ha.iloc[0])
    x = cohorts.year.to_numpy()
    fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True, constrained_layout=True)
    fig.suptitle("Lizard 1991-control cohort progression, 1996–2015")

    axes[0].plot(x, cohorts.adults_ha * scale, color="#1f2937", linewidth=2.2, label="Adults, total")
    axes[0].plot(x, cohorts.adult_carryover_ha * scale, color="#2563eb", label="Adult carryover")
    axes[0].plot(x, cohorts.maturation_ha * scale, color="#d97706", label="Newly matured")
    axes[0].set_ylabel("Adults (COTS/tow)")
    axes[0].legend(loc="upper right", ncol=3, frameon=False)

    axes[1].plot(x, cohorts.recruits_ha, color="#059669", label="Settled/recruits N1")
    axes[1].plot(x, cohorts.subadults_ha, color="#9333ea", label="Immature N2")
    axes[1].plot(x, flows.external_immigration, color="#6b7280", linestyle="--", label="External settlement")
    axes[1].set_ylabel("COTS/ha or COTS/ha/year")
    axes[1].legend(loc="upper right", ncol=3, frameon=False)

    axes[2].plot(x, cohorts.adult_food_factor, color="#2563eb", label="Adult food-survival factor")
    axes[2].plot(x, cohorts.subadult_food_factor, color="#9333ea", linestyle="--", label="N2 food-survival factor")
    axes[2].set_ylabel("Realized food-survival factor")
    axes[2].set_ylim(0, 1.08)
    coral_axis = axes[2].twinx()
    coral_axis.plot(x, trajectories.coral_cover, color="#16a34a", alpha=0.85, label="Coral after COTS grazing")
    coral_axis.axhline(0.4980659994966138 * 0.15, color="#16a34a", alpha=0.55,
                       linestyle=":", label="Starvation threshold (pre-grazing)")
    coral_axis.set_ylabel("Coral cover fraction")
    coral_axis.set_ylim(0, 0.16)
    handles, labels = axes[2].get_legend_handles_labels()
    coral_handles, coral_labels = coral_axis.get_legend_handles_labels()
    axes[2].legend(handles + coral_handles, labels + coral_labels, loc="upper right", ncol=2, frameon=False)

    for axis in axes:
        axis.grid(axis="y", alpha=0.2)
        axis.set_xlim(1996, 2015)
    axes[-1].set_xticks(range(1996, 2016, 2))
    axes[-1].set_xlabel("Calendar year")
    output = run_dir / "cohort_floor_diagnostics.png"
    provenance = run_dir / "cohort_floor_plot_provenance.json"
    if output.exists() or provenance.exists():
        raise FileExistsError("Refusing to overwrite an existing cohort plot")
    fig.savefig(output, dpi=180)
    plt.close(fig)
    provenance.write_text(json.dumps({
        "source_sha256": {path.name: sha256(path) for path in input_files},
        "script_sha256": sha256(Path(__file__).resolve()),
        "output_sha256": sha256(output),
        "reef": reef,
        "year_range": [1996, 2015],
        "cpue_scale_tow_per_adult_ha": scale,
        "coral_timing_note": "Plotted coral is post-grazing; food survival used pre-grazing coral.",
    }, indent=2) + "\n", encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()
