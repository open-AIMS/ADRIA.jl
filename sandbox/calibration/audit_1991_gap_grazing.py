"""Audit the default-off 2004-2012 grazing-bypass diagnostic."""

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
REEFS = {
    "Lizard Island Reef": "Lizard Isles",
    "MacGillivray Reef": "Macgillivray Reef",
    "North Direction Reef": "North Direction Island",
    "Eyrie Reef": "Eyrie Reef",
}
CONTROL = "inherited_1991"
TREATMENT = "inherited_grazing_off_2004_2012"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main(run: Path) -> None:
    output = run / "gap_grazing_audit.csv"
    plot = run / "cots_coral_gap_grazing.png"
    metadata = run / "gap_grazing_audit.json"
    if any(path.exists() for path in (output, plot, metadata)):
        raise FileExistsError("Refusing to overwrite existing gap-grazing audit")
    trajectories = pd.read_csv(run / "trajectories.csv")
    peaks = pd.read_csv(run / "peak_trough_audit.csv")
    cots = pd.read_csv(ROOT / "sandbox/data/reef_cots.csv")
    manta = pd.read_csv(ROOT / "sandbox/data/reef_manta.csv")
    rows = []
    fig, axes = plt.subplots(4, 2, figsize=(12, 12), sharex=True)
    for index, (reef, observed) in enumerate(REEFS.items()):
        subset = trajectories[trajectories.reef_name == reef]
        control = subset[subset.treatment == CONTROL].sort_values("year").set_index("year")
        treatment = subset[subset.treatment == TREATMENT].sort_values("year").set_index("year")
        if len(control) != 34 or len(treatment) != 34:
            raise ValueError(f"Incomplete trajectory: {reef}")
        early = slice(1991, 2003)
        prehistory_difference = max(
            np.max(np.abs(control.loc[early, "adults_ha"] - treatment.loc[early, "adults_ha"])),
            np.max(np.abs(control.loc[early, "coral_cover"] - treatment.loc[early, "coral_cover"])),
        )
        if prehistory_difference > 1e-12:
            raise ValueError(f"Diagnostic changed pre-2004 history at {reef}")
        cpeak = peaks[(peaks.reef_name == reef) & (peaks.treatment == CONTROL)].iloc[0]
        tpeak = peaks[(peaks.reef_name == reef) & (peaks.treatment == TREATMENT)].iloc[0]
        rows.append({
            "reef_name": reef,
            "pre_2004_max_abs_adult_or_coral_difference": prehistory_difference,
            "control_mean_coral_2004_2012": control.loc[2004:2012, "coral_cover"].mean(),
            "grazing_off_mean_coral_2004_2012": treatment.loc[2004:2012, "coral_cover"].mean(),
            "control_coral_2012": control.loc[2012, "coral_cover"],
            "grazing_off_coral_2012": treatment.loc[2012, "coral_cover"],
            "control_mean_adults_ha_2004_2012": control.loc[2004:2012, "adults_ha"].mean(),
            "grazing_off_mean_adults_ha_2004_2012": treatment.loc[2004:2012, "adults_ha"].mean(),
            "control_second_peak_cpue": cpeak.second_peak_cpue,
            "grazing_off_second_peak_cpue": tpeak.second_peak_cpue,
            "control_trough_ratio": cpeak.trough_to_smaller_peak,
            "grazing_off_trough_ratio": tpeak.trough_to_smaller_peak,
        })
        ax_cots, ax_coral = axes[index]
        for name, series, color, style in (
            ("Control", control, "#475569", "--"),
            ("Grazing off 2004–12", treatment, "#dc2626", "-"),
        ):
            ax_cots.plot(series.index[1:], series.simulated_cpue.iloc[1:],
                         color=color, linestyle=style, linewidth=1.5,
                         label=name if index == 0 else None)
            ax_coral.plot(series.index, series.coral_cover,
                          color=color, linestyle=style, linewidth=1.5)
        obs_cots = cots[(cots.reef_name == observed) & cots.year.between(1992, 2024)]
        ax_cots.scatter(obs_cots.year, obs_cots.cotsptow, s=17, color="#111827",
                        alpha=0.65, label="LTMP COTS/tow" if index == 0 else None)
        obs_coral = manta[(manta.reef_name == observed) &
                          (manta.data_type == "manta") &
                          (manta.domain_category == "reef") &
                          (manta.variable == "HC") &
                          (manta.purpose == "MANTA") &
                          (manta.project_code == "LTMP") &
                          (manta.depth == 9) &
                          manta.report_year.between(1992, 2024)]
        ax_coral.scatter(obs_coral.report_year, obs_coral["median"],
                         s=17, color="#111827", alpha=0.6,
                         label="LTMP 9 m coral" if index == 0 else None)
        for ax in (ax_cots, ax_coral):
            ax.axvspan(2004, 2012, color="#fef3c7", alpha=0.35)
            ax.set_xlim(1991, 2024)
            ax.set_ylim(bottom=0)
            ax.grid(alpha=0.15)
        ax_cots.set_ylabel(f"{reef}\nCOTS/tow")
        ax_coral.set_ylabel("Hard-coral fraction")
    axes[0, 0].set_title("COTS: fixed observation scale")
    axes[0, 1].set_title("Coral: model whole reef vs 9 m manta")
    axes[-1, 0].set_xlabel("Year")
    axes[-1, 1].set_xlabel("Year")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    h2, l2 = axes[0, 1].get_legend_handles_labels()
    fig.legend(handles + h2, labels + l2, loc="lower center", ncol=4,
               fontsize=8, bbox_to_anchor=(0.5, -0.01))
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(plot, dpi=170, bbox_inches="tight")
    plt.close(fig)
    pd.DataFrame(rows).to_csv(output, index=False)
    metadata.write_text(json.dumps({
        "status": "diagnostic_not_promoted",
        "interpretation": "Bypass updates COTS demography on a coral copy; it tests direct grazing feedback, not an ecological mechanism or independent source history.",
        "inputs": {name: sha256(path) for name, path in {
            "trajectories": run / "trajectories.csv",
            "peak_trough_audit": run / "peak_trough_audit.csv",
            "cots": ROOT / "sandbox/data/reef_cots.csv",
            "manta": ROOT / "sandbox/data/reef_manta.csv",
            "script": Path(__file__),
        }.items()},
        "outputs": {output.name: sha256(output), plot.name: sha256(plot)},
    }, indent=2) + "\n")
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    main(parser.parse_args().run)
