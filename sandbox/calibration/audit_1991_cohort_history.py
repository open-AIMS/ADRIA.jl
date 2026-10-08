"""Audit a bounded 1991 immature-cohort screen and coral-only counterfactual."""

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
PREVIOUS = ROOT / "sandbox/calibration/runs/20261007T199105_owen_1991_initialization"
REEFS = {
    "Lizard Island Reef": "Lizard Isles",
    "MacGillivray Reef": "Macgillivray Reef",
    "North Direction Reef": "North Direction Island",
    "Eyrie Reef": "Eyrie Reef",
}
TREATMENTS = {
    "inherited_1991": ("Inherited 1991", "#475569", "--"),
    "idw_alpha015_old_coral": ("IDW adults only", "#ca8a04", "-"),
    "idw_alpha015_cohort_half": ("IDW adults + half 6:3:1", "#0d9488", "-"),
    "idw_alpha015_cohort_full": ("IDW adults + full 6:3:1", "#dc2626", "-"),
    "coral_only_no_cots": ("No COTS (coral only)", "#2563eb", ":"),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main(run: Path) -> None:
    output = run / "cohort_history_audit.csv"
    plot = run / "cots_coral_cohort_history.png"
    metadata = run / "cohort_history_audit.json"
    if any(p.exists() for p in (output, plot, metadata)):
        raise FileExistsError("Refusing to overwrite existing cohort audit")
    current = pd.read_csv(run / "trajectories.csv")
    previous = pd.read_csv(PREVIOUS / "trajectories.csv")
    scores = pd.read_csv(run / "scores.csv")
    flows = pd.read_csv(run / "flows.csv")
    previous_flows = pd.read_csv(PREVIOUS / "flows.csv")
    manta = pd.read_csv(ROOT / "sandbox/data/reef_manta.csv")
    cots = pd.read_csv(ROOT / "sandbox/data/reef_cots.csv")
    previous = previous[previous.treatment == "idw_alpha015_old_coral"]
    trajectories = pd.concat([current, previous], ignore_index=True)
    flows = pd.concat([flows, previous_flows[
        previous_flows.treatment == "idw_alpha015_old_coral"]], ignore_index=True)
    expected = len(REEFS) * 34
    for name in TREATMENTS:
        if len(trajectories[trajectories.treatment == name]) != expected:
            raise ValueError(f"Incomplete treatment: {name}")
    old_control = pd.read_csv(PREVIOUS / "trajectories.csv")
    old_control = old_control[old_control.treatment == "inherited_1991"]
    new_control = current[current.treatment == "inherited_1991"]
    control_delta = np.max(np.abs(old_control.coral_cover.to_numpy() -
                                  new_control.coral_cover.to_numpy()))
    if control_delta > 1e-12:
        raise ValueError("Inherited control does not replay prior 1991 coral")
    no_cots = current[current.treatment == "coral_only_no_cots"]
    if no_cots.adults_ha.abs().max() != 0.0:
        raise ValueError("No-COTS counterfactual has nonzero adults")

    rows = []
    fig, axes = plt.subplots(4, 2, figsize=(13, 13), sharex=True)
    for row, (reef, observed) in enumerate(REEFS.items()):
        panels = axes[row]
        subset = trajectories[trajectories.reef_name == reef]
        for name, (label, color, style) in TREATMENTS.items():
            series = subset[subset.treatment == name].sort_values("year")
            panels[1].plot(series.year, series.coral_cover, color=color,
                           linestyle=style, linewidth=1.4, label=label if row == 0 else None)
            if name != "coral_only_no_cots":
                panels[0].plot(series.year.iloc[1:], series.simulated_cpue.iloc[1:],
                               color=color, linestyle=style, linewidth=1.4,
                               label=label if row == 0 else None)
        observed_cots = cots[(cots.reef_name == observed) & cots.year.between(1992, 2024)]
        panels[0].scatter(observed_cots.year, observed_cots.cotsptow,
                          s=18, color="#111827", alpha=0.7,
                          label="LTMP COTS/tow" if row == 0 else None)
        coral_obs = manta[(manta.reef_name == observed) &
                          (manta.data_type == "manta") &
                          (manta.domain_category == "reef") &
                          (manta.variable == "HC") &
                          (manta.purpose == "MANTA") &
                          (manta.project_code == "LTMP") &
                          (manta.depth == 9) &
                          manta.report_year.between(1992, 2024)]
        panels[1].scatter(coral_obs.report_year, coral_obs["median"],
                          s=18, color="#111827", alpha=0.6,
                          label="LTMP 9 m manta coral" if row == 0 else None)
        panels[0].axhline(0.22, color="#92400e", linestyle=":", linewidth=0.8)
        panels[0].set_ylabel(f"{reef}\nCOTS/tow")
        panels[1].set_ylabel("Hard-coral fraction")
        for panel in panels:
            panel.set_xlim(1991, 2024)
            panel.set_ylim(bottom=0)
            panel.grid(alpha=0.15)

        coral_only = subset[subset.treatment == "coral_only_no_cots"].set_index("year")
        for name in TREATMENTS:
            if name == "coral_only_no_cots":
                continue
            series = subset[subset.treatment == name].set_index("year")
            gap = series.loc[2004:2012]
            coral_gap = coral_only.loc[2004:2012]
            reef_flows = flows[(flows.treatment == name) & (flows.reef_name == reef) &
                               flows.year.between(2004, 2012)]
            rows.append({
                "treatment": name,
                "reef_name": reef,
                "mean_adults_ha_2004_2012": gap.adults_ha.mean(),
                "mean_model_coral_2004_2012": gap.coral_cover.mean(),
                "mean_no_cots_coral_2004_2012": coral_gap.coral_cover.mean(),
                "coral_deficit_vs_no_cots_2004_2012":
                    (coral_gap.coral_cover - gap.coral_cover).mean(),
                "model_coral_2012": series.loc[2012, "coral_cover"],
                "no_cots_coral_2012": coral_only.loc[2012, "coral_cover"],
                "mean_external_immigration_2004_2012":
                    reef_flows.external_immigration.mean(),
                "mean_internal_immigration_2004_2012":
                    reef_flows.internal_immigration.mean(),
            })
    axes[0, 0].set_title("COTS: fixed observation scale; 1991 log row excluded")
    axes[0, 1].set_title("Coral: whole-reef model vs 9 m LTMP manta")
    axes[-1, 0].set_xlabel("Year")
    axes[-1, 1].set_xlabel("Year")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    handles2, labels2 = axes[0, 1].get_legend_handles_labels()
    legend = dict(zip(labels + labels2, handles + handles2))
    fig.legend(legend.values(), legend.keys(), loc="lower center",
               ncol=3, fontsize=8, bbox_to_anchor=(0.5, -0.01))
    fig.tight_layout(rect=(0, 0.045, 1, 1))
    fig.savefig(plot, dpi=170, bbox_inches="tight")
    plt.close(fig)
    pd.DataFrame(rows).to_csv(output, index=False)
    with metadata.open("w") as stream:
        json.dump({
            "interpretation": "No-COTS control tests aggregate COTS effect, not residual-only grazing; 6:3:1 is a model-derived ratio, not observed 1991 cohorts.",
            "gap_years": [2004, 2012],
            "inherited_coral_replay_max_abs_difference": control_delta,
            "no_cots_max_abs_adults_ha": float(no_cots.adults_ha.abs().max()),
            "inputs": {name: sha256(path) for name, path in {
                "current_trajectories": run / "trajectories.csv",
                "previous_trajectories": PREVIOUS / "trajectories.csv",
                "current_flows": run / "flows.csv",
                "previous_flows": PREVIOUS / "flows.csv",
                "current_scores": run / "scores.csv",
                "manta": ROOT / "sandbox/data/reef_manta.csv",
                "cots": ROOT / "sandbox/data/reef_cots.csv",
                "script": Path(__file__),
            }.items()},
            "outputs": {"cohort_history_audit.csv": sha256(output),
                        "cots_coral_cohort_history.png": sha256(plot)},
        }, stream, indent=2)
        stream.write("\n")
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    main(parser.parse_args().run)
