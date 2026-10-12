"""Audit the preregistered benthic habitable fraction x historical cyclone screen
on the condition model (crash_mechanism_protocol.md, 2026-10-12)."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from audit_1991_switching_repair import REEFS, TWO_WAVE, exact_replay, sha256  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
PRIOR = ROOT / "sandbox/calibration/runs/20261012T_04_habitat_disturbance_1991_v2"
BASE = "k05_nocyc"
PRIOR_BASE = "k05_bleaching"
ARMS = {
    BASE: ("k = 0.5, no cyclones (base)", "#475569", "-"),
    "k05_cyc": ("k = 0.5, RRN cyclones", "#475569", "--"),
    "ktotal_nocyc": ("k benthic total", "#0d9488", "-"),
    "ktotal_cyc": ("k benthic total + cyclones", "#0d9488", "--"),
    "kmapped_nocyc": ("k benthic mapped", "#dc2626", "-"),
    "kmapped_cyc": ("k benthic mapped + cyclones", "#dc2626", "--"),
}
GAP = (2004, 2010)
ONSET = (2008, 2011)
SECOND_WINDOW = (2010, 2017)
N_YEARS = 34


def main(run: Path) -> None:
    paths = {name: run / name for name in (
        "benthic_cyclone_audit.csv", "benthic_cyclone_audit.json",
        "cots_coral_benthic_cyclone.png")}
    if any(path.exists() for path in paths.values()):
        raise FileExistsError("Refusing to overwrite benthic/cyclone audit")
    tables = {name: pd.read_csv(run / f"{name}.csv")
              for name in ("trajectories", "flows", "mechanisms", "scores", "failures")}
    if len(tables["failures"]):
        raise ValueError("Failed evaluations in run")
    trajectories, flows, mechanisms = (tables[k] for k in ("trajectories", "flows", "mechanisms"))
    if set(trajectories.treatment) != set(ARMS):
        raise ValueError("Missing or extra treatment")
    if any(len(t) != len(ARMS) * len(REEFS) * N_YEARS for t in (trajectories, flows, mechanisms)):
        raise ValueError("Incomplete reef-year outputs")
    replay = {}
    for name in ("trajectories", "flows", "mechanisms"):
        prior = pd.read_csv(PRIOR / f"{name}.csv")
        columns = [c for c in prior if c not in ("treatment", "reef_name", "year")]
        replay[name] = exact_replay(tables[name][tables[name].treatment == BASE],
                                    prior[prior.treatment == PRIOR_BASE], columns)
    reference = flows[flows.treatment == BASE].sort_values(["reef_name", "year"])
    external = {}
    for arm in ARMS:
        current = flows[flows.treatment == arm].sort_values(["reef_name", "year"])
        external[arm] = float(np.max(np.abs(current.external_immigration.to_numpy() -
                                            reference.external_immigration.to_numpy())))
        if external[arm] > 1e-9:
            raise ValueError(f"External immigration changed in {arm}")
    numeric = mechanisms[mechanisms.year >= 1992][[
        "recruits_ha", "subadults_ha", "adults_ha", "fast_coral", "total_coral",
        "food_survival", "body_condition"]].to_numpy()
    if not np.isfinite(numeric).all() or (numeric < -1e-12).any():
        raise ValueError("Nonfinite or negative state/diagnostic")

    observations = pd.read_csv(ROOT / "sandbox/data/reef_cots.csv")
    manta = pd.read_csv(ROOT / "sandbox/data/reef_manta.csv")
    manta = manta[(manta.data_type == "manta") & (manta.domain_category == "reef") &
                  (manta.variable == "HC") & (manta.purpose == "MANTA") &
                  (manta.project_code == "LTMP") & (manta.depth == 9)]
    rows = []
    for score in tables["scores"][tables["scores"].observation_treatment == "smoothed_3y"].itertuples(index=False):
        series = trajectories[(trajectories.treatment == score.treatment) &
                              (trajectories.reef_name == score.reef_name) &
                              (trajectories.year >= 1992)].sort_values("year")
        mech = mechanisms[(mechanisms.treatment == score.treatment) &
                          (mechanisms.reef_name == score.reef_name)]
        obs = observations[observations.reef_name == REEFS[score.reef_name]]
        obs_coral = manta[(manta.reef_name == REEFS[score.reef_name]) &
                          manta.report_year.between(2010, 2014)]
        years = [int(v) for v in str(score.sim_peak_years).split(";") if v and v != "nan"]
        heights = [float(v) for v in str(score.sim_peak_heights_cpue).split(";") if v and v != "nan"]
        ratio = float("nan")
        if len(years) >= 2:
            smooth = series.simulated_cpue.rolling(3, center=True, min_periods=1).mean()
            ratio = float(smooth[series.year.between(years[0], years[1])].min() / min(heights[:2]))
        gap = mech[mech.year.between(*GAP)]
        rows.append({
            "treatment": score.treatment, "reef_name": score.reef_name,
            "first_peak_year": years[0] if years else float("nan"),
            "second_peak_year": years[1] if len(years) > 1 else float("nan"),
            "first_peak_cpue": heights[0] if heights else float("nan"),
            "second_peak_cpue": heights[1] if len(heights) > 1 else float("nan"),
            "max_cpue_2010_2017": float(series[series.year.between(*SECOND_WINDOW)].simulated_cpue.max()),
            "observed_max_cpue_2010_2017": float(obs[obs.year.between(*SECOND_WINDOW)].cotsptow.max()),
            "trough_to_smaller_peak": ratio,
            "matched_peaks": int(score.matched_peaks),
            "smoothed_loss": float(score.loss),
            "mean_gap_adults_ha": gap.adults_ha.mean(),
            "mean_fast_coral_onset_2008_2011": mech[mech.year.between(*ONSET)].fast_coral.mean(),
            "model_coral_2012": float(series[series.year == 2012].coral_cover.iloc[0]),
            "model_coral_2016": float(series[series.year == 2016].coral_cover.iloc[0]),
            "model_coral_site_fraction_2012": float(mech[mech.year == 2012].total_coral_site_fraction.iloc[0]),
            "model_coral_site_fraction_2016": float(mech[mech.year == 2016].total_coral_site_fraction.iloc[0]),
            "observed_manta_coral_2016_2018": float(manta[(manta.reef_name == REEFS[score.reef_name]) &
                manta.report_year.between(2016, 2018)]["median"].mean()),
            "observed_manta_coral_2010_2014": float(obs_coral["median"].mean()) if len(obs_coral) else float("nan"),
        })
    result = pd.DataFrame(rows)
    result.to_csv(paths["benthic_cyclone_audit.csv"], index=False)
    two_wave = result[result.reef_name.isin(TWO_WAVE)]
    summary = two_wave.groupby("treatment").agg(
        two_wave_reefs_with_two_peaks=("second_peak_year", "count"),
        worst_available_trough_ratio=("trough_to_smaller_peak", "max"),
        mean_max_cpue_2010_2017=("max_cpue_2010_2017", "mean"),
        mean_fast_coral_onset=("mean_fast_coral_onset_2008_2011", "mean"),
        mean_gap_adults_ha=("mean_gap_adults_ha", "mean"),
        mean_coral_2012=("model_coral_2012", "mean"),
        mean_coral_2016=("model_coral_2016", "mean"),
        mean_coral_site_2012=("model_coral_site_fraction_2012", "mean"),
        mean_coral_site_2016=("model_coral_site_fraction_2016", "mean"),
    )
    summary["matched_observed_peaks_of_seven"] = result.groupby("treatment").matched_peaks.sum()
    summary["passes_qualitative_gate"] = (
        (summary.two_wave_reefs_with_two_peaks == len(TWO_WAVE)) &
        (summary.matched_observed_peaks_of_seven == 7) &
        (summary.worst_available_trough_ratio <= .10))

    fig, axes = plt.subplots(4, 2, figsize=(14, 13), sharex=True)
    for index, (reef, observed) in enumerate(REEFS.items()):
        cots_ax, coral_ax = axes[index]
        for arm, (label, color, style) in ARMS.items():
            values = trajectories[(trajectories.treatment == arm) &
                                  (trajectories.reef_name == reef)].sort_values("year")
            fast = mechanisms[(mechanisms.treatment == arm) &
                              (mechanisms.reef_name == reef)].sort_values("year")
            cots_ax.plot(values.year.iloc[1:], values.simulated_cpue.iloc[1:], color=color,
                         linestyle=style, linewidth=1.4, label=label if index == 0 else None)
            coral_ax.plot(values.year, values.coral_cover, color=color, linestyle=style, linewidth=1.4)
            coral_ax.plot(fast.year, fast.fast_coral, color=color, linestyle=":", linewidth=.9)
        obs_cots = observations[(observations.reef_name == observed) &
                                observations.year.between(1992, 2024)]
        cots_ax.scatter(obs_cots.year, obs_cots.cotsptow, s=16, color="#111827", alpha=.7,
                        label="LTMP COTS/tow" if index == 0 else None)
        obs_coral = manta[(manta.reef_name == observed) & manta.report_year.between(1992, 2024)]
        coral_ax.scatter(obs_coral.report_year, obs_coral["median"], s=16, color="#111827", alpha=.6)
        cots_ax.set_ylabel(f"{reef}\nCOTS/tow")
        coral_ax.set_ylabel("Coral fraction (solid total, dotted fast)")
        cots_ax.axhline(.22, color="#92400e", linewidth=.8, linestyle=":")
        for axis in (cots_ax, coral_ax):
            axis.set_xlim(1991, 2024)
            axis.set_ylim(bottom=0)
            axis.grid(alpha=.15)
    axes[0, 0].set_title("COTS at fixed q = 0.015 COTS/tow per COTS/ha")
    axes[0, 1].set_title("Model coral (habitable fraction) / LTMP 9 m manta")
    axes[-1, 0].set_xlabel("Year")
    axes[-1, 1].set_xlabel("Year")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, fontsize=8, bbox_to_anchor=(.5, -.01))
    fig.tight_layout(rect=(0, .05, 1, 1))
    fig.savefig(paths["cots_coral_benthic_cyclone.png"], dpi=170, bbox_inches="tight")
    plt.close(fig)

    report = {
        "status": "bounded_diagnostic_not_promoted",
        "base_replay_of_habitat_v2_base_max_abs_delta": replay,
        "external_immigration_max_abs_delta_vs_base": external,
        "treatments": json.loads(summary.to_json(orient="index")),
        "criteria": "2 peaks on each of Lizard/MacGillivray/North Direction; 7/7 matched; "
                    "worst trough ratio <=0.10; absolute CPUE under q=0.015",
        "inputs": {name: sha256(run / name) for name in (
            "metadata.toml", "scores.csv", "trajectories.csv", "flows.csv", "mechanisms.csv")},
        "script_sha256": sha256(Path(__file__)),
        "outputs": {name: sha256(path) for name, path in paths.items()
                    if name != "benthic_cyclone_audit.json"},
    }
    paths["benthic_cyclone_audit.json"].write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    pd.set_option("display.width", 250)
    print(summary.round(3).to_string())
    print("Base replay max delta:", max(replay.values()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    main(parser.parse_args().run)
