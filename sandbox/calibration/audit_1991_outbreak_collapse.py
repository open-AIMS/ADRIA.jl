"""Audit the preregistered post-outbreak collapse x fast-coral growth screen on the
v0.3 leading candidate (Ext x5 + P(Q=30) + gradual Allee): ReefMod-style collapse
{none, D2, D2J, D4} x coral calibration group {g1, grme, grmefast}
(crash_mechanism_protocol.md, 2026-10-12)."""

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
from audit_1991_initiation import exact_replay  # noqa: E402
from audit_1991_switching_repair import REEFS, TWO_WAVE, sha256  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
PRIOR = ROOT / "sandbox/calibration/runs/20261012T_11_predation_capacity_1991_v1"
REFERENCE, PRIOR_REFERENCE = "g1_none", "ext5_PA_q30"
CORALS = ("g1", "grme", "grmefast")
COLLAPSES = ("none", "D2", "D2J", "D4")
ARMS = [f"{coral}_{collapse}" for coral in CORALS for collapse in COLLAPSES]
COLORS = {"none": "#475569", "D2": "#dc2626", "D2J": "#7c3aed", "D4": "#0284c7"}
STYLES = {"g1": "-", "grme": "--", "grmefast": ":"}
GAP = (2004, 2010)
ONSET = (2008, 2013)
SECOND_WINDOW = (2010, 2017)
N_YEARS = 34


def main(run: Path) -> None:
    paths = {name: run / name for name in (
        "outbreak_collapse_audit.csv", "outbreak_collapse_audit.json",
        "cots_coral_outbreak_collapse.png", "fast_coral_outbreak_collapse.png")}
    if any(path.exists() for path in paths.values()):
        raise FileExistsError("Refusing to overwrite outbreak-collapse audit")
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
        replay[name] = exact_replay(tables[name][tables[name].treatment == REFERENCE],
                                    prior[prior.treatment == PRIOR_REFERENCE], columns)
    prior_flows = pd.read_csv(PRIOR / "flows.csv")
    reference = prior_flows[prior_flows.treatment == PRIOR_REFERENCE].sort_values(["reef_name", "year"])
    external = {}
    for arm in ARMS:
        current = flows[flows.treatment == arm].sort_values(["reef_name", "year"])
        external[arm] = float(np.max(np.abs(current.external_immigration.to_numpy() -
                                            reference.external_immigration.to_numpy())))
        if external[arm] > 1e-9:
            raise ValueError(f"External immigration differs from Ext x5 in {arm}")
    numeric = mechanisms[mechanisms.year >= 1992][[
        "recruits_ha", "subadults_ha", "adults_ha", "fast_coral", "total_coral",
        "food_survival", "body_condition", "collapse_removed_ha"]].to_numpy()
    if not np.isfinite(numeric).all() or (numeric < -1e-12).any():
        raise ValueError("Nonfinite or negative state/diagnostic")
    removed = mechanisms.groupby("treatment").collapse_removed_ha.sum()
    if (removed[[a for a in ARMS if a.endswith("_none")]] != 0.0).any():
        raise ValueError("Collapse removals in an arm with the switch off")
    if (removed[[a for a in ARMS if "_D2" in a]] <= 0.0).any():
        raise ValueError("D2 collapse never triggered (wiring failure)")

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
                          (mechanisms.reef_name == score.reef_name)].sort_values("year")
        obs = observations[observations.reef_name == REEFS[score.reef_name]]
        reef_flows = flows[(flows.treatment == score.treatment) &
                           (flows.reef_name == score.reef_name)]
        years = [int(v) for v in str(score.sim_peak_years).split(";") if v and v != "nan"]
        heights = [float(v) for v in str(score.sim_peak_heights_cpue).split(";") if v and v != "nan"]
        ratio = float("nan")
        if len(years) >= 2:
            smooth = series.simulated_cpue.rolling(3, center=True, min_periods=1).mean()
            ratio = float(smooth[series.year.between(years[0], years[1])].min() / min(heights[:2]))
        first_wave = series[series.year.between(1995, 2006)].reset_index(drop=True)
        cpue = first_wave.simulated_cpue.to_numpy()
        steps = cpue[1:] / np.maximum(cpue[:-1], 1e-12)
        drop_index = int(np.argmin(steps))
        collapse_years = mech[mech.collapse_removed_ha > 0.0].year.tolist()
        gap = mech[mech.year.between(*GAP)]
        onset = mech[mech.year.between(*ONSET)]
        gap_flows = reef_flows[reef_flows.year.between(*GAP)]
        rows.append({
            "treatment": score.treatment, "coral": score.treatment.split("_")[0],
            "collapse": score.treatment.split("_")[1], "reef_name": score.reef_name,
            "first_peak_year": years[0] if years else float("nan"),
            "second_peak_year": years[1] if len(years) > 1 else float("nan"),
            "first_peak_cpue": heights[0] if heights else float("nan"),
            "second_peak_cpue": heights[1] if len(heights) > 1 else float("nan"),
            "max_cpue_2010_2017": float(series[series.year.between(*SECOND_WINDOW)].simulated_cpue.max()),
            "observed_max_cpue_2010_2017": float(obs[obs.year.between(*SECOND_WINDOW)].cotsptow.max()),
            "trough_to_smaller_peak": ratio,
            "matched_peaks": int(score.matched_peaks),
            "smoothed_loss": float(score.loss),
            "collapse_years": ";".join(str(y) for y in collapse_years),
            "collapse_removed_ha": float(mech.collapse_removed_ha.sum()),
            "steepest_one_year_drop_ratio_1995_2006": float(steps[drop_index]),
            "steepest_drop_year": int(first_wave.year[drop_index + 1]),
            "min_smoothed_adults_1999_2012": float(mech[mech.year.between(1999, 2012)].adults_ha.rolling(3, center=True, min_periods=1).mean().min()),
            "mean_gap_adults_ha": gap.adults_ha.mean(),
            "mean_gap_settlement_ha": float(gap_flows.settled_recruits.mean()),
            "peak_fast_coral_2008_2013": float(onset.fast_coral.max()),
            "max_condition_2008_2013": float(onset.body_condition.max()),
            "model_coral_site_fraction_2012": float(mech[mech.year == 2012].total_coral_site_fraction.iloc[0]),
            "model_coral_site_fraction_2016": float(mech[mech.year == 2016].total_coral_site_fraction.iloc[0]),
            "observed_manta_coral_2010_2014": float(manta[(manta.reef_name == REEFS[score.reef_name]) &
                manta.report_year.between(2010, 2014)]["median"].mean()),
            "coral_site_rmse_vs_ltmp": float(np.sqrt(np.nanmean([
                (float(mech[mech.year == y].total_coral_site_fraction.iloc[0]) - v) ** 2
                for y, v in manta[(manta.reef_name == REEFS[score.reef_name]) &
                                  manta.report_year.between(1992, 2024)]
                .groupby("report_year")["median"].mean().items()]))),
        })
    result = pd.DataFrame(rows)
    result["order"] = result.treatment.map({a: i for i, a in enumerate(ARMS)})
    result = result.sort_values(["order", "reef_name"]).drop(columns="order")
    result.to_csv(paths["outbreak_collapse_audit.csv"], index=False)

    two_wave = result[result.reef_name.isin(TWO_WAVE)]
    summary = result.groupby("treatment").agg(
        summed_loss=("smoothed_loss", "sum"),
        reefs_with_two_peaks=("second_peak_year", "count"),
        matched_observed_peaks_of_seven=("matched_peaks", "sum"),
        mean_first_peak_year=("first_peak_year", "mean"),
        mean_first_peak_cpue=("first_peak_cpue", "mean"),
        mean_steepest_drop_ratio=("steepest_one_year_drop_ratio_1995_2006", "mean"),
        reefs_collapsing=("collapse_removed_ha", lambda v: int((v > 0).sum())),
        mean_min_trough_adults_ha=("min_smoothed_adults_1999_2012", "mean"),
        mean_gap_adults_ha=("mean_gap_adults_ha", "mean"),
        mean_max_cpue_2010_2017=("max_cpue_2010_2017", "mean"),
        mean_peak_fast_coral_2008_2013=("peak_fast_coral_2008_2013", "mean"),
        mean_coral_site_rmse=("coral_site_rmse_vs_ltmp", "mean"),
    ).reindex(ARMS)
    gate = two_wave.groupby("treatment").agg(
        two_wave_reefs_with_two_peaks=("second_peak_year", "count"),
        worst_available_trough_ratio=("trough_to_smaller_peak", "max"))
    summary = summary.join(gate)
    summary["passes_qualitative_gate"] = (
        (summary.two_wave_reefs_with_two_peaks == len(TWO_WAVE)) &
        (summary.matched_observed_peaks_of_seven == 7) &
        (summary.worst_available_trough_ratio <= .10))

    fig, axes = plt.subplots(4, 2, figsize=(14, 13), sharex=True)
    fast_fig, fast_axes = plt.subplots(4, 1, figsize=(8, 12), sharex=True)
    for index, (reef, observed) in enumerate(REEFS.items()):
        cots_ax, coral_ax = axes[index]
        for arm in ARMS:
            coral, collapse = arm.split("_")
            values = trajectories[(trajectories.treatment == arm) &
                                  (trajectories.reef_name == reef)].sort_values("year")
            mech = mechanisms[(mechanisms.treatment == arm) &
                              (mechanisms.reef_name == reef)].sort_values("year")
            label = f"{coral} / {collapse}" if index == 0 else None
            width = 1.8 if arm == REFERENCE else 1.2
            cots_ax.plot(values.year.iloc[1:], values.simulated_cpue.iloc[1:], color=COLORS[collapse],
                         linestyle=STYLES[coral], linewidth=width, label=label)
            coral_ax.plot(mech.year, mech.total_coral_site_fraction, color=COLORS[collapse],
                          linestyle=STYLES[coral], linewidth=width)
            fast_axes[index].plot(mech.year, mech.fast_coral, color=COLORS[collapse],
                                  linestyle=STYLES[coral], linewidth=width, label=label)
        obs_cots = observations[(observations.reef_name == observed) &
                                observations.year.between(1992, 2024)]
        cots_ax.scatter(obs_cots.year, obs_cots.cotsptow, s=16, color="#111827", alpha=.7, zorder=5,
                        label="LTMP COTS/tow" if index == 0 else None)
        obs_coral = manta[(manta.reef_name == observed) & manta.report_year.between(1992, 2024)]
        coral_ax.scatter(obs_coral.report_year, obs_coral["median"], s=16, color="#111827", alpha=.6, zorder=5)
        cots_ax.set_ylabel(f"{reef}\nCOTS/tow")
        coral_ax.set_ylabel("Coral fraction of site area")
        cots_ax.axhline(.22, color="#92400e", linewidth=.8, linestyle=":")
        cots_ax.axhline(67.5 * 0.015, color="#dc2626", linewidth=.6, linestyle="-.", alpha=.6)
        fast_axes[index].axhline(0.05, color="#92400e", linewidth=.8, linestyle=":")
        fast_axes[index].axhline(0.10, color="#16a34a", linewidth=.8, linestyle=":")
        fast_axes[index].set_ylabel(f"{reef}\nfast coral (habitable)")
        for axis in (cots_ax, coral_ax, fast_axes[index]):
            axis.set_xlim(1991, 2024)
            axis.set_ylim(bottom=0)
            axis.grid(alpha=.15)
    axes[0, 0].set_title("COTS at q = 0.015 (red dash-dot: collapse threshold 67.5/ha)")
    axes[0, 1].set_title("Model coral (site-area fraction) / LTMP 9 m manta")
    fast_axes[0].set_title("Fast coral (habitable fraction); dotted: starvation 0.05, condition saturation 0.10")
    axes[-1, 0].set_xlabel("Year")
    axes[-1, 1].set_xlabel("Year")
    fast_axes[-1].set_xlabel("Year")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=7, fontsize=8, bbox_to_anchor=(.5, -.01))
    fig.tight_layout(rect=(0, .05, 1, 1))
    fig.savefig(paths["cots_coral_outbreak_collapse.png"], dpi=170, bbox_inches="tight")
    handles, labels = fast_axes[0].get_legend_handles_labels()
    fast_fig.legend(handles, labels, loc="lower center", ncol=4, fontsize=7, bbox_to_anchor=(.5, -.01))
    fast_fig.tight_layout(rect=(0, .06, 1, 1))
    fast_fig.savefig(paths["fast_coral_outbreak_collapse.png"], dpi=170, bbox_inches="tight")
    plt.close("all")

    report = {
        "status": "bounded_diagnostic_not_promoted",
        "reference_replay_max_rel_delta": replay,
        "external_immigration_max_abs_delta_vs_reference": external,
        "collapse_removed_ha_total_by_arm": json.loads(removed.to_json()),
        "treatments": json.loads(summary.to_json(orient="index")),
        "criteria": "2 peaks on each of Lizard/MacGillivray/North Direction; 7/7 matched; "
                    "worst trough ratio <=0.10; absolute CPUE under q=0.015",
        "inputs": {name: sha256(run / name) for name in (
            "metadata.toml", "scores.csv", "trajectories.csv", "flows.csv", "mechanisms.csv")},
        "script_sha256": sha256(Path(__file__)),
        "outputs": {name: sha256(path) for name, path in paths.items()
                    if name != "outbreak_collapse_audit.json"},
    }
    paths["outbreak_collapse_audit.json"].write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    pd.set_option("display.width", 250)
    print(summary.round(3).to_string())
    print("Reference replay max delta:", max(replay.values()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    main(parser.parse_args().run)
