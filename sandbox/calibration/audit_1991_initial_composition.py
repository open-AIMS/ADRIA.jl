"""Audit the preregistered data-derived 1991 coral composition screen
(crash_mechanism_protocol.md, 2026-10-12): inherited 2023 composition versus the
observed 1993-95 LTMP photo-transect fast share, with Q re-selected under the observed
composition and a fast-growth multiplier sensitivity."""

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
from audit_1991_growth_starvation import PHOTO, PHOTO_REEFS, observed_fast, slope  # noqa: E402
from audit_1991_initiation import exact_replay  # noqa: E402
from audit_1991_switching_repair import REEFS, TWO_WAVE, sha256  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
PRIOR = ROOT / "sandbox/calibration/runs/20261012T_13_growth_starvation_1991_v1"
REFERENCE, PRIOR_REFERENCE = "inh_m110_q30", "m110_ref"
Q_GRID = (20, 30, 45, 70, 100)
ARMS = [REFERENCE] + [f"obs_m110_q{q}" for q in Q_GRID] + ["obs_m100_q30", "obs_m120_q30"]
SECOND_WINDOW = (2010, 2017)
N_YEARS = 34


def main(run: Path) -> None:
    paths = {name: run / name for name in (
        "initial_composition_audit.csv", "initial_composition_audit.json",
        "cots_fast_coral_initial_composition.png")}
    if any(path.exists() for path in paths.values()):
        raise FileExistsError("Refusing to overwrite initial-composition audit")
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
    numeric = mechanisms[mechanisms.year >= 1991][[
        "recruits_ha", "subadults_ha", "adults_ha", "fast_coral", "total_coral",
        "food_survival", "body_condition"]].to_numpy()
    if not np.isfinite(numeric).all() or (numeric < -1e-12).any():
        raise ValueError("Nonfinite or negative state/diagnostic")
    # Totals preserved: the 1991 total cover must match between compositions.
    start = mechanisms[mechanisms.year == 1991].pivot(index="reef_name", columns="treatment",
                                                       values="total_coral")
    total_delta = float((start[ARMS[1:]].sub(start[REFERENCE], axis=0)).abs().to_numpy().max())
    if total_delta > 1e-9:
        raise ValueError(f"Initial total cover differs between compositions: {total_delta}")

    observed = observed_fast()
    cots_obs = pd.read_csv(ROOT / "sandbox/data/reef_cots.csv")
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
                          (mechanisms.reef_name == score.reef_name)].sort_values("year").set_index("year")
        obs = cots_obs[cots_obs.reef_name == REEFS[score.reef_name]]
        years = [int(v) for v in str(score.sim_peak_years).split(";") if v and v != "nan"]
        heights = [float(v) for v in str(score.sim_peak_heights_cpue).split(";") if v and v != "nan"]
        ratio = float("nan")
        if len(years) >= 2:
            smooth = series.simulated_cpue.rolling(3, center=True, min_periods=1).mean()
            ratio = float(smooth[series.year.between(years[0], years[1])].min() / min(heights[:2]))
        first_wave = series[series.year.between(1995, 2006)].reset_index(drop=True)
        cpue = first_wave.simulated_cpue.to_numpy()
        steps = cpue[1:] / np.maximum(cpue[:-1], 1e-12)
        manta_reef = manta[(manta.reef_name == REEFS[score.reef_name]) &
                           manta.report_year.between(1992, 2024)].groupby("report_year")["median"].mean()
        row = {
            "treatment": score.treatment, "reef_name": score.reef_name,
            "composition": score.treatment.split("_")[0],
            "multiplier": int(score.treatment.split("_")[1][1:]) / 100,
            "Q": float(score.treatment.split("_q")[1]),
            "first_peak_year": years[0] if years else float("nan"),
            "second_peak_year": years[1] if len(years) > 1 else float("nan"),
            "first_peak_cpue": heights[0] if heights else float("nan"),
            "second_peak_cpue": heights[1] if len(heights) > 1 else float("nan"),
            "cpue_1994": float(series[series.year == 1994].simulated_cpue.iloc[0]),
            "cpue_1998": float(series[series.year == 1998].simulated_cpue.iloc[0]),
            "max_cpue_2010_2017": float(series[series.year.between(*SECOND_WINDOW)].simulated_cpue.max()),
            "observed_max_cpue_2010_2017": float(obs[obs.year.between(*SECOND_WINDOW)].cotsptow.max()),
            "trough_to_smaller_peak": ratio,
            "matched_peaks": int(score.matched_peaks),
            "smoothed_loss": float(score.loss),
            "steepest_one_year_drop_ratio_1995_2006": float(steps.min()),
            "min_smoothed_adults_1999_2012": float(mech.loc[1999:2012].adults_ha.rolling(3, center=True, min_periods=1).mean().min()),
            "fast_coral_1991": float(mech.loc[1991, "fast_coral"]),
            "fast_share_1991": float(mech.loc[1991, "fast_coral"] / mech.loc[1991, "total_coral"]),
            "min_fast_coral_1998_2003": float(mech.loc[1998:2003].fast_coral.min()),
            "peak_fast_coral_2008_2013": float(mech.loc[2008:2013].fast_coral.max()),
            "coral_site_rmse_vs_ltmp": float(np.sqrt(np.mean(
                (mech.total_coral_site_fraction.loc[manta_reef.index] - manta_reef) ** 2))),
        }
        if score.reef_name in PHOTO_REEFS:
            fast = observed[(observed.model_reef == score.reef_name) & observed.year.between(1991, 2024)]
            late = fast[fast.year.between(2017, 2024)]
            row.update({
                "fast_rmse_vs_photo": float(np.sqrt(np.mean(
                    (mech.fast_coral.loc[fast.year].to_numpy() - fast.observed_fast.to_numpy()) ** 2))),
                "observed_slope_2017_2024": slope(late.year, late.observed_fast),
                "model_slope_2017_2024": slope(range(2017, 2025), mech.fast_coral.loc[2017:2024]),
            })
        rows.append(row)
    result = pd.DataFrame(rows)
    result["order"] = result.treatment.map({a: i for i, a in enumerate(ARMS)})
    result = result.sort_values(["order", "reef_name"]).drop(columns="order")
    result.to_csv(paths["initial_composition_audit.csv"], index=False)

    two_wave = result[result.reef_name.isin(TWO_WAVE)]
    summary = result.groupby("treatment").agg(
        summed_loss=("smoothed_loss", "sum"),
        reefs_with_first_peak=("first_peak_year", "count"),
        reefs_with_two_peaks=("second_peak_year", "count"),
        matched_observed_peaks_of_seven=("matched_peaks", "sum"),
        mean_first_peak_year=("first_peak_year", "mean"),
        mean_first_peak_cpue=("first_peak_cpue", "mean"),
        mean_steepest_drop_ratio=("steepest_one_year_drop_ratio_1995_2006", "mean"),
        mean_min_trough_adults_ha=("min_smoothed_adults_1999_2012", "mean"),
        mean_max_cpue_2010_2017=("max_cpue_2010_2017", "mean"),
        mean_fast_share_1991=("fast_share_1991", "mean"),
        mean_min_fast_1998_2003=("min_fast_coral_1998_2003", "mean"),
        mean_fast_rmse_vs_photo=("fast_rmse_vs_photo", "mean"),
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
    grid = summary.loc[[f"obs_m110_q{q}" for q in Q_GRID]]
    eligible = grid[grid.reefs_with_first_peak == len(REEFS)]
    q_star = None
    if len(eligible):
        best = eligible.summed_loss.min()
        q_star = min(float(arm.split("_q")[1]) for arm in eligible.index
                     if eligible.loc[arm, "summed_loss"] == best)
    selected = f"obs_m110_q{int(q_star)}" if q_star is not None else None

    shown = [REFERENCE, "obs_m110_q30"] + ([selected] if selected and selected != "obs_m110_q30" else []) + \
        ["obs_m100_q30", "obs_m120_q30"]
    colors = {REFERENCE: "#475569", "obs_m110_q30": "#dc2626", "obs_m100_q30": "#f59e0b",
              "obs_m120_q30": "#7c3aed"}
    if selected:
        colors.setdefault(selected, "#0d9488")
    fig, axes = plt.subplots(4, 2, figsize=(14, 13), sharex=True)
    for index, (reef, observed_name) in enumerate(REEFS.items()):
        cots_ax, fast_ax = axes[index]
        for arm in shown:
            values = trajectories[(trajectories.treatment == arm) &
                                  (trajectories.reef_name == reef)].sort_values("year")
            mech = mechanisms[(mechanisms.treatment == arm) & (mechanisms.reef_name == reef)].sort_values("year")
            label = arm + (" (selected Q)" if arm == selected else "") if index == 0 else None
            width = 2.2 if arm in (REFERENCE, selected) else 1.2
            cots_ax.plot(values.year.iloc[1:], values.simulated_cpue.iloc[1:], color=colors[arm],
                         linewidth=width, label=label)
            fast_ax.plot(mech.year, mech.fast_coral, color=colors[arm], linewidth=width)
        obs_cots = cots_obs[(cots_obs.reef_name == observed_name) & cots_obs.year.between(1992, 2024)]
        cots_ax.scatter(obs_cots.year, obs_cots.cotsptow, s=16, color="#111827", alpha=.7, zorder=5,
                        label="LTMP COTS/tow" if index == 0 else None)
        if reef in PHOTO_REEFS:
            fast = observed[observed.model_reef == reef]
            fast_ax.scatter(fast.year, fast.observed_fast, s=16, color="#111827", alpha=.7, zorder=5)
        cots_ax.set_ylabel(f"{reef}\nCOTS/tow")
        fast_ax.set_ylabel("Fast coral (habitable)")
        cots_ax.axhline(.22, color="#92400e", linewidth=.8, linestyle=":")
        fast_ax.axhline(.05, color="#92400e", linewidth=.8, linestyle=":")
        for axis in (cots_ax, fast_ax):
            axis.set_xlim(1991, 2024)
            axis.set_ylim(bottom=0)
            axis.grid(alpha=.15)
    axes[0, 0].set_title("COTS at q = 0.015: inherited 2023 vs observed 1993-95 composition")
    axes[0, 1].set_title("Fast coral / LTMP photo transects (dotted: starvation threshold 0.05)")
    axes[-1, 0].set_xlabel("Year")
    axes[-1, 1].set_xlabel("Year")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=6, fontsize=8, bbox_to_anchor=(.5, -.01))
    fig.tight_layout(rect=(0, .04, 1, 1))
    fig.savefig(paths["cots_fast_coral_initial_composition.png"], dpi=170, bbox_inches="tight")
    plt.close(fig)

    report = {
        "status": "bounded_diagnostic_not_promoted",
        "reference_replay_max_rel_delta": replay,
        "external_immigration_max_abs_delta_vs_reference": external,
        "initial_total_cover_max_abs_delta": total_delta,
        "selected_Q_under_observed_composition": q_star, "selected_arm": selected,
        "treatments": json.loads(summary.to_json(orient="index")),
        "criteria": "2 peaks on each of Lizard/MacGillivray/North Direction; 7/7 matched; "
                    "worst trough ratio <=0.10; absolute CPUE under q=0.015",
        "inputs": {name: sha256(run / name) for name in (
            "metadata.toml", "scores.csv", "trajectories.csv", "flows.csv", "mechanisms.csv")} |
                  {"reef_photo_transect.csv": sha256(PHOTO)},
        "script_sha256": sha256(Path(__file__)),
        "outputs": {name: sha256(path) for name, path in paths.items()
                    if name != "initial_composition_audit.json"},
    }
    paths["initial_composition_audit.json"].write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    pd.set_option("display.width", 250)
    print(summary.round(3).to_string())
    print("Selected Q under observed composition:", q_star)
    print("Reference replay max delta:", max(replay.values()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    main(parser.parse_args().run)
