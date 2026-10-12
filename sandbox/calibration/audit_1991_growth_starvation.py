"""Audit the preregistered fast-coral growth fit x faster starvation screen
(crash_mechanism_protocol.md, 2026-10-12): per-reef RME calibration groups (group 5
for the calibration reefs), fast-group linear-extension multiplier m x starvation
{ref, M0 (memory 0), PG (post-grazing preferred food)}.

Stage 1 selects m* from the `ref` arms by the frozen coral criterion: summed squared
error of the 2017-2024 least-squares fast-cover slope against LTMP photo transects
(Acropora + Pocilloporidae, 9 m). Stage 2 reports COTS metrics at m*."""

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
PRIOR = ROOT / "sandbox/calibration/runs/20261012T_12_outbreak_collapse_1991_v1"
PHOTO = ROOT / "sandbox/data/reef_photo_transect.csv"
REFERENCE, PRIOR_REFERENCE = "m100_ref", "grme_none"
MULTIPLIERS = (100, 105, 110, 115, 120)
STARVATION = ("ref", "M0", "PG")
ARMS = [f"m{m}_{s}" for m in MULTIPLIERS for s in STARVATION]
PHOTO_REEFS = {"Lizard Island Reef": "Lizard Isles", "MacGillivray Reef": "Macgillivray Reef",
               "North Direction Reef": "North Direction Island"}
FAST_TAXA = ("Acropora", "Pocilloporidae")
RECOVERY = (2017, 2024)
EARLY_RECOVERY = (2001, 2013)
ONSET = (2008, 2013)
SECOND_WINDOW = (2010, 2017)
N_YEARS = 34
COLORS = {"ref": "#475569", "M0": "#0284c7", "PG": "#dc2626"}


def slope(years, values) -> float:
    years, values = np.asarray(years, float), np.asarray(values, float)
    return float(np.polyfit(years, values, 1)[0]) if len(years) >= 2 else float("nan")


def observed_fast() -> pd.DataFrame:
    photo = pd.read_csv(PHOTO, low_memory=False)
    photo = photo[(photo.purpose == "COMPOSITION") & (photo.variable == "HARD CORAL") &
                  (photo.depth == 9) & photo.reef_name.isin(PHOTO_REEFS.values()) &
                  photo.reefpage_category.isin(FAST_TAXA)]
    fast = photo.groupby(["reef_name", "report_year"])["median"].sum().reset_index()
    fast["model_reef"] = fast.reef_name.map({v: k for k, v in PHOTO_REEFS.items()})
    return fast.rename(columns={"report_year": "year", "median": "observed_fast"})


def main(run: Path) -> None:
    paths = {name: run / name for name in (
        "growth_starvation_audit.csv", "growth_starvation_coral_fit.csv",
        "growth_starvation_audit.json", "fast_coral_growth_fit.png",
        "cots_coral_growth_starvation.png")}
    if any(path.exists() for path in paths.values()):
        raise FileExistsError("Refusing to overwrite growth-starvation audit")
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
        "food_survival", "body_condition"]].to_numpy()
    if not np.isfinite(numeric).all() or (numeric < -1e-12).any():
        raise ValueError("Nonfinite or negative state/diagnostic")

    # Stage 1: coral fit -----------------------------------------------------------
    observed = observed_fast()
    fit_rows = []
    for arm in ARMS:
        for reef, photo_reef in PHOTO_REEFS.items():
            model = mechanisms[(mechanisms.treatment == arm) & (mechanisms.reef_name == reef)].set_index("year")
            obs = observed[observed.model_reef == reef]
            late = obs[obs.year.between(*RECOVERY)]
            early = obs[obs.year.between(*EARLY_RECOVERY)]
            model_years = list(range(RECOVERY[0], RECOVERY[1] + 1))
            early_years = list(range(EARLY_RECOVERY[0], EARLY_RECOVERY[1] + 1))
            matched = obs[obs.year.between(1992, 2024)]
            fit_rows.append({
                "treatment": arm, "multiplier": int(arm[1:4]) / 100, "starvation": arm.split("_")[1],
                "reef_name": reef,
                "observed_slope_2017_2024": slope(late.year, late.observed_fast),
                "model_slope_2017_2024": slope(model_years, model.loc[model_years, "fast_coral"]),
                "observed_slope_2001_2013": slope(early.year, early.observed_fast),
                "model_slope_2001_2013": slope(early_years, model.loc[early_years, "fast_coral"]),
                "fast_rmse_vs_photo": float(np.sqrt(np.mean(
                    (model.loc[matched.year, "fast_coral"].to_numpy() - matched.observed_fast.to_numpy()) ** 2))),
                "model_mean_adults_2017_2024": float(model.loc[model_years, "adults_ha"].mean()),
            })
    fit = pd.DataFrame(fit_rows)
    fit["slope_error_sq"] = (fit.model_slope_2017_2024 - fit.observed_slope_2017_2024) ** 2
    fit["early_slope_error_sq"] = (fit.model_slope_2001_2013 - fit.observed_slope_2001_2013) ** 2
    fit.to_csv(paths["growth_starvation_coral_fit.csv"], index=False)
    criterion = fit[fit.starvation == "ref"].groupby("multiplier").slope_error_sq.sum()
    best = criterion.min()
    m_star = float(min(criterion[criterion == best].index))
    boundary = m_star in (min(criterion.index), max(criterion.index))
    selected = [f"m{round(100 * m_star)}_{s}" for s in STARVATION]

    # Stage 2: COTS metrics ----------------------------------------------------------
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
                          (mechanisms.reef_name == score.reef_name)].sort_values("year")
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
        drop_index = int(np.argmin(steps))
        onset = mech[mech.year.between(*ONSET)]
        manta_reef = manta[(manta.reef_name == REEFS[score.reef_name]) &
                           manta.report_year.between(1992, 2024)].groupby("report_year")["median"].mean()
        site = mech.set_index("year").total_coral_site_fraction
        rows.append({
            "treatment": score.treatment, "multiplier": int(score.treatment[1:4]) / 100,
            "starvation": score.treatment.split("_")[1], "reef_name": score.reef_name,
            "first_peak_year": years[0] if years else float("nan"),
            "second_peak_year": years[1] if len(years) > 1 else float("nan"),
            "first_peak_cpue": heights[0] if heights else float("nan"),
            "second_peak_cpue": heights[1] if len(heights) > 1 else float("nan"),
            "max_cpue_2010_2017": float(series[series.year.between(*SECOND_WINDOW)].simulated_cpue.max()),
            "observed_max_cpue_2010_2017": float(obs[obs.year.between(*SECOND_WINDOW)].cotsptow.max()),
            "trough_to_smaller_peak": ratio,
            "matched_peaks": int(score.matched_peaks),
            "smoothed_loss": float(score.loss),
            "steepest_one_year_drop_ratio_1995_2006": float(steps[drop_index]),
            "steepest_drop_year": int(first_wave.year[drop_index + 1]),
            "min_smoothed_adults_1999_2012": float(mech[mech.year.between(1999, 2012)].adults_ha.rolling(3, center=True, min_periods=1).mean().min()),
            "mean_gap_adults_2004_2010": float(mech[mech.year.between(2004, 2010)].adults_ha.mean()),
            "peak_fast_coral_2008_2013": float(onset.fast_coral.max()),
            "max_condition_2008_2013": float(onset.body_condition.max()),
            "coral_site_rmse_vs_ltmp": float(np.sqrt(np.mean((site.loc[manta_reef.index] - manta_reef) ** 2))),
        })
    result = pd.DataFrame(rows)
    result["order"] = result.treatment.map({a: i for i, a in enumerate(ARMS)})
    result = result.sort_values(["order", "reef_name"]).drop(columns="order")
    result.to_csv(paths["growth_starvation_audit.csv"], index=False)
    two_wave = result[result.reef_name.isin(TWO_WAVE)]
    summary = result.groupby("treatment").agg(
        summed_loss=("smoothed_loss", "sum"),
        reefs_with_two_peaks=("second_peak_year", "count"),
        matched_observed_peaks_of_seven=("matched_peaks", "sum"),
        mean_first_peak_year=("first_peak_year", "mean"),
        mean_first_peak_cpue=("first_peak_cpue", "mean"),
        mean_steepest_drop_ratio=("steepest_one_year_drop_ratio_1995_2006", "mean"),
        mean_min_trough_adults_ha=("min_smoothed_adults_1999_2012", "mean"),
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
    coral_summary = fit.groupby("treatment").agg(
        slope_error_sq_2017_2024=("slope_error_sq", "sum"),
        slope_error_sq_2001_2013=("early_slope_error_sq", "sum"),
        mean_fast_rmse_vs_photo=("fast_rmse_vs_photo", "mean"))
    summary = summary.join(coral_summary)

    # Plots ------------------------------------------------------------------------
    cmap = plt.get_cmap("viridis")
    fig, axes = plt.subplots(3, 1, figsize=(9, 10), sharex=True)
    for axis, (reef, photo_reef) in zip(axes, PHOTO_REEFS.items()):
        for i, m in enumerate(MULTIPLIERS):
            arm = f"m{m}_ref"
            model = mechanisms[(mechanisms.treatment == arm) & (mechanisms.reef_name == reef)].sort_values("year")
            axis.plot(model.year, model.fast_coral, color=cmap(i / (len(MULTIPLIERS) - 1)),
                      linewidth=2.4 if m / 100 == m_star else 1.1, label=f"x{m / 100:.2f}" +
                      (" (selected)" if m / 100 == m_star else ""))
        obs = observed[observed.model_reef == reef]
        axis.scatter(obs.year, obs.observed_fast, color="#111827", s=18, zorder=5,
                     label="LTMP photo transect Acropora + Pocilloporidae")
        axis.axvspan(*RECOVERY, color="#16a34a", alpha=.08)
        axis.set_ylabel(f"{reef}\nfast coral (habitable)")
        axis.set_ylim(bottom=0)
        axis.grid(alpha=.15)
    axes[0].set_title("Fast coral, ref starvation arms by fast linear-extension multiplier\n(shaded: 2017-2024 fit window)")
    axes[0].legend(fontsize=7, ncol=2)
    axes[-1].set_xlabel("Year")
    fig.tight_layout()
    fig.savefig(paths["fast_coral_growth_fit.png"], dpi=170)
    plt.close(fig)

    fig, axes = plt.subplots(4, 2, figsize=(14, 13), sharex=True)
    for index, (reef, observed_name) in enumerate(REEFS.items()):
        cots_ax, fast_ax = axes[index]
        for arm in selected + ["m100_ref"]:
            starvation = arm.split("_")[1]
            values = trajectories[(trajectories.treatment == arm) &
                                  (trajectories.reef_name == reef)].sort_values("year")
            mech = mechanisms[(mechanisms.treatment == arm) & (mechanisms.reef_name == reef)].sort_values("year")
            style = ":" if arm == "m100_ref" and arm not in selected else "-"
            label = (f"x{int(arm[1:4]) / 100:.2f} / {starvation}") if index == 0 else None
            cots_ax.plot(values.year.iloc[1:], values.simulated_cpue.iloc[1:], color=COLORS[starvation],
                         linestyle=style, linewidth=1.5, label=label)
            fast_ax.plot(mech.year, mech.fast_coral, color=COLORS[starvation], linestyle=style, linewidth=1.5)
        obs_cots = cots_obs[(cots_obs.reef_name == observed_name) & cots_obs.year.between(1992, 2024)]
        cots_ax.scatter(obs_cots.year, obs_cots.cotsptow, s=16, color="#111827", alpha=.7, zorder=5,
                        label="LTMP COTS/tow" if index == 0 else None)
        if reef in PHOTO_REEFS:
            obs = observed[observed.model_reef == reef]
            fast_ax.scatter(obs.year, obs.observed_fast, s=16, color="#111827", alpha=.7, zorder=5)
        cots_ax.set_ylabel(f"{reef}\nCOTS/tow")
        fast_ax.set_ylabel("Fast coral (habitable)")
        cots_ax.axhline(.22, color="#92400e", linewidth=.8, linestyle=":")
        fast_ax.axhline(.05, color="#92400e", linewidth=.8, linestyle=":")
        for axis in (cots_ax, fast_ax):
            axis.set_xlim(1991, 2024)
            axis.set_ylim(bottom=0)
            axis.grid(alpha=.15)
    axes[0, 0].set_title(f"COTS at q = 0.015; selected fast multiplier x{m_star:.2f}")
    axes[0, 1].set_title("Fast coral / LTMP photo transects (Acropora + Pocilloporidae)")
    axes[-1, 0].set_xlabel("Year")
    axes[-1, 1].set_xlabel("Year")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=5, fontsize=8, bbox_to_anchor=(.5, -.01))
    fig.tight_layout(rect=(0, .04, 1, 1))
    fig.savefig(paths["cots_coral_growth_starvation.png"], dpi=170, bbox_inches="tight")
    plt.close(fig)

    report = {
        "status": "bounded_diagnostic_not_promoted",
        "reference_replay_max_rel_delta": replay,
        "external_immigration_max_abs_delta_vs_reference": external,
        "coral_criterion_by_multiplier": json.loads(criterion.to_json()),
        "selected_multiplier": m_star, "selected_at_grid_boundary": bool(boundary),
        "selected_arms": selected,
        "treatments": json.loads(summary.to_json(orient="index")),
        "criteria": "2 peaks on each of Lizard/MacGillivray/North Direction; 7/7 matched; "
                    "worst trough ratio <=0.10; absolute CPUE under q=0.015",
        "inputs": {name: sha256(run / name) for name in (
            "metadata.toml", "scores.csv", "trajectories.csv", "flows.csv", "mechanisms.csv")} |
                  {"reef_photo_transect.csv": sha256(PHOTO)},
        "script_sha256": sha256(Path(__file__)),
        "outputs": {name: sha256(path) for name, path in paths.items()
                    if name != "growth_starvation_audit.json"},
    }
    paths["growth_starvation_audit.json"].write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    pd.set_option("display.width", 250)
    print("Coral criterion (ref arms):", criterion.round(8).to_dict(), "selected", m_star,
          "boundary" if boundary else "")
    print(fit[fit.starvation == "ref"][["multiplier", "reef_name", "observed_slope_2017_2024",
                                         "model_slope_2017_2024", "observed_slope_2001_2013",
                                         "model_slope_2001_2013", "fast_rmse_vs_photo",
                                         "model_mean_adults_2017_2024"]].round(4).to_string())
    print(summary.round(3).to_string())
    print("Reference replay max delta:", max(replay.values()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    main(parser.parse_args().run)
