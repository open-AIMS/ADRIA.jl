"""Audit the preregistered first-wave initiation screen on the adopted v0.3 domain:
external supply x5 and/or 1991 juvenile preload x5, crossed with none / A / P+A
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
from audit_1991_switching_repair import REEFS, TWO_WAVE, sha256  # noqa: E402


def exact_replay(current: pd.DataFrame, old: pd.DataFrame, columns: list[str]) -> float:
    """Relative replay check (1e-9). The v0.3 package stores cyclone mortality
    computed with Python's exp, whereas the prior run computed it with Julia's exp
    at run time; these differ in the last bit, giving ~1e-15 relative deltas."""
    joined = current.merge(old, on=["reef_name", "year"], suffixes=("_new", "_old"),
                           validate="one_to_one")
    if len(joined) != len(current) or len(joined) != len(old):
        raise ValueError("Reference and prior run do not cover the same reef-years")
    worst = 0.0
    for column in columns:
        new, ref = joined[f"{column}_new"], joined[f"{column}_old"]
        both_nan = new.isna() & ref.isna()
        rel = np.abs(new - ref) / np.maximum(1.0, np.abs(ref))
        worst = max(worst, float(np.nanmax(np.where(both_nan, 0.0, rel))))
    if worst > 1e-9:
        raise ValueError(f"Reference arm did not replay prior run: {worst}")
    return worst

ROOT = Path(__file__).resolve().parents[2]
PRIOR = ROOT / "sandbox/calibration/runs/20261012T_trough_mechanisms_1991_v1"
BASE = "none_none"
PRIOR_BASE = "bgeo_none"
EXTERNAL_MULTIPLIER = 5.0
COLORS = {"none": "#475569", "ext5": "#0d9488", "juv5": "#d97706", "both": "#dc2626"}
STYLES = {"none": "-", "A": "--", "PA": ":"}
ARMS = {f"{init}_{mech}": (f"{init} / {mech}", COLORS[init], STYLES[mech])
        for init in ("none", "ext5", "juv5", "both") for mech in ("none", "A", "PA")}
GAP = (2004, 2010)
ONSET = (2008, 2011)
SECOND_WINDOW = (2010, 2017)
N_YEARS = 34


def main(run: Path) -> None:
    paths = {name: run / name for name in (
        "initiation_audit.csv", "initiation_audit.json",
        "cots_coral_initiation.png")}
    if any(path.exists() for path in paths.values()):
        raise FileExistsError("Refusing to overwrite initiation audit")
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
        factor = EXTERNAL_MULTIPLIER if arm.startswith(("ext5", "both")) else 1.0
        external[arm] = float(np.max(np.abs(current.external_immigration.to_numpy() -
                                            factor * reference.external_immigration.to_numpy())))
        if external[arm] > 1e-9:
            raise ValueError(f"External immigration differs from x{factor} in {arm}")
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
        reef_flows = flows[(flows.treatment == score.treatment) &
                           (flows.reef_name == score.reef_name)]
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
            "mean_gap_juvenile_predation_ha": gap.juvenile_predation_ha.mean(),
            "mean_gap_subadults_ha": gap.subadults_ha.mean(),
            "mean_gap_maturation_fraction": gap.maturation_fraction.mean(),
            "min_smoothed_adults_1999_2012": float(mech[mech.year.between(1999, 2012)].adults_ha.rolling(3, center=True, min_periods=1).mean().min()),
            "settlement_1992_1995": float(reef_flows[reef_flows.year.between(1992, 1995)].settled_recruits.mean()),
            "external_share_1992_1995": float(reef_flows[reef_flows.year.between(1992, 1995)].external_immigration.sum() /
                                              max(reef_flows[reef_flows.year.between(1992, 1995)].settled_recruits.sum(), 1e-12)),
            "cpue_1994": float(series[series.year == 1994].simulated_cpue.iloc[0]),
            "cpue_1998": float(series[series.year == 1998].simulated_cpue.iloc[0]),
            "gap_external_share": float(reef_flows[reef_flows.year.between(*GAP)].external_immigration.sum() /
                                        max(reef_flows[reef_flows.year.between(*GAP)].settled_recruits.sum(), 1e-12)),
            "observed_manta_coral_2016_2018": float(manta[(manta.reef_name == REEFS[score.reef_name]) &
                manta.report_year.between(2016, 2018)]["median"].mean()),
            "observed_manta_coral_2010_2014": float(obs_coral["median"].mean()) if len(obs_coral) else float("nan"),
        })
    result = pd.DataFrame(rows)
    result.to_csv(paths["initiation_audit.csv"], index=False)
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
        mean_gap_subadults_ha=("mean_gap_subadults_ha", "mean"),
        mean_min_trough_adults_ha=("min_smoothed_adults_1999_2012", "mean"),
        mean_first_peak_year=("first_peak_year", "mean"),
        mean_cpue_1998=("cpue_1998", "mean"),
        mean_settlement_1992_1995=("settlement_1992_1995", "mean"),
        mean_gap_external_share=("gap_external_share", "mean"),
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
            coral_ax.plot(fast.year, fast.total_coral_site_fraction, color=color, linestyle=style, linewidth=1.4)
        obs_cots = observations[(observations.reef_name == observed) &
                                observations.year.between(1992, 2024)]
        cots_ax.scatter(obs_cots.year, obs_cots.cotsptow, s=16, color="#111827", alpha=.7,
                        label="LTMP COTS/tow" if index == 0 else None)
        obs_coral = manta[(manta.reef_name == observed) & manta.report_year.between(1992, 2024)]
        coral_ax.scatter(obs_coral.report_year, obs_coral["median"], s=16, color="#111827", alpha=.6)
        cots_ax.set_ylabel(f"{reef}\nCOTS/tow")
        coral_ax.set_ylabel("Coral fraction of site area")
        cots_ax.axhline(.22, color="#92400e", linewidth=.8, linestyle=":")
        for axis in (cots_ax, coral_ax):
            axis.set_xlim(1991, 2024)
            axis.set_ylim(bottom=0)
            axis.grid(alpha=.15)
    axes[0, 0].set_title("COTS at fixed q = 0.015 COTS/tow per COTS/ha")
    axes[0, 1].set_title("Model coral (site-area fraction) / LTMP 9 m manta")
    axes[-1, 0].set_xlabel("Year")
    axes[-1, 1].set_xlabel("Year")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, fontsize=8, bbox_to_anchor=(.5, -.01))
    fig.tight_layout(rect=(0, .05, 1, 1))
    fig.savefig(paths["cots_coral_initiation.png"], dpi=170, bbox_inches="tight")
    plt.close(fig)

    report = {
        "status": "bounded_diagnostic_not_promoted",
        "base_replay_of_trough_bgeo_none_max_abs_delta": replay,
        "external_immigration_max_abs_delta_vs_expected_multiple": external,
        "treatments": json.loads(summary.to_json(orient="index")),
        "criteria": "2 peaks on each of Lizard/MacGillivray/North Direction; 7/7 matched; "
                    "worst trough ratio <=0.10; absolute CPUE under q=0.015",
        "inputs": {name: sha256(run / name) for name in (
            "metadata.toml", "scores.csv", "trajectories.csv", "flows.csv", "mechanisms.csv")},
        "script_sha256": sha256(Path(__file__)),
        "outputs": {name: sha256(path) for name, path in paths.items()
                    if name != "initiation_audit.json"},
    }
    paths["initiation_audit.json"].write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    pd.set_option("display.width", 250)
    print(summary.round(3).to_string())
    print("Base replay max delta:", max(replay.values()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    main(parser.parse_args().run)
