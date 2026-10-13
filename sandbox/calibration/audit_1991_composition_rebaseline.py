"""Audit the preregistered first-wave composition range on the re-baselined model
(crash_mechanism_protocol.md, 2026-10-13): 1991 fast share {acrpoc, halfother, other}
x bleaching beta {1, 1.5} on the working coral configuration with Q = 20. The acrpoc
arms must replay run 11; the adoption rule is frozen in the protocol."""

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
from audit_1991_coral_calibration import observed_photo  # noqa: E402
from audit_1991_growth_starvation import PHOTO, PHOTO_REEFS  # noqa: E402
from audit_1991_initiation import exact_replay  # noqa: E402
from audit_1991_joint_coral import objective, reef_errors  # noqa: E402
from audit_1991_switching_repair import REEFS, TWO_WAVE, sha256  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
PRIOR = ROOT / "sandbox/calibration/runs/20261013T_11_cots_rebaseline_1991_v1"
BETA, LEVELS = (100, 150), ("acrpoc", "halfother", "other")
ARMS = [f"wc_b{b}_q20_{level}" for b in BETA for level in LEVELS]
REFERENCES = {f"wc_b{b}_q20_acrpoc": f"wc_b{b}_q20" for b in BETA}
COLORS = {"acrpoc": "#dc2626", "halfother": "#d97706", "other": "#0284c7"}
N_YEARS = 34


def main(run: Path) -> None:
    paths = {name: run / name for name in (
        "composition_rebaseline_audit.csv", "composition_rebaseline_audit.json",
        "composition_rebaseline.png")}
    if any(path.exists() for path in paths.values()):
        raise FileExistsError("Refusing to overwrite composition audit")
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
        for arm, prior_arm in REFERENCES.items():
            replay[f"{arm}/{name}"] = exact_replay(tables[name][tables[name].treatment == arm],
                                                   prior[prior.treatment == prior_arm], columns)
    reference = flows[flows.treatment == ARMS[0]].sort_values(["reef_name", "year"])
    external = {arm: float(np.max(np.abs(
        flows[flows.treatment == arm].sort_values(["reef_name", "year"]).external_immigration.to_numpy() -
        reference.external_immigration.to_numpy()))) for arm in ARMS}
    if max(external.values()) > 1e-9:
        raise ValueError("External immigration differs between arms")
    numeric = mechanisms[["recruits_ha", "subadults_ha", "adults_ha", "fast_coral", "total_coral"]].to_numpy()
    if not np.isfinite(numeric).all() or (numeric < -1e-12).any():
        raise ValueError("Nonfinite or negative state")
    start = mechanisms[mechanisms.year == 1991].pivot(index="reef_name", columns="treatment", values="total_coral")
    total_delta = float(start[ARMS].sub(start[ARMS[0]], axis=0).abs().to_numpy().max())
    if total_delta > 1e-9:
        raise ValueError(f"Initial total cover differs between compositions: {total_delta}")

    photo = observed_photo()
    cots_obs = pd.read_csv(ROOT / "sandbox/data/reef_cots.csv")
    manta_raw = pd.read_csv(ROOT / "sandbox/data/reef_manta.csv")
    manta_raw = manta_raw[(manta_raw.data_type == "manta") & (manta_raw.domain_category == "reef") &
                          (manta_raw.variable == "HC") & (manta_raw.purpose == "MANTA") &
                          (manta_raw.project_code == "LTMP") & (manta_raw.depth == 9)]
    rows = []
    for score in tables["scores"][tables["scores"].observation_treatment == "smoothed_3y"].itertuples(index=False):
        series = trajectories[(trajectories.treatment == score.treatment) &
                              (trajectories.reef_name == score.reef_name) &
                              (trajectories.year >= 1992)].sort_values("year")
        mech = mechanisms[(mechanisms.treatment == score.treatment) &
                          (mechanisms.reef_name == score.reef_name)].sort_values("year").set_index("year")
        mech = mech.assign(slow_coral=mech.total_coral - mech.fast_coral)
        obs = cots_obs[cots_obs.reef_name == REEFS[score.reef_name]]
        years = [int(v) for v in str(score.sim_peak_years).split(";") if v and v != "nan"]
        heights = [float(v) for v in str(score.sim_peak_heights_cpue).split(";") if v and v != "nan"]
        ratio = float("nan")
        if len(years) >= 2:
            smooth = series.simulated_cpue.rolling(3, center=True, min_periods=1).mean()
            ratio = float(smooth[series.year.between(years[0], years[1])].min() / min(heights[:2]))
        manta = manta_raw[(manta_raw.reef_name == REEFS[score.reef_name]) &
                          manta_raw.report_year.between(1992, 2024)].groupby("report_year")["median"].mean()
        _, b, _, level = score.treatment.split("_")
        row = {
            "treatment": score.treatment, "reef_name": score.reef_name, "composition": level,
            "bleaching_scalar": int(b[1:]) / 100,
            "first_peak_year": years[0] if years else float("nan"),
            "first_peak_cpue": heights[0] if heights else float("nan"),
            "second_peak_year": years[1] if len(years) > 1 else float("nan"),
            "second_peak_cpue": heights[1] if len(heights) > 1 else float("nan"),
            "max_cpue_1995_2000": float(series[series.year.between(1995, 2000)].simulated_cpue.max()),
            "observed_max_cpue_1995_2000": float(obs[obs.year.between(1995, 2000)].cotsptow.max()),
            "max_cpue_2010_2017": float(series[series.year.between(2010, 2017)].simulated_cpue.max()),
            "observed_max_cpue_2010_2017": float(obs[obs.year.between(2010, 2017)].cotsptow.max()),
            "trough_to_smaller_peak": ratio, "matched_peaks": int(score.matched_peaks),
            "smoothed_loss": float(score.loss),
            "min_smoothed_adults_1999_2012": float(mech.loc[1999:2012].adults_ha.rolling(3, center=True, min_periods=1).mean().min()),
            "fast_share_1991": float(mech.loc[1991, "fast_coral"] / mech.loc[1991, "total_coral"]),
            "mean_fast_coral_1992_1996": float(mech.loc[1992:1996].fast_coral.mean()),
            "max_fast_coral_1992_1996": float(mech.loc[1992:1996].fast_coral.max()),
            "mean_fast_coral_2000_2013": float(mech.loc[2000:2013].fast_coral.mean()),
            "max_subadults_ha": float(mech.loc[1992:2024].subadults_ha.max()),
        }
        row |= reef_errors(mech, score.reef_name, photo, manta, (1992, 2024))
        rows.append(row)
    result = pd.DataFrame(rows)
    result["order"] = result.treatment.map({a: i for i, a in enumerate(ARMS)})
    result = result.sort_values(["order", "reef_name"]).drop(columns="order")
    result.to_csv(paths["composition_rebaseline_audit.csv"], index=False)

    summary = result.groupby("treatment").agg(
        bleaching_scalar=("bleaching_scalar", "first"), composition=("composition", "first"),
        summed_loss=("smoothed_loss", "sum"),
        matched_observed_peaks_of_seven=("matched_peaks", "sum"),
        mean_first_peak_year=("first_peak_year", "mean"),
        mean_first_peak_cpue=("first_peak_cpue", "mean"),
        mean_min_trough_adults_ha=("min_smoothed_adults_1999_2012", "mean"),
        mean_rmse_fast=("rmse_fast_vs_photo", "mean"),
    ).reindex(ARMS)
    for reef, short in (("Lizard Island Reef", "liz"), ("MacGillivray Reef", "macg"),
                        ("North Direction Reef", "nd"), ("Eyrie Reef", "eyr")):
        r = result[result.reef_name == reef].set_index("treatment")
        summary[f"{short}_first_peak"] = r.first_peak_cpue
    gate = result[result.reef_name.isin(TWO_WAVE)].groupby("treatment").agg(
        two_wave_reefs_with_two_peaks=("second_peak_year", "count"),
        worst_available_trough_ratio=("trough_to_smaller_peak", "max"))
    summary = summary.join(gate)
    summary["passes_qualitative_gate"] = (
        (summary.two_wave_reefs_with_two_peaks == len(TWO_WAVE)) &
        (summary.matched_observed_peaks_of_seven == 7) &
        (summary.worst_available_trough_ratio <= .10))
    summary["coral_J"] = objective(result).reindex(ARMS)
    adoption = {}
    for b in BETA:
        base = summary.loc[f"wc_b{b}_q20_acrpoc"]
        for level in LEVELS[1:]:
            arm = summary.loc[f"wc_b{b}_q20_{level}"]
            checks = {
                "lizard_and_macgillivray_first_peaks_rise": bool(arm.liz_first_peak > base.liz_first_peak and
                                                                  arm.macg_first_peak > base.macg_first_peak),
                "mean_trough_le_4": bool(arm.mean_min_trough_adults_ha <= 4.0),
                "loss_worsens_le_0.5": bool(arm.summed_loss - base.summed_loss <= 0.5),
                "fast_rmse_worsens_le_0.01": bool(arm.mean_rmse_fast - base.mean_rmse_fast <= 0.01),
            }
            adoption[f"wc_b{b}_q20_{level}"] = checks | {"meets_rule": all(checks.values())}

    fig, axes = plt.subplots(4, 2, figsize=(14, 13), sharex=True)
    for index, (reef, observed_name) in enumerate(REEFS.items()):
        cots_ax, fast_ax = axes[index]
        for arm in ARMS:
            _, b, _, level = arm.split("_")
            values = trajectories[(trajectories.treatment == arm) & (trajectories.reef_name == reef)].sort_values("year")
            mech = mechanisms[(mechanisms.treatment == arm) & (mechanisms.reef_name == reef)].sort_values("year")
            kwargs = dict(color=COLORS[level], linestyle="-" if b == "b150" else ":", linewidth=1.6,
                          label=f"{level} / beta {int(b[1:]) / 100}" if index == 0 else None)
            cots_ax.plot(values.year.iloc[1:], values.simulated_cpue.iloc[1:], **kwargs)
            fast_ax.plot(mech.year, mech.fast_coral, **kwargs)
        obs_cots = cots_obs[(cots_obs.reef_name == observed_name) & cots_obs.year.between(1992, 2024)]
        cots_ax.scatter(obs_cots.year, obs_cots.cotsptow, s=16, color="#111827", alpha=.7, zorder=5,
                        label="LTMP observations" if index == 0 else None)
        if reef in PHOTO_REEFS:
            p = photo[photo.model_reef == reef]
            fast_ax.scatter(p.year, p.obs_fast, s=16, color="#111827", alpha=.7, zorder=5)
        cots_ax.set_ylabel(f"{reef}\nCOTS/tow (q = 0.015)")
        fast_ax.set_ylabel("fast coral (habitable)")
        fast_ax.axhline(.05, color="#92400e", linewidth=.8, linestyle=":")
        for axis in (cots_ax, fast_ax):
            axis.set_xlim(1991, 2024)
            axis.set_ylim(bottom=0)
            axis.grid(alpha=.15)
    axes[0, 0].set_title("COTS: 1991 fast share x bleaching beta (solid: 1.5), working coral, Q 20")
    axes[0, 1].set_title("Fast coral vs photo transects (dotted: theta 0.05)")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, fontsize=8, bbox_to_anchor=(.5, -.01))
    fig.tight_layout(rect=(0, .04, 1, 1))
    fig.savefig(paths["composition_rebaseline.png"], dpi=160, bbox_inches="tight")
    plt.close(fig)

    report = {
        "status": "bounded_diagnostic_not_promoted",
        "reference_replay_max_rel_delta": replay,
        "external_immigration_max_abs_delta_vs_first_arm": external,
        "initial_total_cover_max_abs_delta": total_delta,
        "adoption_rule_checks": adoption,
        "treatments": json.loads(summary.to_json(orient="index")),
        "inputs": {name: sha256(run / name) for name in (
            "metadata.toml", "scores.csv", "trajectories.csv", "flows.csv", "mechanisms.csv")} |
                  {"reef_photo_transect.csv": sha256(PHOTO)},
        "script_sha256": sha256(Path(__file__)),
        "outputs": {name: sha256(path) for name, path in paths.items()
                    if name != "composition_rebaseline_audit.json"},
    }
    paths["composition_rebaseline_audit.json"].write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 40)
    print(summary.round(3).to_string())
    print("Adoption checks:", json.dumps(adoption, indent=1))
    print("Reference replay max delta:", max(replay.values()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    main(parser.parse_args().run.resolve())
