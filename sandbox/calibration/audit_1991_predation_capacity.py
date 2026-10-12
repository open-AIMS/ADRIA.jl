"""Audit the preregistered predator-pit capacity tuning grid on the v0.3 domain:
External x5 (working assumption) + gradual Allee, with P off or Q on a fixed
log-spaced grid (crash_mechanism_protocol.md, 2026-10-12). Q is a tuned parameter."""

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
PRIOR = ROOT / "sandbox/calibration/runs/20261012T_initiation_1991_v3"
CAPACITIES = (20, 30, 45, 70, 100, 150, 220, 330, 500)
ARMS = ["ext5_A"] + [f"ext5_PA_q{q}" for q in CAPACITIES]
REFERENCES = {"ext5_A": "ext5_A", "ext5_PA_q20": "ext5_PA"}
GAP = (2004, 2010)
SECOND_WINDOW = (2010, 2017)
N_YEARS = 34
PLATEAU = 0.10


def capacity(arm: str) -> float:
    return float(arm.split("_q")[1]) if "_q" in arm else 0.0


def select(losses: pd.DataFrame, reefs: list[str]) -> str:
    """Frozen rule: lowest summed loss over `reefs` among arms with a first peak
    on every reef; ties to the smaller Q."""
    subset = losses[losses.reef_name.isin(reefs)]
    eligible = subset.groupby("treatment").first_peak_year.apply(lambda v: v.notna().all())
    total = subset.groupby("treatment").smoothed_loss.sum()[eligible[eligible].index]
    best = total.min()
    return min(total[total == best].index, key=capacity)


def main(run: Path) -> None:
    paths = {name: run / name for name in (
        "predation_capacity_audit.csv", "predation_capacity_audit.json",
        "cots_coral_predation_capacity.png", "predation_capacity_loss.png")}
    if any(path.exists() for path in paths.values()):
        raise FileExistsError("Refusing to overwrite predation-capacity audit")
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
    for arm, prior_arm in REFERENCES.items():
        for name in ("trajectories", "flows", "mechanisms"):
            prior = pd.read_csv(PRIOR / f"{name}.csv")
            columns = [c for c in prior if c not in ("treatment", "reef_name", "year")]
            replay[f"{arm}/{name}"] = exact_replay(tables[name][tables[name].treatment == arm],
                                                   prior[prior.treatment == prior_arm], columns)
    prior_flows = pd.read_csv(PRIOR / "flows.csv")
    reference = prior_flows[prior_flows.treatment == "ext5_A"].sort_values(["reef_name", "year"])
    external = {}
    for arm in ARMS:
        current = flows[flows.treatment == arm].sort_values(["reef_name", "year"])
        external[arm] = float(np.max(np.abs(current.external_immigration.to_numpy() -
                                            reference.external_immigration.to_numpy())))
        if external[arm] > 1e-9:
            raise ValueError(f"External immigration differs from Ext x5 in {arm}")
    numeric = mechanisms[mechanisms.year >= 1992][[
        "recruits_ha", "subadults_ha", "adults_ha", "fast_coral", "total_coral",
        "food_survival", "body_condition", "juvenile_predation_ha"]].to_numpy()
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
        years = [int(v) for v in str(score.sim_peak_years).split(";") if v and v != "nan"]
        heights = [float(v) for v in str(score.sim_peak_heights_cpue).split(";") if v and v != "nan"]
        ratio = float("nan")
        if len(years) >= 2:
            smooth = series.simulated_cpue.rolling(3, center=True, min_periods=1).mean()
            ratio = float(smooth[series.year.between(years[0], years[1])].min() / min(heights[:2]))
        gap = mech[mech.year.between(*GAP)]
        gap_flows = reef_flows[reef_flows.year.between(*GAP)]
        rows.append({
            "treatment": score.treatment, "Q": capacity(score.treatment),
            "reef_name": score.reef_name,
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
            "min_smoothed_adults_1999_2012": float(mech[mech.year.between(1999, 2012)].adults_ha.rolling(3, center=True, min_periods=1).mean().min()),
            "mean_gap_settlement_ha": float(gap_flows.settled_recruits.mean()),
            "mean_gap_juvenile_predation_ha": gap.juvenile_predation_ha.mean(),
            "gap_external_share": float(gap_flows.external_immigration.sum() /
                                        max(gap_flows.settled_recruits.sum(), 1e-12)),
            "cpue_1994": float(series[series.year == 1994].simulated_cpue.iloc[0]),
            "cpue_1998": float(series[series.year == 1998].simulated_cpue.iloc[0]),
            "model_coral_site_fraction_2012": float(mech[mech.year == 2012].total_coral_site_fraction.iloc[0]),
            "model_coral_site_fraction_2016": float(mech[mech.year == 2016].total_coral_site_fraction.iloc[0]),
            "observed_manta_coral_2010_2014": float(manta[(manta.reef_name == REEFS[score.reef_name]) &
                manta.report_year.between(2010, 2014)]["median"].mean()),
        })
    result = pd.DataFrame(rows).sort_values(["Q", "reef_name"])
    result.to_csv(paths["predation_capacity_audit.csv"], index=False)

    two_wave = result[result.reef_name.isin(TWO_WAVE)]
    summary = result.groupby("treatment").agg(
        Q=("Q", "first"),
        summed_loss=("smoothed_loss", "sum"),
        reefs_with_first_peak=("first_peak_year", "count"),
        reefs_with_two_peaks=("second_peak_year", "count"),
        matched_observed_peaks_of_seven=("matched_peaks", "sum"),
        mean_first_peak_year=("first_peak_year", "mean"),
        mean_max_cpue_2010_2017=("max_cpue_2010_2017", "mean"),
        mean_min_trough_adults_ha=("min_smoothed_adults_1999_2012", "mean"),
        mean_gap_adults_ha=("mean_gap_adults_ha", "mean"),
        mean_gap_settlement_ha=("mean_gap_settlement_ha", "mean"),
        mean_gap_juvenile_predation_ha=("mean_gap_juvenile_predation_ha", "mean"),
        mean_coral_site_2012=("model_coral_site_fraction_2012", "mean"),
    ).sort_values("Q")
    gate = two_wave.groupby("treatment").agg(
        two_wave_reefs_with_two_peaks=("second_peak_year", "count"),
        worst_available_trough_ratio=("trough_to_smaller_peak", "max"))
    summary = summary.join(gate)
    summary["passes_qualitative_gate"] = (
        (summary.two_wave_reefs_with_two_peaks == len(TWO_WAVE)) &
        (summary.matched_observed_peaks_of_seven == 7) &
        (summary.worst_available_trough_ratio <= .10))

    reefs = list(REEFS)
    selected = select(result, reefs)
    best_loss = summary.loc[selected, "summed_loss"]
    plateau = summary[summary.summed_loss <= (1 + PLATEAU) * best_loss].Q.tolist()
    held_out = {}
    for reef in reefs:
        chosen = select(result, [r for r in reefs if r != reef])
        row = result[(result.treatment == chosen) & (result.reef_name == reef)].iloc[0]
        ref = result[(result.treatment == "ext5_A") & (result.reef_name == reef)].iloc[0]
        held_out[reef] = {
            "selected_on_other_reefs": chosen, "Q": capacity(chosen),
            "held_out_loss": float(row.smoothed_loss), "reference_ext5_A_loss": float(ref.smoothed_loss),
            "held_out_peak_years": [int(y) for y in (row.first_peak_year, row.second_peak_year) if y == y],
            "held_out_matched_peaks": int(row.matched_peaks),
        }

    cmap = plt.get_cmap("viridis")
    colors = {arm: ("#475569" if arm == "ext5_A" else cmap(i / (len(CAPACITIES) - 1)))
              for i, arm in enumerate(ARMS[1:])} | {"ext5_A": "#475569"}
    fig, axes = plt.subplots(4, 2, figsize=(14, 13), sharex=True)
    for index, (reef, observed) in enumerate(REEFS.items()):
        cots_ax, coral_ax = axes[index]
        for arm in ARMS:
            values = trajectories[(trajectories.treatment == arm) &
                                  (trajectories.reef_name == reef)].sort_values("year")
            fast = mechanisms[(mechanisms.treatment == arm) &
                              (mechanisms.reef_name == reef)].sort_values("year")
            style = "--" if arm == "ext5_A" else "-"
            width = 2.4 if arm == selected else 1.1
            label = ("A only (no P)" if arm == "ext5_A" else f"Q={capacity(arm):g}") + \
                (" (selected)" if arm == selected else "")
            cots_ax.plot(values.year.iloc[1:], values.simulated_cpue.iloc[1:], color=colors[arm],
                         linestyle=style, linewidth=width, label=label if index == 0 else None)
            coral_ax.plot(fast.year, fast.total_coral_site_fraction, color=colors[arm],
                          linestyle=style, linewidth=width)
        obs_cots = observations[(observations.reef_name == observed) &
                                observations.year.between(1992, 2024)]
        cots_ax.scatter(obs_cots.year, obs_cots.cotsptow, s=16, color="#111827", alpha=.7, zorder=5,
                        label="LTMP COTS/tow" if index == 0 else None)
        obs_coral = manta[(manta.reef_name == observed) & manta.report_year.between(1992, 2024)]
        coral_ax.scatter(obs_coral.report_year, obs_coral["median"], s=16, color="#111827", alpha=.6, zorder=5)
        cots_ax.set_ylabel(f"{reef}\nCOTS/tow")
        coral_ax.set_ylabel("Coral fraction of site area")
        cots_ax.axhline(.22, color="#92400e", linewidth=.8, linestyle=":")
        for axis in (cots_ax, coral_ax):
            axis.set_xlim(1991, 2024)
            axis.set_ylim(bottom=0)
            axis.grid(alpha=.15)
    axes[0, 0].set_title("COTS at q = 0.015, External x5 + gradual Allee, predator-pit capacity Q")
    axes[0, 1].set_title("Model coral (site-area fraction) / LTMP 9 m manta")
    axes[-1, 0].set_xlabel("Year")
    axes[-1, 1].set_xlabel("Year")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=6, fontsize=8, bbox_to_anchor=(.5, -.01))
    fig.tight_layout(rect=(0, .05, 1, 1))
    fig.savefig(paths["cots_coral_predation_capacity.png"], dpi=170, bbox_inches="tight")
    plt.close(fig)

    pits = summary[summary.Q > 0]
    fig, (loss_ax, trough_ax) = plt.subplots(1, 2, figsize=(11, 4))
    for reef in reefs:
        per = result[(result.reef_name == reef) & (result.Q > 0)]
        loss_ax.plot(per.Q, per.smoothed_loss, marker="o", linewidth=1, label=reef)
        trough_ax.plot(per.Q, per.min_smoothed_adults_1999_2012, marker="o", linewidth=1, label=reef)
    loss_ax.plot(pits.Q, pits.summed_loss / len(reefs), color="black", linewidth=2.2, label="mean")
    loss_ax.axhline(summary.loc["ext5_A", "summed_loss"] / len(reefs), color="#475569",
                    linestyle="--", label="A only (mean)")
    trough_ax.axhspan(0, 3, color="#16a34a", alpha=.12, label="observed trough ~0-3")
    for axis, ylabel in ((loss_ax, "Frozen smoothed loss"), (trough_ax, "Trough min adults/ha (1999-2012)")):
        axis.set_xscale("log")
        axis.set_xlabel("Q (COTS/ha/yr)")
        axis.set_ylabel(ylabel)
        axis.axvline(capacity(selected), color="#dc2626", linewidth=.8)
        axis.grid(alpha=.2)
    loss_ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(paths["predation_capacity_loss.png"], dpi=170)
    plt.close(fig)

    report = {
        "status": "bounded_tuning_grid_not_promoted",
        "Q_status": "tuned parameter (not evidence-bounded)",
        "reference_replay_max_rel_delta": replay,
        "external_immigration_max_abs_delta_vs_ext5_A": external,
        "selected_arm": selected, "selected_Q": capacity(selected),
        "selection_rule": "min summed smoothed_3y fixed_q0015 loss over 4 reefs; first peak on every reef; ties to smaller Q",
        "plateau_within_10pct_Q": plateau,
        "leave_one_reef_out": held_out,
        "treatments": json.loads(summary.to_json(orient="index")),
        "criteria": "2 peaks on each of Lizard/MacGillivray/North Direction; 7/7 matched; "
                    "worst trough ratio <=0.10; absolute CPUE under q=0.015",
        "inputs": {name: sha256(run / name) for name in (
            "metadata.toml", "scores.csv", "trajectories.csv", "flows.csv", "mechanisms.csv")},
        "script_sha256": sha256(Path(__file__)),
        "outputs": {name: sha256(path) for name, path in paths.items()
                    if name != "predation_capacity_audit.json"},
    }
    paths["predation_capacity_audit.json"].write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    pd.set_option("display.width", 250)
    print(summary.round(3).to_string())
    print("Selected:", selected, "plateau Q:", plateau)
    print("Leave-one-reef-out:", json.dumps(held_out, indent=1))
    print("Reference replay max delta:", max(replay.values()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    main(parser.parse_args().run)
