"""Audit the preregistered Lizard Island local-retention screen
(crash_mechanism_protocol.md, 2026-10-13): source-conserving retention in the Lizard
complex (rho 0.25 / 0.5 / 0.9) or on every reef (0.5) on the adopted eb/lin baseline.
Gate G2' is evaluated separately with evaluate_gate_g2prime.py on the run's re-score."""

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
from audit_1991_switching_repair import REEFS, sha256  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
PRIOR = ROOT / "sandbox/calibration/runs/20261013T_17_initial_state_rebaseline_1991_v1"
LEVELS = ("owen", "liz25", "liz50", "liz90", "all50")
ARMS = [f"wc_b150_q20_eb_lin_{level}" for level in LEVELS]
REFERENCE, PRIOR_REFERENCE = ARMS[0], "wc_b150_q20_eb_lin"
COLORS = {"owen": "#52514e", "liz25": "#86b6ef", "liz50": "#3987e5", "liz90": "#0d366b", "all50": "#eb6834"}
N_YEARS = 34
Q = 0.015


def main(run: Path) -> None:
    paths = {name: run / name for name in (
        "local_retention_audit.csv", "local_retention_audit.json", "local_retention.png")}
    if any(path.exists() for path in paths.values()):
        raise FileExistsError("Refusing to overwrite local-retention audit")
    tables = {name: pd.read_csv(run / f"{name}.csv")
              for name in ("trajectories", "flows", "mechanisms", "scores", "failures")}
    if len(tables["failures"]):
        raise ValueError("Failed evaluations in run")
    trajectories, flows, mechanisms = (tables[k] for k in ("trajectories", "flows", "mechanisms"))
    if set(trajectories.treatment) != set(ARMS) or len(trajectories) != len(ARMS) * 4 * N_YEARS:
        raise ValueError("Missing, extra or incomplete treatments")
    replay = {}
    for name in ("trajectories", "flows", "mechanisms"):
        prior = pd.read_csv(PRIOR / f"{name}.csv")
        columns = [c for c in prior if c not in ("treatment", "reef_name", "year")]
        replay[name] = exact_replay(tables[name][tables[name].treatment == REFERENCE],
                                    prior[prior.treatment == PRIOR_REFERENCE], columns)
    reference = flows[flows.treatment == REFERENCE].sort_values(["reef_name", "year"])
    external = {arm: float(np.max(np.abs(flows[flows.treatment == arm].sort_values(["reef_name", "year"])
                                         .external_immigration.to_numpy() - reference.external_immigration.to_numpy())))
                for arm in ARMS}
    if max(external.values()) > 1e-9:
        raise ValueError("External immigration changed: retention must only redistribute internal larvae")
    if not np.isfinite(mechanisms[["recruits_ha", "subadults_ha", "adults_ha"]].to_numpy()).all():
        raise ValueError("Nonfinite state")
    metadata = (run / "metadata.toml").read_text(encoding="utf-8")
    cots_obs = pd.read_csv(ROOT / "sandbox/data/reef_cots.csv")

    rows = []
    smooth = tables["scores"][tables["scores"].observation_treatment == "smoothed_3y"]
    for s in smooth.itertuples(index=False):
        traj = trajectories[(trajectories.treatment == s.treatment) & (trajectories.reef_name == s.reef_name)].set_index("year")
        flow = flows[(flows.treatment == s.treatment) & (flows.reef_name == s.reef_name)].set_index("year")
        years = [int(v) for v in str(s.sim_peak_years).split(";") if v and v != "nan"]
        heights = [float(v) for v in str(s.sim_peak_heights_cpue).split(";") if v and v != "nan"]
        rows.append({
            "treatment": s.treatment, "retention": s.treatment.split("_")[-1], "reef_name": s.reef_name,
            "smoothed_loss": float(s.loss), "matched_smoothed": int(s.matched_peaks),
            "first_peak_year": years[0] if years else np.nan, "first_peak_cpue": heights[0] if heights else np.nan,
            "second_peak_year": years[1] if len(years) > 1 else np.nan,
            "second_peak_cpue": heights[1] if len(heights) > 1 else np.nan,
            "max_cpue_1995_2000": float(traj.loc[1995:2000].simulated_cpue.max()),
            "max_cpue_2010_2017": float(traj.loc[2010:2017].simulated_cpue.max()),
            "mean_cpue_1999_2011": float(traj.loc[1999:2011].simulated_cpue.mean()),
            "min_adults_ha_1999_2011": float(traj.loc[1999:2011].adults_ha.min()),
            "local_retention_1992_2015": float(flow.loc[1992:2015].local_retention.sum()),
            "internal_immigration_1992_2015": float(flow.loc[1992:2015].internal_immigration.sum()),
            "external_immigration_1992_2015": float(flow.loc[1992:2015].external_immigration.sum()),
            "local_fecundity_1992_2015": float(flow.loc[1992:2015].local_fecundity.sum()),
        })
    result = pd.DataFrame(rows)
    result["order"] = result.treatment.map({a: i for i, a in enumerate(ARMS)})
    result = result.sort_values(["order", "reef_name"]).drop(columns="order")
    result.to_csv(paths["local_retention_audit.csv"], index=False)

    raw = tables["scores"][tables["scores"].observation_treatment == "raw"].groupby("treatment").matched_peaks.sum()
    summary = result.groupby("treatment").agg(summed_loss=("smoothed_loss", "sum"),
                                              matched_smoothed=("matched_smoothed", "sum")).reindex(ARMS)
    summary["matched_raw"] = raw.reindex(ARMS)
    for reef, short in (("Lizard Island Reef", "liz"), ("MacGillivray Reef", "macg"),
                        ("North Direction Reef", "nd"), ("Eyrie Reef", "eyr")):
        r = result[result.reef_name == reef].set_index("treatment")
        summary[f"{short}_first"] = r.max_cpue_1995_2000
        summary[f"{short}_second"] = r.max_cpue_2010_2017
        summary[f"{short}_trough_cpue"] = r.mean_cpue_1999_2011
    liz = result[result.reef_name == "Lizard Island Reef"].set_index("treatment")
    # Within-complex exchange between Lizard sub-reefs is logged as internal immigration,
    # so the domain-sourced share (local + internal) is the comparable retention metric.
    summary["liz_domain_sourced_share_of_settlers"] = (liz.local_retention_1992_2015 + liz.internal_immigration_1992_2015) / (
        liz.local_retention_1992_2015 + liz.internal_immigration_1992_2015 + liz.external_immigration_1992_2015)

    fig, axes = plt.subplots(5, 1, figsize=(11, 15), sharex=True)
    for index, (reef, observed_name) in enumerate(REEFS.items()):
        axis = axes[index]
        for arm, level in zip(ARMS, LEVELS):
            x = trajectories[(trajectories.treatment == arm) & (trajectories.reef_name == reef) & (trajectories.year >= 1992)]
            axis.plot(x.year, x.simulated_cpue, color=COLORS[level], linewidth=2,
                      linestyle="--" if level in ("owen", "all50") else "-",
                      label=level if index == 0 else None)
        obs = cots_obs[(cots_obs.reef_name == observed_name) & cots_obs.year.between(1992, 2024)]
        axis.scatter(obs.year, obs.cotsptow, s=18, color="#0b0b0b", zorder=5, label="LTMP COTS/tow" if index == 0 else None)
        axis.axhline(.22, color="#52514e", linestyle=":", linewidth=.8)
        axis.set_ylabel(f"{reef}\nCOTS/tow (linear)")
        axis.set_ylim(bottom=0)
        axis.grid(color="#e4e3df", linewidth=.8)
    flux_ax = axes[4]
    for arm, level in zip(ARMS, LEVELS):
        f = flows[(flows.treatment == arm) & (flows.reef_name == "Lizard Island Reef") & (flows.year >= 1992)]
        flux_ax.plot(f.year, f.local_retention + f.internal_immigration, color=COLORS[level], linewidth=2,
                     linestyle="--" if level in ("owen", "all50") else "-")
    f = flows[(flows.treatment == REFERENCE) & (flows.reef_name == "Lizard Island Reef") & (flows.year >= 1992)]
    flux_ax.plot(f.year, f.external_immigration, color="#1baf7a", linewidth=1.5, linestyle=":", label="external (all arms)")
    flux_ax.set_ylabel("Lizard settlers /ha/yr\n(local + internal)")
    flux_ax.set_xlim(1991, 2024)
    flux_ax.grid(color="#e4e3df", linewidth=.8)
    flux_ax.legend(fontsize=8, frameon=False, loc="upper right")
    axes[0].set_title("Lizard-complex retention (rho of surviving larvae) vs Owen; all50 = every reef", loc="left", fontsize=10)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=6, fontsize=8, frameon=False, bbox_to_anchor=(.5, -.01))
    fig.tight_layout(rect=(0, .03, 1, 1))
    fig.savefig(paths["local_retention.png"], dpi=160, bbox_inches="tight", facecolor="#fcfcfb")
    plt.close(fig)

    report = {
        "status": "bounded_diagnostic_not_promoted",
        "reference_replay_max_rel_delta": replay,
        "external_immigration_max_abs_delta_vs_reference": external,
        "row_sum_check_in_metadata": "max_abs_row_sum_change" in metadata,
        "summary": json.loads(summary.reset_index().to_json(orient="records")),
        "inputs": {name: sha256(run / name) for name in (
            "metadata.toml", "scores.csv", "trajectories.csv", "flows.csv", "mechanisms.csv")},
        "script_sha256": sha256(Path(__file__)),
        "outputs": {name: sha256(path) for name, path in paths.items() if name != "local_retention_audit.json"},
    }
    paths["local_retention_audit.json"].write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    pd.set_option("display.width", 260)
    pd.set_option("display.max_columns", 30)
    print(summary.round(3).to_string())
    print(result[result.reef_name == "Lizard Island Reef"][[
        "retention", "local_retention_1992_2015", "internal_immigration_1992_2015",
        "external_immigration_1992_2015", "local_fecundity_1992_2015", "min_adults_ha_1999_2011"]].round(1).to_string(index=False))
    print("Reference replay max delta:", max(replay.values()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    main(parser.parse_args().run.resolve())
