"""Audit the preregistered Lizard hard-coral initialisation screen (crash_mechanism_protocol.md,
2026-10-14): Lizard 1991 cover x c_L {1, 1.5, 2} x observation q {0.015, 0.030} x chronic mu
{0.03, 0.08}. Each arm is read under its own q: q = 0.015 arms from the q-0.015 re-score / G2'
evaluation, q = 0.030 arms from the q-0.030 ones (rescore_detectability.jl --q=0.030)."""

from __future__ import annotations

import argparse
import json
import sys
from itertools import product
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from audit_1991_coral_calibration import observed_photo  # noqa: E402
from audit_1991_growth_starvation import PHOTO  # noqa: E402
from audit_1991_initiation import exact_replay  # noqa: E402
from audit_1991_joint_coral import objective, reef_errors  # noqa: E402
from audit_1991_switching_repair import REEFS, sha256  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
PRIOR = ROOT / "sandbox/calibration/runs/20261013T_17_initial_state_rebaseline_1991_v1"
QS, MUS, CS = (15, 30), (30, 80), (100, 150, 200)
ARMS = [f"wc_b150_q20_eb_lin_obs{q:03d}_mu{m:03d}_c{c}" for q, m, c in product(QS, MUS, CS)]
REFERENCE, PRIOR_REFERENCE = "wc_b150_q20_eb_lin_obs015_mu030_c100", "wc_b150_q20_eb_lin"
PRATCHETT_1996_COVER, PRATCHETT_PEAK_DENSITY = 0.322, (40.0, 55.0)
N_YEARS = 34


def main(run: Path, rescores: dict, gates: dict) -> None:
    paths = {name: run / name for name in ("lizard_coral_audit.csv", "lizard_coral_audit.json", "lizard_coral.png")}
    if any(path.exists() for path in paths.values()):
        raise FileExistsError("Refusing to overwrite Lizard-coral audit")
    tables = {name: pd.read_csv(run / f"{name}.csv") for name in ("trajectories", "flows", "mechanisms", "failures")}
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
    external = max(float(np.max(np.abs(flows[flows.treatment == a].sort_values(["reef_name", "year"]).external_immigration.to_numpy()
                                       - reference.external_immigration.to_numpy()))) for a in ARMS)
    if external > 1e-9:
        raise ValueError("External immigration differs between arms")
    if not np.isfinite(mechanisms[["recruits_ha", "subadults_ha", "adults_ha", "fast_coral", "total_coral"]].to_numpy()).all():
        raise ValueError("Nonfinite state")

    photo = observed_photo()
    manta_raw = pd.read_csv(ROOT / "sandbox/data/reef_manta.csv")
    manta_raw = manta_raw[(manta_raw.data_type == "manta") & (manta_raw.domain_category == "reef") &
                          (manta_raw.variable == "HC") & (manta_raw.purpose == "MANTA") &
                          (manta_raw.project_code == "LTMP") & (manta_raw.depth == 9)]
    scores = {q: pd.read_csv(rescores[q] / "detectability_scores.csv") for q in QS}
    series = {q: pd.read_csv(rescores[q] / "detectability_series.csv") for q in QS}
    g2 = {q: pd.read_csv(gates[q] / "gate_g2prime.csv").set_index(["treatment", "a_half"]) for q in QS}
    cots_obs = pd.read_csv(ROOT / "sandbox/data/reef_cots.csv")
    rows, coral_rows = [], []
    for arm in ARMS:
        q = int(arm.split("_obs")[1][:3])
        s = scores[q][(scores[q].treatment == arm) & (scores[q].a_half == 0.0)]
        for reef, observed_name in REEFS.items():
            sm = s[(s.reef_name == reef) & (s.observation_treatment == "smoothed_3y")].iloc[0]
            raw = s[(s.reef_name == reef) & (s.observation_treatment == "raw")].iloc[0]
            x = series[q][(series[q].treatment == arm) & (series[q].reef_name == reef) & (series[q].a_half == 0.0)].set_index("year")
            mech = mechanisms[(mechanisms.treatment == arm) & (mechanisms.reef_name == reef)].set_index("year")
            mech = mech.assign(slow_coral=mech.total_coral - mech.fast_coral)
            manta = manta_raw[(manta_raw.reef_name == observed_name) &
                              manta_raw.report_year.between(1992, 2024)].groupby("report_year")["median"].mean()
            row = {"treatment": arm, "reef_name": reef, "q": q / 1000, "mu": int(arm.split("_mu")[1][:3]) / 1000,
                   "c_lizard": int(arm.split("_c")[-1]) / 100, "loss_own_q": float(sm.loss),
                   "matched_smoothed": int(sm.matched_peaks), "matched_raw": int(raw.matched_peaks),
                   "first_wave_cpue": float(x.loc[1995:2000].detected_cpue.max()),
                   "second_wave_cpue": float(x.loc[2010:2017].detected_cpue.max()),
                   "peak_adults_1995_2000": float(mech.loc[1995:2000].adults_ha.max()),
                   "total_habitable_1996": float(mech.loc[1996, "total_coral"]),
                   "total_site_1996": float(mech.loc[1996, "total_coral_site_fraction"]),
                   "fast_1992_1996": float(mech.loc[1992:1996].fast_coral.mean())}
            row |= reef_errors(mech, reef, photo, manta, (1992, 2024))
            rows.append(row)
    result = pd.DataFrame(rows)
    result.to_csv(paths["lizard_coral_audit.csv"], index=False)
    summary = result.groupby("treatment").agg(
        q=("q", "first"), mu=("mu", "first"), c_lizard=("c_lizard", "first"),
        summed_loss_own_q=("loss_own_q", "sum"), matched_raw=("matched_raw", "sum"),
        matched_smoothed=("matched_smoothed", "sum"), rmse_manta=("rmse_total_vs_manta", "mean")).reindex(ARMS)
    summary["coral_J"] = objective(result).reindex(ARMS)
    for reef, short in (("Lizard Island Reef", "liz"), ("MacGillivray Reef", "macg"),
                        ("North Direction Reef", "nd"), ("Eyrie Reef", "eyr")):
        r = result[result.reef_name == reef].set_index("treatment")
        summary[f"{short}_first"] = r.first_wave_cpue
        summary[f"{short}_second"] = r.second_wave_cpue
    liz = result[result.reef_name == "Lizard Island Reef"].set_index("treatment")
    summary["liz_peak_adults"] = liz.peak_adults_1995_2000
    summary["liz_total_habitable_1996"] = liz.total_habitable_1996
    summary["liz_rmse_manta"] = liz.rmse_total_vs_manta
    for a_half in (0.0, 10.0):
        summary[f"g2_pass_ahalf{a_half:g}"] = [bool(g2[int(a.split("_obs")[1][:3])].loc[(a, a_half), "passes_g2prime"]) for a in ARMS]
        summary[f"liz_trough_ahalf{a_half:g}"] = [float(g2[int(a.split("_obs")[1][:3])].loc[(a, a_half), "liz_trough_model_mean_cpue"]) for a in ARMS]
        summary[f"macg_trough_ahalf{a_half:g}"] = [float(g2[int(a.split("_obs")[1][:3])].loc[(a, a_half), "macg_trough_model_mean_cpue"]) for a in ARMS]

    fig, axes = plt.subplots(2, 2, figsize=(14, 8), sharex=True)
    colors = {100: "#86b6ef", 150: "#2a78d6", 200: "#0d366b"}
    for col, q in enumerate(QS):
        cots_ax, coral_ax = axes[0, col], axes[1, col]
        for c in CS:
            arm = f"wc_b150_q20_eb_lin_obs{q:03d}_mu030_c{c}"
            x = series[q][(series[q].treatment == arm) & (series[q].reef_name == "Lizard Island Reef") &
                          (series[q].a_half == 0.0) & (series[q].year >= 1992)]
            mech = mechanisms[(mechanisms.treatment == arm) & (mechanisms.reef_name == "Lizard Island Reef")]
            cots_ax.plot(x.year, x.detected_cpue, color=colors[c], linewidth=2, label=f"c_L {c / 100:g}")
            coral_ax.plot(mech.year, mech.total_coral, color=colors[c], linewidth=2)
        obs = cots_obs[(cots_obs.reef_name == "Lizard Isles") & cots_obs.year.between(1992, 2024)]
        cots_ax.scatter(obs.year, obs.cotsptow, s=16, color="#0b0b0b", zorder=5, label="LTMP COTS/tow")
        m = manta_raw[(manta_raw.reef_name == "Lizard Isles") & manta_raw.report_year.between(1992, 2024)]
        coral_ax.scatter(m.report_year, m["median"], s=16, color="#0b0b0b", zorder=5, label="LTMP manta HC")
        coral_ax.scatter([1996.8, 1999.0], [0.322, 0.229], s=40, marker="D", color="#eb6834", zorder=6,
                         label="Pratchett 2005 reef-wide HC")
        cots_ax.set_title(f"Lizard, observation q = {q / 1000:g} (mu 0.03)", loc="left", fontsize=10)
        cots_ax.set_ylabel("COTS/tow")
        coral_ax.set_ylabel("total hard coral (habitable)")
        for axis in (cots_ax, coral_ax):
            axis.grid(color="#e4e3df", linewidth=.8)
            axis.set_xlim(1991, 2024)
            axis.set_ylim(bottom=0)
            axis.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(paths["lizard_coral.png"], dpi=160, facecolor="#fcfcfb")
    plt.close(fig)

    report = {"status": "bounded_diagnostic_not_promoted", "reference_replay_max_rel_delta": replay,
              "external_immigration_max_abs_delta": external,
              "pratchett_reference": {"lizard_hc_1996": PRATCHETT_1996_COVER, "peak_density_ha": PRATCHETT_PEAK_DENSITY},
              "summary": json.loads(summary.reset_index().to_json(orient="records")),
              "inputs": {name: sha256(run / name) for name in ("metadata.toml", "trajectories.csv", "flows.csv", "mechanisms.csv")}
                        | {f"rescore_q{q}": sha256(rescores[q] / "detectability_scores.csv") for q in QS}
                        | {f"g2_q{q}": sha256(gates[q] / "gate_g2prime.csv") for q in QS},
              "script_sha256": sha256(Path(__file__)),
              "outputs": {name: sha256(path) for name, path in paths.items() if name != "lizard_coral_audit.json"}}
    paths["lizard_coral_audit.json"].write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    pd.set_option("display.width", 300)
    pd.set_option("display.max_columns", 40)
    print(summary.drop(columns=["q", "mu", "c_lizard"]).round(3).to_string())
    print("Reference replay max delta:", max(replay.values()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("--rescore015", type=Path, required=True)
    parser.add_argument("--rescore030", type=Path, required=True)
    parser.add_argument("--gate015", type=Path, required=True)
    parser.add_argument("--gate030", type=Path, required=True)
    a = parser.parse_args()
    main(a.run.resolve(), {15: a.rescore015.resolve(), 30: a.rescore030.resolve()},
         {15: a.gate015.resolve(), 30: a.gate030.resolve()})
