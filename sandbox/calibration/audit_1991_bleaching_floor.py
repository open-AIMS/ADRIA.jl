"""Audit the preregistered bleaching-floor x starting fast-share screen
(crash_mechanism_protocol.md, 2026-10-14): ADRIA_BLEACHING_DHW_FLOOR {1, 4} x beta {1, 1.5}
x fast-share multiplier m {1, 2, 3} on the eb/lin baseline. COTS metrics, coral J
components (frozen) and the frozen m-selection rule; G2' comes from evaluate_gate_g2prime.py
on the run's re-score (path given)."""

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
from audit_1991_growth_starvation import PHOTO, PHOTO_REEFS  # noqa: E402
from audit_1991_initiation import exact_replay  # noqa: E402
from audit_1991_joint_coral import objective, reef_errors  # noqa: E402
from audit_1991_switching_repair import REEFS, sha256  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
PRIOR = ROOT / "sandbox/calibration/runs/20261013T_17_initial_state_rebaseline_1991_v1"
FLOORS, BETAS, MS = (100, 400), (100, 150), (100, 200, 300)
ARMS = [f"wc_b{b}_q20_eb_lin_f{f}_m{m}" for f, b, m in product(FLOORS, BETAS, MS)]
REFERENCE, PRIOR_REFERENCE = "wc_b150_q20_eb_lin_f100_m100", "wc_b150_q20_eb_lin"
N_YEARS = 34


def main(run: Path, gate: Path) -> None:
    paths = {name: run / name for name in ("bleaching_floor_audit.csv", "bleaching_floor_audit.json",
                                            "bleaching_floor.png")}
    if any(path.exists() for path in paths.values()):
        raise FileExistsError("Refusing to overwrite bleaching-floor audit")
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
    if not np.isfinite(mechanisms[["recruits_ha", "subadults_ha", "adults_ha", "fast_coral", "total_coral"]].to_numpy()).all():
        raise ValueError("Nonfinite state")
    start = mechanisms[mechanisms.year == 1991].pivot(index="reef_name", columns="treatment", values="total_coral")
    total_delta = float(start[ARMS].sub(start[REFERENCE], axis=0).abs().to_numpy().max())
    if total_delta > 1e-9:
        raise ValueError(f"Initial total cover differs between arms: {total_delta}")
    g2 = pd.read_csv(gate / "gate_g2prime.csv")
    g2 = g2[g2.run == run.name].set_index(["treatment", "a_half"])

    photo = observed_photo()
    manta_raw = pd.read_csv(ROOT / "sandbox/data/reef_manta.csv")
    manta_raw = manta_raw[(manta_raw.data_type == "manta") & (manta_raw.domain_category == "reef") &
                          (manta_raw.variable == "HC") & (manta_raw.purpose == "MANTA") &
                          (manta_raw.project_code == "LTMP") & (manta_raw.depth == 9)]
    cots_obs = pd.read_csv(ROOT / "sandbox/data/reef_cots.csv")
    rows = []
    for s in tables["scores"][tables["scores"].observation_treatment == "smoothed_3y"].itertuples(index=False):
        traj = trajectories[(trajectories.treatment == s.treatment) & (trajectories.reef_name == s.reef_name)].set_index("year")
        mech = mechanisms[(mechanisms.treatment == s.treatment) & (mechanisms.reef_name == s.reef_name)].set_index("year")
        mech = mech.assign(slow_coral=mech.total_coral - mech.fast_coral)
        manta = manta_raw[(manta_raw.reef_name == REEFS[s.reef_name]) &
                          manta_raw.report_year.between(1992, 2024)].groupby("report_year")["median"].mean()
        _, b, _, _, _, f, m = s.treatment.split("_")
        row = {"treatment": s.treatment, "reef_name": s.reef_name, "floor": int(f[1:]) / 100,
               "beta": int(b[1:]) / 100, "m": int(m[1:]) / 100, "smoothed_loss": float(s.loss),
               "matched_smoothed": int(s.matched_peaks),
               "first_wave_1995_2000": float(traj.loc[1995:2000].simulated_cpue.max()),
               "second_wave_2010_2017": float(traj.loc[2010:2017].simulated_cpue.max()),
               "inter_wave_cpue_1999_2011": float(traj.loc[1999:2011].simulated_cpue.mean()),
               "fast_1991": float(mech.loc[1991, "fast_coral"]),
               "fast_1993": float(mech.loc[1993, "fast_coral"]),
               "fast_1994": float(mech.loc[1994, "fast_coral"]),
               "mean_fast_1992_1996": float(mech.loc[1992:1996].fast_coral.mean()),
               "maturation_1994_1998": float(mech.loc[1994:1998].maturation_ha.sum()),
               "peak_adults_1995_2000": float(mech.loc[1995:2000].adults_ha.max()),
               "min_food_survival_1995_1998": float(mech.loc[1995:1998].food_survival.min()),
               "mean_fast_2005_2013": float(mech.loc[2005:2013].fast_coral.mean())}
        row |= reef_errors(mech, s.reef_name, photo, manta, (1992, 2024))
        rows.append(row)
    result = pd.DataFrame(rows)
    result.to_csv(paths["bleaching_floor_audit.csv"], index=False)

    raw = tables["scores"][tables["scores"].observation_treatment == "raw"].groupby("treatment").matched_peaks.sum()
    summary = result.groupby("treatment").agg(
        floor=("floor", "first"), beta=("beta", "first"), m=("m", "first"),
        summed_loss=("smoothed_loss", "sum"), matched_smoothed=("matched_smoothed", "sum"),
        rmse_fast=("rmse_fast_vs_photo", "mean"), rmse_slow=("rmse_slow_vs_photo", "mean"),
        rmse_manta_total=("rmse_total_vs_manta", "mean")).reindex(ARMS)
    summary["matched_raw"] = raw.reindex(ARMS)
    summary["coral_J"] = objective(result).reindex(ARMS)
    for reef, short in (("Lizard Island Reef", "liz"), ("MacGillivray Reef", "macg"),
                        ("North Direction Reef", "nd"), ("Eyrie Reef", "eyr")):
        r = result[result.reef_name == reef].set_index("treatment")
        summary[f"{short}_first"] = r.first_wave_1995_2000
        summary[f"{short}_second"] = r.second_wave_2010_2017
    liz = result[result.reef_name == "Lizard Island Reef"].set_index("treatment")
    summary["liz_fast_1993_to_1994"] = liz.fast_1993.round(3).astype(str) + "->" + liz.fast_1994.round(3).astype(str)
    summary["liz_peak_adults"] = liz.peak_adults_1995_2000
    summary["liz_maturation_1994_1998"] = liz.maturation_1994_1998
    for a_half in (0.0, 10.0):
        summary[f"g2_pass_ahalf{a_half:g}"] = [bool(g2.loc[(arm, a_half), "passes_g2prime"]) for arm in ARMS]
        summary[f"g2_liz_trough_ahalf{a_half:g}"] = [float(g2.loc[(arm, a_half), "liz_trough_model_mean_cpue"]) for arm in ARMS]
        summary[f"g2_macg_trough_ahalf{a_half:g}"] = [float(g2.loc[(arm, a_half), "macg_trough_model_mean_cpue"]) for arm in ARMS]

    # Frozen m-selection: per (floor, beta), eligible if G2' Lizard/MacG trough (linear) no worse
    # than the reference and manta total RMSE <= 1.1 x reference; minimum summed loss.
    ref = summary.loc[REFERENCE]
    selection = {}
    for f, b in product(FLOORS, BETAS):
        arms = summary.loc[[f"wc_b{b}_q20_eb_lin_f{f}_m{m}" for m in MS]]
        eligible = arms[(arms.g2_liz_trough_ahalf0 <= ref.g2_liz_trough_ahalf0 + 1e-12) &
                        (arms.g2_macg_trough_ahalf0 <= ref.g2_macg_trough_ahalf0 + 1e-12) &
                        (arms.rmse_manta_total <= 1.1 * ref.rmse_manta_total)]
        selection[f"floor{f}_beta{b}"] = {
            "eligible": list(eligible.index),
            "selected": eligible.summed_loss.idxmin() if len(eligible) else None}

    fig, axes = plt.subplots(4, 2, figsize=(14, 13), sharex=True)
    shown = {REFERENCE: ("#52514e", "--"), "wc_b150_q20_eb_lin_f400_m100": ("#2a78d6", "-"),
             "wc_b150_q20_eb_lin_f400_m200": ("#1c5cab", "-"), "wc_b150_q20_eb_lin_f400_m300": ("#0d366b", "-")}
    for index, (reef, observed_name) in enumerate(REEFS.items()):
        cots_ax, fast_ax = axes[index]
        for arm, (color, style) in shown.items():
            x = trajectories[(trajectories.treatment == arm) & (trajectories.reef_name == reef) & (trajectories.year >= 1992)]
            mm = mechanisms[(mechanisms.treatment == arm) & (mechanisms.reef_name == reef)]
            label = arm.replace("wc_b150_q20_eb_lin_", "") if index == 0 else None
            cots_ax.plot(x.year, x.simulated_cpue, color=color, linestyle=style, linewidth=1.8, label=label)
            fast_ax.plot(mm.year, mm.fast_coral, color=color, linestyle=style, linewidth=1.8)
        obs = cots_obs[(cots_obs.reef_name == observed_name) & cots_obs.year.between(1992, 2024)]
        cots_ax.scatter(obs.year, obs.cotsptow, s=16, color="#0b0b0b", zorder=5, label="LTMP" if index == 0 else None)
        if reef in PHOTO_REEFS:
            p = photo[photo.model_reef == reef]
            fast_ax.scatter(p.year, p.obs_fast, s=16, color="#0b0b0b", zorder=5)
        cots_ax.set_ylabel(f"{reef}\nCOTS/tow (linear q)")
        fast_ax.set_ylabel("fast coral (habitable)")
        for axis in (cots_ax, fast_ax):
            axis.set_xlim(1991, 2024)
            axis.set_ylim(bottom=0)
            axis.grid(color="#e4e3df", linewidth=.8)
    axes[0, 0].set_title("beta 1.5: floor 1 (ref) vs floor 4 x starting fast share m", loc="left", fontsize=10)
    axes[0, 1].set_title("Fast coral vs photo transects", loc="left", fontsize=10)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=5, fontsize=8, frameon=False, bbox_to_anchor=(.5, -.01))
    fig.tight_layout(rect=(0, .04, 1, 1))
    fig.savefig(paths["bleaching_floor.png"], dpi=160, bbox_inches="tight", facecolor="#fcfcfb")
    plt.close(fig)

    report = {"status": "bounded_diagnostic_not_promoted", "reference_replay_max_rel_delta": replay,
              "initial_total_cover_max_abs_delta": total_delta, "m_selection": selection,
              "summary": json.loads(summary.reset_index().to_json(orient="records")),
              "inputs": {name: sha256(run / name) for name in ("metadata.toml", "scores.csv", "trajectories.csv",
                                                                 "flows.csv", "mechanisms.csv")}
                        | {"gate_g2prime.csv": sha256(gate / "gate_g2prime.csv"), "reef_photo_transect.csv": sha256(PHOTO)},
              "script_sha256": sha256(Path(__file__)),
              "outputs": {name: sha256(path) for name, path in paths.items() if name != "bleaching_floor_audit.json"}}
    paths["bleaching_floor_audit.json"].write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    pd.set_option("display.width", 280)
    pd.set_option("display.max_columns", 40)
    print(summary.drop(columns=["floor", "beta", "m"]).round(3).to_string())
    print("m selection:", json.dumps(selection, indent=1))
    print("Reference replay max delta:", max(replay.values()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("gate", type=Path)
    args = parser.parse_args()
    main(args.run.resolve(), args.gate.resolve())
