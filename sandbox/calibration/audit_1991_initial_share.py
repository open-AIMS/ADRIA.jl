"""Audit the starting fast-share partial calibration (crash_mechanism_protocol.md, 2026-10-14):
cv2 coral configuration x 1991 fast-share level {min, q25, q50, q75, max} of each reef's
healthy-year photo-transect range, plus `ref` (old coral, replays run 17) and `cv2_acrpoc`.
Stage 2b picks each calibration reef's level by its own smoothed frozen loss and prints the
OWEN_INITIAL_SHARE_PICKS string; a later `cv2_pick` arm (if present) is the confirmation."""

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
from audit_1991_switching_repair import REEFS, sha256  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
PRIOR = ROOT / "sandbox/calibration/runs/20261013T_17_initial_state_rebaseline_1991_v1"
LEVELS = ("min", "q25", "q50", "q75", "max")
BASE = "wc_b150_q20_eb_lin"
REFERENCE, PRIOR_REFERENCE = f"{BASE}_ref", BASE
CALIBRATION_REEFS = ("Lizard Island Reef", "MacGillivray Reef", "North Direction Reef")
N_YEARS = 34


def main(run: Path, gate: Path) -> None:
    paths = {name: run / name for name in ("initial_share_audit.csv", "initial_share_audit.json", "initial_share.png")}
    if any(path.exists() for path in paths.values()):
        raise FileExistsError("Refusing to overwrite initial-share audit")
    tables = {name: pd.read_csv(run / f"{name}.csv")
              for name in ("trajectories", "flows", "mechanisms", "scores", "failures")}
    if len(tables["failures"]):
        raise ValueError("Failed evaluations in run")
    trajectories, flows, mechanisms = (tables[k] for k in ("trajectories", "flows", "mechanisms"))
    arms = sorted(set(trajectories.treatment), key=lambda a: (a != REFERENCE, a))
    expected = {REFERENCE, f"{BASE}_cv2_acrpoc"} | {f"{BASE}_cv2_{lv}" for lv in LEVELS}
    if not expected <= set(arms) or len(trajectories) != len(arms) * 4 * N_YEARS:
        raise ValueError("Missing or incomplete treatments")
    replay = {}
    for name in ("trajectories", "flows", "mechanisms"):
        prior = pd.read_csv(PRIOR / f"{name}.csv")
        columns = [c for c in prior if c not in ("treatment", "reef_name", "year")]
        replay[name] = exact_replay(tables[name][tables[name].treatment == REFERENCE],
                                    prior[prior.treatment == PRIOR_REFERENCE], columns)
    ext_ref = flows[flows.treatment == REFERENCE].sort_values(["reef_name", "year"]).external_immigration.to_numpy()
    external = max(float(np.max(np.abs(flows[flows.treatment == a].sort_values(["reef_name", "year"])
                                       .external_immigration.to_numpy() - ext_ref))) for a in arms)
    if external > 1e-9:
        raise ValueError("External immigration differs between arms")
    start = mechanisms[mechanisms.year == 1991].pivot(index="reef_name", columns="treatment", values="total_coral")
    total_delta = float(start[arms].sub(start[REFERENCE], axis=0).abs().to_numpy().max())
    if total_delta > 1e-9:
        raise ValueError(f"1991 total cover differs between arms: {total_delta}")

    photo = observed_photo()
    manta_raw = pd.read_csv(ROOT / "sandbox/data/reef_manta.csv")
    manta_raw = manta_raw[(manta_raw.data_type == "manta") & (manta_raw.domain_category == "reef") &
                          (manta_raw.variable == "HC") & (manta_raw.purpose == "MANTA") &
                          (manta_raw.project_code == "LTMP") & (manta_raw.depth == 9)]
    g2 = pd.read_csv(gate / "gate_g2prime.csv")
    g2 = g2[g2.run == run.name].set_index(["treatment", "a_half"])
    cots_obs = pd.read_csv(ROOT / "sandbox/data/reef_cots.csv")
    raw = tables["scores"][tables["scores"].observation_treatment == "raw"].set_index(["treatment", "reef_name"])
    rows = []
    for s in tables["scores"][tables["scores"].observation_treatment == "smoothed_3y"].itertuples(index=False):
        traj = trajectories[(trajectories.treatment == s.treatment) & (trajectories.reef_name == s.reef_name)].set_index("year")
        mech = mechanisms[(mechanisms.treatment == s.treatment) & (mechanisms.reef_name == s.reef_name)].set_index("year")
        mech = mech.assign(slow_coral=mech.total_coral - mech.fast_coral)
        manta = manta_raw[(manta_raw.reef_name == REEFS[s.reef_name]) &
                          manta_raw.report_year.between(1992, 2024)].groupby("report_year")["median"].mean()
        row = {"treatment": s.treatment, "level": s.treatment.split("_")[-1], "reef_name": s.reef_name,
               "smoothed_loss": float(s.loss), "matched_smoothed": int(s.matched_peaks),
               "matched_raw": int(raw.loc[(s.treatment, s.reef_name)].matched_peaks),
               "fast_share_1991": float(mech.loc[1991, "fast_coral"] / mech.loc[1991, "total_coral"]),
               "first_wave_cpue": float(traj.loc[1995:2000].simulated_cpue.max()),
               "second_wave_cpue": float(traj.loc[2010:2017].simulated_cpue.max()),
               "peak_adults_1995_2000": float(traj.loc[1995:2000].adults_ha.max()),
               "fast_1992_1996": float(mech.loc[1992:1996].fast_coral.mean())}
        row |= reef_errors(mech, s.reef_name, photo, manta, (1992, 2024))
        rows.append(row)
    result = pd.DataFrame(rows)
    result.to_csv(paths["initial_share_audit.csv"], index=False)

    summary = result.groupby("treatment").agg(level=("level", "first"), summed_loss=("smoothed_loss", "sum"),
                                              matched_raw=("matched_raw", "sum"),
                                              matched_smoothed=("matched_smoothed", "sum")).reindex(arms)
    summary["coral_J"] = objective(result).reindex(arms)
    for reef, short in (("Lizard Island Reef", "liz"), ("MacGillivray Reef", "macg"),
                        ("North Direction Reef", "nd"), ("Eyrie Reef", "eyr")):
        r = result[result.reef_name == reef].set_index("treatment")
        summary[f"{short}_loss"] = r.smoothed_loss
        summary[f"{short}_share"] = r.fast_share_1991
        summary[f"{short}_first"] = r.first_wave_cpue
        summary[f"{short}_second"] = r.second_wave_cpue
    summary["liz_peak_adults"] = result[result.reef_name == "Lizard Island Reef"].set_index("treatment").peak_adults_1995_2000
    for a_half in (0.0, 10.0):
        summary[f"g2_pass_ahalf{a_half:g}"] = [bool(g2.loc[(a, a_half), "passes_g2prime"]) for a in arms]
        for short, col in (("liz", "liz_trough_model_mean_cpue"), ("macg", "macg_trough_model_mean_cpue")):
            summary[f"{short}_trough_ahalf{a_half:g}"] = [float(g2.loc[(a, a_half), col]) for a in arms]

    # Stage 2b: per-reef pick by its own smoothed loss among the five levels (ties to the median).
    level_arms = [f"{BASE}_cv2_{lv}" for lv in LEVELS]
    picks = {}
    for reef in CALIBRATION_REEFS:
        losses = result[(result.reef_name == reef) & result.treatment.isin(level_arms)].set_index("level").smoothed_loss
        best = losses.min()
        tied = [lv for lv in LEVELS if abs(losses[lv] - best) <= 1e-9]
        picks[reef] = "q50" if "q50" in tied else tied[0]
    pick_string = ";".join(f"{reef}={level}" for reef, level in picks.items())
    confirmation = {}
    if f"{BASE}_cv2_pick" in arms:
        p, a = summary.loc[f"{BASE}_cv2_pick"], summary.loc[f"{BASE}_cv2_acrpoc"]
        confirmation = {"summed_loss": float(p.summed_loss), "acrpoc_loss": float(a.summed_loss),
                        "liz_trough_ahalf10": float(p.liz_trough_ahalf10), "macg_trough_ahalf10": float(p.macg_trough_ahalf10),
                        "g2_ahalf10": bool(p.g2_pass_ahalf10), "acrpoc_g2_ahalf10": bool(a.g2_pass_ahalf10),
                        "g2_not_worse_than_acrpoc": bool(p.g2_pass_ahalf10 or not a.g2_pass_ahalf10)}

    fig, axes = plt.subplots(4, 1, figsize=(11, 13), sharex=True)
    colors = {"ref": "#52514e", "acrpoc": "#eb6834", "min": "#cde2fb", "q25": "#86b6ef", "q50": "#3987e5",
              "q75": "#1c5cab", "max": "#0d366b", "pick": "#1baf7a"}
    for index, (reef, observed_name) in enumerate(REEFS.items()):
        axis = axes[index]
        for arm in arms:
            level = arm.split("_")[-1]
            x = trajectories[(trajectories.treatment == arm) & (trajectories.reef_name == reef) & (trajectories.year >= 1992)]
            axis.plot(x.year, x.simulated_cpue, color=colors[level], linewidth=2.2 if level == "pick" else 1.6,
                      linestyle="--" if level == "ref" else "-", label=level if index == 0 else None)
        obs = cots_obs[(cots_obs.reef_name == observed_name) & cots_obs.year.between(1992, 2024)]
        axis.scatter(obs.year, obs.cotsptow, s=16, color="#0b0b0b", zorder=5, label="LTMP" if index == 0 else None)
        axis.set_ylabel(f"{reef}\nCOTS/tow (q 0.015)")
        axis.set_xlim(1991, 2024)
        axis.set_ylim(bottom=0)
        axis.grid(color="#e4e3df", linewidth=.8)
    axes[0].set_title("cv2 coral configuration x 1991 fast-share level (ref: old coral, acrpoc)", loc="left", fontsize=10)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=5, fontsize=8, frameon=False, bbox_to_anchor=(.5, -.01))
    fig.tight_layout(rect=(0, .04, 1, 1))
    fig.savefig(paths["initial_share.png"], dpi=160, bbox_inches="tight", facecolor="#fcfcfb")
    plt.close(fig)

    report = {"status": "bounded_diagnostic_not_promoted", "reference_replay_max_rel_delta": replay,
              "external_immigration_max_abs_delta": external, "initial_total_cover_max_abs_delta": total_delta,
              "per_reef_picks": picks, "picks_env": pick_string, "confirmation": confirmation,
              "summary": json.loads(summary.reset_index().to_json(orient="records")),
              "inputs": {name: sha256(run / name) for name in ("metadata.toml", "scores.csv", "trajectories.csv",
                                                                 "flows.csv", "mechanisms.csv")}
                        | {"gate_g2prime.csv": sha256(gate / "gate_g2prime.csv"), "reef_photo_transect.csv": sha256(PHOTO)},
              "script_sha256": sha256(Path(__file__)),
              "outputs": {name: sha256(path) for name, path in paths.items() if name != "initial_share_audit.json"}}
    paths["initial_share_audit.json"].write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    pd.set_option("display.width", 300)
    pd.set_option("display.max_columns", 40)
    print(summary.round(3).to_string())
    print("Per-reef picks:", picks)
    print("OWEN_INITIAL_SHARE_PICKS=" + pick_string)
    if confirmation:
        print("Confirmation:", json.dumps(confirmation, indent=1))
    print("Reference replay max delta:", max(replay.values()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("gate", type=Path)
    args = parser.parse_args()
    main(args.run.resolve(), args.gate.resolve())
