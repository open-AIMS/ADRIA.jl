"""Audit the intermediate coral settings (crash_mechanism_protocol.md, "Preregistered
intermediate coral settings", 2026-10-14): working configuration x chronic mu {0.03, 0.06,
0.09, 0.12} x bleaching multiplier below 4 DHW {1, 0.5, 0.25}, plus the run-17 coral control.

    python audit_1991_coral_mix.py <run> <rescore_q015> <gate_q015>
"""

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
from audit_1991_growth_starvation import PHOTO  # noqa: E402
from audit_1991_initiation import exact_replay  # noqa: E402
from audit_1991_joint_coral import objective, reef_errors  # noqa: E402
from audit_1991_switching_repair import REEFS, sha256  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
PRIOR = ROOT / "sandbox/calibration/runs/20261014T_37_candidates_1991_v1"
BASE = "wc_b150_q20_eb_lin_obs015_s1"
REPLAYS = {f"{BASE}_m120_l025_cv2_pnd": f"{BASE}_cv2_pnd", f"{BASE}_ref": f"{BASE}_ref"}
TROUGH_BOUND = None  # read from the gate output
J_MAX, LOSS_SLACK = 0.30, 0.5
TIMING_REEFS = ("Lizard Island Reef", "MacGillivray Reef", "North Direction Reef")
SHORT = {"Lizard Island Reef": "liz", "MacGillivray Reef": "macg", "North Direction Reef": "nd", "Eyrie Reef": "eyr"}


def years(text) -> list[int]:
    return [int(v) for v in str(text).split(";") if v and v != "nan"]


def main(run: Path, rescore: Path, gate: Path) -> None:
    paths = {n: run / n for n in ("coral_mix_audit.csv", "coral_mix_audit.json", "coral_mix.png")}
    if any(p.exists() for p in paths.values()):
        raise FileExistsError("Refusing to overwrite coral-mix audit")
    tables = {n: pd.read_csv(run / f"{n}.csv") for n in ("trajectories", "flows", "mechanisms", "failures")}
    if len(tables["failures"]):
        raise ValueError("Failed evaluations in run")
    trajectories, mechanisms = tables["trajectories"], tables["mechanisms"]
    arms = sorted(set(trajectories.treatment))
    if len(arms) != 13:
        raise ValueError("Expected 13 arms")
    replay = {}
    for current, prior_arm in REPLAYS.items():
        for name in ("trajectories", "flows", "mechanisms"):
            prior = pd.read_csv(PRIOR / f"{name}.csv")
            columns = [c for c in prior if c not in ("treatment", "reef_name", "year")]
            replay[f"{current}/{name}"] = exact_replay(tables[name][tables[name].treatment == current],
                                                       prior[prior.treatment == prior_arm], columns)

    photo = observed_photo()
    manta_raw = pd.read_csv(ROOT / "sandbox/data/reef_manta.csv")
    manta_raw = manta_raw[(manta_raw.data_type == "manta") & (manta_raw.domain_category == "reef") &
                          (manta_raw.variable == "HC") & (manta_raw.purpose == "MANTA") &
                          (manta_raw.project_code == "LTMP") & (manta_raw.depth == 9)]
    coral_rows = []
    for (name, reef), mech in mechanisms.groupby(["treatment", "reef_name"]):
        mech = mech.set_index("year").assign(slow_coral=lambda m: m.total_coral - m.fast_coral)
        manta = manta_raw[(manta_raw.reef_name == REEFS[reef]) &
                          manta_raw.report_year.between(1992, 2024)].groupby("report_year")["median"].mean()
        coral_rows.append({"treatment": name, "reef_name": reef} | reef_errors(mech, reef, photo, manta, (1992, 2024)))
    coral_j = objective(pd.DataFrame(coral_rows))

    scores = pd.read_csv(rescore / "detectability_scores.csv")
    series = pd.read_csv(rescore / "detectability_series.csv")
    gate_json = json.loads((gate / "gate_g2prime.json").read_text())
    bounds = {r: v["poisson_upper_95"] for r, v in gate_json["observed_troughs"].items()}
    g2 = pd.read_csv(gate / "gate_g2prime.csv").set_index(["treatment", "a_half"])
    rows = []
    for name in arms:
        s = scores[(scores.run == run.name) & (scores.treatment == name) & (scores.a_half == 0.0)]
        primary = s[s.score_v2_primary.astype(str).str.lower() == "true"].set_index("reef_name")
        parts = name.split("_")
        mu = int(parts[7][1:]) / 1000 if name.endswith("_cv2_pnd") else 0.03
        low = int(parts[8][1:]) / 100 if name.endswith("_cv2_pnd") else 1.0
        row = {"treatment": name, "label": "C5 (run-17 coral, acrpoc shares)" if name.endswith("_ref") else f"mu {mu:.2f}, low-DHW x{low:g}",
               "mu": mu, "low_dhw_scale": low, "control": name.endswith("_ref"),
               "summed_loss_v2": float(primary.loss_v2.sum()), "matched_v2": int(primary.matched_peaks.sum()),
               "coral_J": float(coral_j[name]),
               "liz_peak_adults_1995_2000": float(trajectories[(trajectories.treatment == name) &
                   (trajectories.reef_name == "Lizard Island Reef") & trajectories.year.between(1995, 2000)].adults_ha.max())}
        for a_half in (0.0, 10.0, 30.0):
            gate_row = g2.loc[(name, a_half)]
            tag = f"ahalf{a_half:g}"
            row |= {f"liz_trough_{tag}": float(gate_row.liz_trough_model_mean_cpue),
                    f"macg_trough_{tag}": float(gate_row.macg_trough_model_mean_cpue),
                    f"g2_pass_{tag}": bool(gate_row.passes_g2prime),
                    f"g2_c1_c3_{tag}": bool(gate_row.c1_two_waves and gate_row.c2_peaks_matched_raw_and_smoothed and
                                            gate_row.c3_distinct_waves)}
            row[f"worst_trough_ratio_{tag}"] = max(row[f"liz_trough_{tag}"] / bounds["Lizard Island Reef"],
                                                   row[f"macg_trough_{tag}"] / bounds["MacGillivray Reef"])
        errors = []
        for reef, short in SHORT.items():
            sim, obs = years(primary.loc[reef, "sim_peak_years"]), years(primary.loc[reef, "obs_peak_years"])
            cpue = series[(series.run == run.name) & (series.treatment == name) & (series.a_half == 0.0) &
                          (series.reef_name == reef)].set_index("year").detected_cpue
            row |= {f"{short}_loss_v2": float(primary.loc[reef, "loss_v2"]), f"{short}_first_peak": sim[0] if sim else np.nan,
                    f"{short}_first_wave_cpue": float(cpue.loc[1993:2003].max()),
                    f"{short}_second_wave_cpue": float(cpue.loc[2008:2018].max())}
            if reef in TIMING_REEFS and sim and obs:
                errors.append(abs(sim[0] - obs[0]))
        row["mean_abs_first_peak_error"] = float(np.mean(errors)) if len(errors) == 3 else np.nan
        rows.append(row)
    result = pd.DataFrame(rows).sort_values(["control", "mu", "low_dhw_scale"]).reset_index(drop=True)

    w = result[result.treatment == f"{BASE}_m120_l025_cv2_pnd"].iloc[0]
    eligible = result[(~result.control) & (result.coral_J <= J_MAX) &
                      (result.summed_loss_v2 <= w.summed_loss_v2 + LOSS_SLACK)]
    strict = eligible[eligible.worst_trough_ratio_ahalf10 <= 1.0]
    weak = eligible[eligible.worst_trough_ratio_ahalf30 <= 1.0]
    # Trade-off front: arms not dominated on (coral J, worst trough ratio at A_half 10).
    front = [r.treatment for r in result.itertuples() if not any(
        (o.coral_J <= r.coral_J and o.worst_trough_ratio_ahalf10 <= r.worst_trough_ratio_ahalf10 and
         (o.coral_J < r.coral_J or o.worst_trough_ratio_ahalf10 < r.worst_trough_ratio_ahalf10))
        for o in result.itertuples())]
    result["on_front"] = result.treatment.isin(front)
    result.to_csv(paths["coral_mix_audit.csv"], index=False)

    fig, ax = plt.subplots(figsize=(8, 5.5))
    colors = {1.0: "#0d366b", 0.5: "#3987e5", 0.25: "#86b6ef"}
    for r in result.itertuples():
        color = "#eb6834" if r.control else colors[r.low_dhw_scale]
        ax.scatter(r.coral_J, r.worst_trough_ratio_ahalf10, s=30 + 400 * (r.mu if not r.control else 0.03),
                   color=color, edgecolors="#fcfcfb", linewidths=1.5, zorder=3)
        ax.annotate(r.label.replace("low-DHW ", ""), (r.coral_J, r.worst_trough_ratio_ahalf10), fontsize=7,
                    xytext=(5, 3), textcoords="offset points", color="#52514e")
    ax.axhline(1.0, color="#9a9893", linewidth=1, linestyle="--")
    ax.axvline(J_MAX, color="#9a9893", linewidth=1, linestyle=":")
    ax.set_xlabel("Coral fit J (lower is better)")
    ax.set_ylabel("Worst between-wave COTS / observed bound (A½ 10)")
    ax.set_title("Coral fit vs between-wave COTS (marker size: chronic background mu)", loc="left", fontsize=10)
    ax.grid(color="#e4e3df", linewidth=.8)
    fig.tight_layout()
    fig.savefig(paths["coral_mix.png"], dpi=160, facecolor="#fcfcfb")
    plt.close(fig)

    report = {"status": "bounded_diagnostic_not_promoted", "reference_replay_max_rel_delta": replay,
              "trough_bounds": bounds, "selection": {
                  "rule": f"coral J <= {J_MAX}, loss v2 <= W + {LOSS_SLACK}; troughs pass at A_half 10 (strict) or 30 (weak)",
                  "eligible": eligible.treatment.tolist(), "strict": strict.treatment.tolist(),
                  "weak": weak.treatment.tolist(), "front_A_half_10": front},
              "summary": json.loads(result.to_json(orient="records")),
              "inputs": {n: sha256(run / n) for n in ("metadata.toml", "trajectories.csv", "flows.csv", "mechanisms.csv")}
                        | {"detectability_scores.csv": sha256(rescore / "detectability_scores.csv"),
                           "gate_g2prime.csv": sha256(gate / "gate_g2prime.csv"), "reef_photo_transect.csv": sha256(PHOTO)},
              "script_sha256": sha256(Path(__file__)),
              "outputs": {n: sha256(p) for n, p in paths.items() if n != "coral_mix_audit.json"}}
    paths["coral_mix_audit.json"].write_text(json.dumps(report, indent=2, sort_keys=True, default=str) + "\n")
    pd.set_option("display.width", 320)
    pd.set_option("display.max_columns", 60)
    print(result[["label", "summed_loss_v2", "coral_J", "worst_trough_ratio_ahalf10", "worst_trough_ratio_ahalf30",
                  "liz_trough_ahalf10", "macg_trough_ahalf10", "g2_pass_ahalf10", "g2_pass_ahalf30",
                  "liz_first_peak", "macg_first_peak", "nd_first_peak", "mean_abs_first_peak_error",
                  "liz_peak_adults_1995_2000", "nd_first_wave_cpue", "macg_second_wave_cpue", "on_front"]]
          .round(3).to_string(index=False))
    print("Bounds:", bounds)
    print("Selection:", json.dumps(report["selection"], indent=1))
    print("Replay max:", max(replay.values()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("rescore", type=Path)
    parser.add_argument("gate", type=Path)
    a = parser.parse_args()
    main(a.run.resolve(), a.rescore.resolve(), a.gate.resolve())
