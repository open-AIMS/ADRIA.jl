"""Audit the candidate set x replicate seeds (crash_mechanism_protocol.md, "Preregistered
candidate set and replicate-seed run", 2026-10-14; close-out B3 + B4 replicates + C3).

Each arm is read at its own observation q (obs015 arms from the q = 0.015 re-score / G2'
evaluation, obs030 arms from the --q=0.030 ones). Primary score: frozen score v2.

    python audit_1991_candidates.py <run> <rescore_q015> <rescore_q030> <gate_q015> <gate_q030>
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
RUNS = ROOT / "sandbox/calibration/runs"
BASE = "wc_b150_q20_eb_lin"
CANDIDATES = {  # id -> (observation q, suffix)
    "W": (0.015, "cv2_pnd"), "C1": (0.015, "cv2_pick"), "C2": (0.030, "cv2_pnd"),
    "C3a": (0.015, "cv2_pndl50"), "C3b": (0.030, "cv2_pndl50"), "C5": (0.015, "ref"),
    # C6 (sixth round): W with chronic mu 0.09; present from run 20261014T_48 onwards.
    "C6": (0.015, "m090_l025_cv2_pnd")}
RUN37 = "20261014T_37_candidates_1991_v1"
REPLAYS = {  # s1 arm -> (prior run, prior arm)
    f"{BASE}_obs015_s1_cv2_pick": ("20261014T_25_initial_share_1991_v2", f"{BASE}_cv2_pick"),
    f"{BASE}_obs015_s1_cv2_pnd": ("20261014T_32_timing_1991_v1", f"{BASE}_obs015_cv2_pnd"),
    f"{BASE}_obs030_s1_cv2_pnd": ("20261014T_32_timing_1991_v1", f"{BASE}_obs030_cv2_pnd"),
    f"{BASE}_obs015_s1_ref": ("20261014T_25_initial_share_1991_v2", f"{BASE}_ref")}
TIMING_REEFS = ("Lizard Island Reef", "MacGillivray Reef", "North Direction Reef")
SHORT = {"Lizard Island Reef": "liz", "MacGillivray Reef": "macg", "North Direction Reef": "nd", "Eyrie Reef": "eyr"}
N_SEEDS, N_YEARS = 5, 34


def arm(cid: str, seed: int) -> str:
    q, suffix = CANDIDATES[cid]
    return f"{BASE}_obs{round(1000 * q):03d}_s{seed}_{suffix}"


def years(text) -> list[int]:
    return [int(v) for v in str(text).split(";") if v and v != "nan"]


def main(run: Path, rescores: dict[float, Path], gates: dict[float, Path]) -> None:
    paths = {n: run / n for n in ("candidates_audit.csv", "candidates_replicates.csv",
                                  "candidates_audit.json", "candidates.png")}
    if any(p.exists() for p in paths.values()):
        raise FileExistsError("Refusing to overwrite candidate audit")
    tables = {n: pd.read_csv(run / f"{n}.csv") for n in ("trajectories", "flows", "mechanisms", "failures")}
    if len(tables["failures"]):
        raise ValueError("Failed evaluations in run")
    trajectories, mechanisms = tables["trajectories"], tables["mechanisms"]
    present = [c for c in CANDIDATES if arm(c, 1) in set(trajectories.treatment)]
    if "C6" in CANDIDATES and "C6" not in present:
        CANDIDATES.pop("C6")
    expected = {arm(c, s) for c in CANDIDATES for s in range(1, N_SEEDS + 1)}
    if set(trajectories.treatment) != expected or len(trajectories) != len(expected) * 4 * N_YEARS:
        raise ValueError("Missing or incomplete treatments")

    replays = dict(REPLAYS)
    if run.name != RUN37:
        # Later candidate runs: every original arm (all seeds) must replay run 37, and C6 at
        # seed pair 1 must replay the coral-mix arm of run 45.
        replays = {name: (RUN37, name) for name in expected if "_m090_" not in name}
        if "C6" in CANDIDATES:
            replays[arm("C6", 1)] = ("20261014T_45_coral_mix_1991_v1", arm("C6", 1))
    replay = {}
    for current, (prior_run, prior_arm) in replays.items():
        for name in ("trajectories", "flows", "mechanisms"):
            prior = pd.read_csv(RUNS / prior_run / f"{name}.csv")
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

    scores = {q: pd.read_csv(r / "detectability_scores.csv") for q, r in rescores.items()}
    series = {q: pd.read_csv(r / "detectability_series.csv") for q, r in rescores.items()}
    g2 = {q: pd.read_csv(g / "gate_g2prime.csv").set_index(["treatment", "a_half"]) for q, g in gates.items()}
    rows, reps = [], []
    for cid, (q, _) in CANDIDATES.items():
        for seed in range(1, N_SEEDS + 1):
            name = arm(cid, seed)
            s = scores[q]
            s = s[(s.run == run.name) & (s.treatment == name) & (s.a_half == 0.0)]
            primary = s[s.score_v2_primary.astype(str).str.lower() == "true"].set_index("reef_name")
            if len(primary) != 4 or len(s) != 8:
                raise ValueError(f"Missing re-score rows for {name}")
            x = series[q]
            x = x[(x.run == run.name) & (x.treatment == name) & (x.a_half == 0.0)]
            gate0, gate10, gate30 = (g2[q].loc[(name, a)] for a in (0.0, 10.0, 30.0))
            row = {"candidate": cid, "seed": seed, "treatment": name, "q": q,
                   "summed_loss_v1": float(s[s.observation_treatment == "smoothed_3y"].loss.sum()),
                   "summed_loss_v2": float(primary.loss_v2.sum()),
                   "matched_v2": int(primary.matched_peaks.sum()),
                   "observed_v2": int(sum(len(years(v)) for v in primary.obs_peak_years)),
                   "coral_J": float(coral_j[name]),
                   "g2_c1_c3_ahalf0": bool(gate0.c1_two_waves and gate0.c2_peaks_matched_raw_and_smoothed and
                                           gate0.c3_distinct_waves),
                   "g2_pass_ahalf10": bool(gate10.passes_g2prime),
                   "liz_trough_ahalf10": float(gate10.liz_trough_model_mean_cpue),
                   "macg_trough_ahalf10": float(gate10.macg_trough_model_mean_cpue),
                   "g2_pass_ahalf30": bool(gate30.passes_g2prime),
                   "liz_trough_ahalf30": float(gate30.liz_trough_model_mean_cpue),
                   "macg_trough_ahalf30": float(gate30.macg_trough_model_mean_cpue),
                   "liz_peak_adults_1995_2000": float(trajectories[(trajectories.treatment == name) &
                       (trajectories.reef_name == "Lizard Island Reef") &
                       trajectories.year.between(1995, 2000)].adults_ha.max())}
            errors = []
            for reef, short in SHORT.items():
                sim, obs = years(primary.loc[reef, "sim_peak_years"]), years(primary.loc[reef, "obs_peak_years"])
                cpue = x[x.reef_name == reef].set_index("year").detected_cpue
                first, second = cpue.loc[1993:2003], cpue.loc[2008:2018]
                y1, y2 = int(first.idxmax()), int(second.idxmax())
                trough = float(cpue.loc[y1:y2].min())
                row |= {f"{short}_loss_v2": float(primary.loc[reef, "loss_v2"]),
                        f"{short}_first_peak": sim[0] if sim else np.nan,
                        f"{short}_first_wave_cpue": float(first.max())}
                if reef in TIMING_REEFS and sim and obs:
                    errors.append(abs(sim[0] - obs[0]))
                reps.append({"candidate": cid, "seed": seed, "reef_name": reef,
                             "scored_first_peak": sim[0] if sim else np.nan,
                             "scored_second_peak": sim[1] if len(sim) > 1 else np.nan,
                             "first_peak_year": y1, "first_peak_cpue": float(first.max()),
                             "second_peak_year": y2, "second_peak_cpue": float(second.max()),
                             "trough_cpue": trough, "second_prominence_cpue": float(second.max()) - trough,
                             "period_years": y2 - y1, "loss_v2": float(primary.loc[reef, "loss_v2"])})
            row["mean_abs_first_peak_error"] = float(np.mean(errors)) if len(errors) == 3 else np.nan
            rows.append(row)
    result = pd.DataFrame(rows)
    replicates = pd.DataFrame(reps)
    result.to_csv(paths["candidates_audit.csv"], index=False)
    replicates.to_csv(paths["candidates_replicates.csv"], index=False)

    metrics = ["summed_loss_v2", "summed_loss_v1", "coral_J", "mean_abs_first_peak_error",
               "liz_peak_adults_1995_2000", "liz_trough_ahalf10", "macg_trough_ahalf10",
               "liz_trough_ahalf30", "macg_trough_ahalf30"]
    summary = result.groupby("candidate")[metrics].agg(["mean", "std"]).reindex(list(CANDIDATES))
    summary[("g2_pass_ahalf10", "seeds")] = result.groupby("candidate").g2_pass_ahalf10.sum()
    summary[("g2_pass_ahalf30", "seeds")] = result.groupby("candidate").g2_pass_ahalf30.sum()
    summary[("g2_c1_c3_ahalf0", "seeds")] = result.groupby("candidate").g2_c1_c3_ahalf0.sum()
    # AGENTS gate 3 (replicated seeds): better than C1 on v2 loss and first-peak timing in
    # the same seed, without worse coral J (tolerance 0.005), in at least 4 of 5 seeds.
    c1 = result[result.candidate == "C1"].set_index("seed")
    gate3 = {}
    for cid in CANDIDATES:
        r = result[result.candidate == cid].set_index("seed")
        better = ((r.summed_loss_v2 < c1.summed_loss_v2) &
                  (r.mean_abs_first_peak_error <= c1.mean_abs_first_peak_error) &
                  (r.coral_J <= c1.coral_J + 0.005))
        gate3[cid] = {"seeds_better_than_C1": int(better.sum()), "passes": bool(better.sum() >= 4)}
    reef_stats = (replicates.groupby(["candidate", "reef_name"])
                  [["first_peak_year", "first_peak_cpue", "second_peak_year", "second_peak_cpue",
                    "second_prominence_cpue", "period_years", "trough_cpue", "loss_v2"]]
                  .agg(["mean", "std"]))

    fig, axes = plt.subplots(4, 1, figsize=(11, 13), sharex=True)
    colors = {"W": "#1baf7a", "C1": "#52514e", "C2": "#3987e5", "C3a": "#eb6834", "C3b": "#c4458f", "C5": "#9a9893",
              "C6": "#0d366b"}
    cots_obs = pd.read_csv(ROOT / "sandbox/data/reef_cots.csv")
    for i, (reef, observed_name) in enumerate(REEFS.items()):
        axis = axes[i]
        for cid, (q, _) in CANDIDATES.items():
            x = series[q]
            x = x[(x.run == run.name) & (x.a_half == 0.0) & (x.reef_name == reef) & (x.year >= 1992) &
                  x.treatment.isin([arm(cid, s) for s in range(1, N_SEEDS + 1)])]
            band = x.groupby("year").detected_cpue.agg(["min", "max", "mean"])
            axis.fill_between(band.index, band["min"], band["max"], color=colors[cid], alpha=.15, linewidth=0)
            axis.plot(band.index, band["mean"], color=colors[cid], linewidth=1.8 if cid == "W" else 1.3,
                      label=f"{cid} (q {q})" if i == 0 else None)
        obs = cots_obs[(cots_obs.reef_name == observed_name) & cots_obs.year.between(1992, 2024)]
        axis.scatter(obs.year, obs.cotsptow, s=16, color="#0b0b0b", zorder=5, label="LTMP" if i == 0 else None)
        axis.set_ylabel(f"{reef}\nCOTS/tow")
        axis.set_xlim(1991, 2024)
        axis.set_ylim(bottom=0)
        axis.grid(color="#e4e3df", linewidth=.8)
    axes[0].set_title("Candidate set: mean and range over 5 seed pairs (each at its own observation q)",
                      loc="left", fontsize=10)
    axes[0].legend(fontsize=8, frameon=False, ncol=4)
    fig.tight_layout()
    fig.savefig(paths["candidates.png"], dpi=160, bbox_inches="tight", facecolor="#fcfcfb")
    plt.close(fig)

    report = {"status": "bounded_diagnostic_not_promoted", "reference_replay_max_rel_delta": replay,
              "gate3_replicated_seeds": gate3,
              "summary": {f"{a}|{b}": v for (a, b), v in summary.to_dict().items()},
              "reef_replicate_stats": {f"{cid}|{reef}": {f"{m}|{s}": float(v) for (m, s), v in r.items()}
                                       for (cid, reef), r in reef_stats.iterrows()},
              "inputs": {n: sha256(run / n) for n in ("metadata.toml", "trajectories.csv", "flows.csv", "mechanisms.csv")}
                        | {f"rescore_q{q}/detectability_scores.csv": sha256(r / "detectability_scores.csv") for q, r in rescores.items()}
                        | {f"gate_q{q}/gate_g2prime.csv": sha256(g / "gate_g2prime.csv") for q, g in gates.items()}
                        | {"reef_photo_transect.csv": sha256(PHOTO)},
              "script_sha256": sha256(Path(__file__)),
              "outputs": {n: sha256(p) for n, p in paths.items() if n != "candidates_audit.json"}}
    paths["candidates_audit.json"].write_text(json.dumps(report, indent=2, sort_keys=True, default=str) + "\n")
    pd.set_option("display.width", 320)
    pd.set_option("display.max_columns", 60)
    print(summary.round(3).to_string())
    s1 = result[result.seed == 1].set_index("candidate")
    print(s1[["summed_loss_v1", "summed_loss_v2", "coral_J", "liz_first_peak", "macg_first_peak", "nd_first_peak",
              "liz_first_wave_cpue", "macg_first_wave_cpue", "nd_first_wave_cpue", "liz_peak_adults_1995_2000",
              "liz_loss_v2", "macg_loss_v2", "nd_loss_v2", "eyr_loss_v2", "g2_pass_ahalf10", "g2_pass_ahalf30",
              "liz_trough_ahalf10", "macg_trough_ahalf10"]].round(3).to_string())
    print(reef_stats.round(2).to_string())
    print("Gate 3:", gate3)
    print("Replay max:", max(replay.values()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("run", "rescore_q015", "rescore_q030", "gate_q015", "gate_q030"):
        parser.add_argument(name, type=Path)
    a = parser.parse_args()
    main(a.run.resolve(), {0.015: a.rescore_q015.resolve(), 0.030: a.rescore_q030.resolve()},
         {0.015: a.gate_q015.resolve(), 0.030: a.gate_q030.resolve()})
