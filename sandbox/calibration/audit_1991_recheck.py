"""Audit the candidate re-check (crash_mechanism_protocol.md, "Preregistered candidate re-check",
2026-10-15, fourteenth round): W, C3a and C6 with the new coral and the 1990s-corrected chl
table, x five seed pairs x both outbreak layers, on domain v0.4 to 2026.

    python audit_1991_recheck.py <run_original> <rescore_original> <gate_original> \\
        <run_corrected> <rescore_corrected> <gate_corrected>
Outputs go to <run_corrected>.
"""

from __future__ import annotations

import argparse
import json
import re
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
RUNS = ROOT / "sandbox/calibration/runs"
TARGETS = ROOT / "sandbox/data/derived/density_targets_2020s/targets_reef_year.csv"
PRIOR = {"original": "20261015T_22_ext_chl2_original_1991_v1", "corrected": "20261015T_23_ext_chl2_corrected_1991_v1"}
REPLAY_ARM = "wc_b100_q20_eb_lin_obs015_s1_t60_e50_aprma_cv2_pnd"
WAVES = {"w3": (1993, 2001), "w4": (2008, 2015), "w5": (2019, 2026)}
SHORT = {"Lizard Island Reef": "liz", "MacGillivray Reef": "macg", "North Direction Reef": "nd", "Eyrie Reef": "eyr"}
ONSET = 3.0
GATE_MIN_PASSES, V3_MARGIN = 4, 1.5


def candidate(arm: str) -> str:
    if arm == REPLAY_ARM:
        return "W_uncorrected_chl"
    if "_m090_" in arm:
        return "C6"
    return "C3a" if arm.endswith("_cv2_pndl50") else "W"


def seed(arm: str) -> int:
    return int(re.search(r"_s(\d+)_", arm).group(1))


def wave(adults: pd.Series, window: tuple[int, int]) -> dict:
    a = adults.loc[window[0]:window[1]]
    above = a[a >= ONSET]
    return {"onset": int(above.index.min()) if len(above) else np.nan, "peak": float(a.max()),
            "peak_year": int(a.idxmax())}


def main(runs: dict[str, tuple[Path, Path, Path]]) -> None:
    out = runs["corrected"][0]
    paths = {n: out / n for n in ("recheck_audit.csv", "recheck_summary.csv", "recheck_audit.json", "recheck.png")}
    if any(p.exists() for p in paths.values()):
        raise FileExistsError("Refusing to overwrite recheck audit")
    targets = pd.read_csv(TARGETS)
    rows, replay, series = [], {}, {}
    for ds, (run, rescore, gate) in runs.items():
        tables = {n: pd.read_csv(run / f"{n}.csv") for n in ("trajectories", "flows", "mechanisms", "failures")}
        if len(tables["failures"]):
            raise ValueError(f"Failed evaluations in {run.name}")
        for name in ("trajectories", "flows", "mechanisms"):
            prior = pd.read_csv(RUNS / PRIOR[ds] / f"{name}.csv")
            columns = [c for c in prior if c not in ("treatment", "reef_name", "year")]
            replay[f"{ds}/{name}"] = exact_replay(tables[name][tables[name].treatment == REPLAY_ARM],
                                                  prior[prior.treatment == REPLAY_ARM], columns)
        scores = pd.read_csv(rescore / "detectability_scores.csv")
        g2 = pd.read_csv(gate / "gate_g2prime.csv").set_index(["treatment", "a_half"])
        traj = tables["trajectories"]
        for arm in sorted(set(traj.treatment)):
            s = scores[(scores.treatment == arm) & (scores.a_half == 0.0)]
            primary = s[s.score_v2_primary.astype(str).str.lower() == "true"].set_index("reef_name")
            row = {"dataset": ds, "candidate": candidate(arm), "seed": seed(arm), "treatment": arm,
                   "loss_v3": float(primary.loss_v3.sum()), "loss_v2_manta": float(primary.loss_v2_manta.sum()),
                   "diver_error": float(primary.diver_error.sum()), "matched_v2": int(primary.matched_peaks.sum()),
                   "g2_pass_ahalf30": bool(g2.loc[(arm, 30.0)].passes_g2prime),
                   "g2_pass_ahalf10": bool(g2.loc[(arm, 10.0)].passes_g2prime)}
            for reef, short in SHORT.items():
                a = traj[(traj.treatment == arm) & (traj.reef_name == reef)].set_index("year").adults_ha
                series[(ds, arm, reef)] = a
                row[f"{short}_diver_error"] = float(primary.loc[reef, "diver_error"])
                for w, window in WAVES.items():
                    for k, v in wave(a, window).items():
                        row[f"{short}_{w}_{k}"] = v
                tg = targets[targets.target_reef == reef].set_index("report_year").target_density_ha
                for y in (2025, 2026):
                    row[f"{short}_{y}_model"] = float(a.get(y, np.nan))
                    row[f"{short}_{y}_target"] = float(tg.get(y, np.nan))
            rows.append(row)
    result = pd.DataFrame(rows).sort_values(["dataset", "candidate", "seed"]).reset_index(drop=True)
    result.to_csv(paths["recheck_audit.csv"], index=False)

    cands = result[result.candidate != "W_uncorrected_chl"]
    summary = cands.groupby(["dataset", "candidate"]).agg(
        v3_mean=("loss_v3", "mean"), v3_min=("loss_v3", "min"), v3_max=("loss_v3", "max"),
        manta_mean=("loss_v2_manta", "mean"), diver_mean=("diver_error", "mean"),
        g2_passes_ahalf30=("g2_pass_ahalf30", "sum"), g2_passes_ahalf10=("g2_pass_ahalf10", "sum"),
        n_seeds=("seed", "nunique")).reset_index()
    summary["v3_spread"] = summary.v3_max - summary.v3_min
    for short in SHORT.values():
        for col in ("w3_onset", "w4_onset", "w5_onset", "w5_peak"):
            summary[f"{short}_{col}_mean"] = cands.groupby(["dataset", "candidate"])[f"{short}_{col}"].mean().values
        for y in (2025, 2026):
            ratio = np.log((cands[f"{short}_{y}_model"] + 0.5) / (cands[f"{short}_{y}_target"] + 0.5))
            summary[f"{short}_{y}_model_over_target_gmean"] = np.exp(ratio.groupby([cands.dataset, cands.candidate]).mean()).values
    gate_ok = summary.groupby("candidate").g2_passes_ahalf30.min() >= GATE_MIN_PASSES
    best = summary.groupby("candidate").v3_mean.mean()
    retained = sorted(c for c in gate_ok.index if gate_ok[c] and best[c] <= best[gate_ok[gate_ok].index].min() + V3_MARGIN) \
        if gate_ok.any() else []
    summary.to_csv(paths["recheck_summary.csv"], index=False)
    w1 = result[(result.candidate == "W") & (result.seed == 1)].set_index("dataset")
    unc = result[result.candidate == "W_uncorrected_chl"].set_index("dataset")
    correction_effect = {ds: {"d_loss_v3": float(w1.loc[ds, "loss_v3"] - unc.loc[ds, "loss_v3"]),
                              "d_loss_v2_manta": float(w1.loc[ds, "loss_v2_manta"] - unc.loc[ds, "loss_v2_manta"]),
                              "d_diver_error": float(w1.loc[ds, "diver_error"] - unc.loc[ds, "diver_error"])}
                         for ds in runs}

    colors = {"W": "#0d366b", "C3a": "#eb6834", "C6": "#3a9e6f"}
    styles = {"W": "-", "C3a": "--", "C6": "-."}
    fig, axes = plt.subplots(4, 2, figsize=(13, 13), sharex=True)
    for col, ds in enumerate(("original", "corrected")):
        for r, reef in enumerate(REEFS):
            ax = axes[r, col]
            for cand in ("W", "C3a", "C6"):
                arms = cands[(cands.dataset == ds) & (cands.candidate == cand)].treatment
                m = pd.concat([series[(ds, a, reef)] for a in arms], axis=1)
                ax.fill_between(m.index, m.min(axis=1), m.max(axis=1), color=colors[cand], alpha=0.15, linewidth=0)
                ax.plot(m.index, m.mean(axis=1), color=colors[cand], linestyle=styles[cand], linewidth=1.8,
                        label=f"{cand} (seed mean, range)" if r == 0 else None)
            tg = targets[targets.target_reef == reef]
            ax.scatter(tg.report_year, tg.target_density_ha, color="#0b0b0b", s=22, zorder=5,
                       label="diver/cull target" if r == 0 else None)
            ax.set_yscale("symlog", linthresh=1.0)
            ax.grid(color="#e4e3df", linewidth=.8)
            ax.spines[["top", "right"]].set_visible(False)
            if col == 0:
                ax.set_ylabel(f"{reef}\nadults / ha")
        axes[0, col].set_title(f"supply: {ds} outbreak layer", loc="left", fontsize=10)
        axes[0, col].legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(paths["recheck.png"], dpi=160, facecolor="#fcfcfb")
    plt.close(fig)

    report = {"status": "bounded_diagnostic_not_promoted", "reference_replay_max_rel_delta": replay,
              "decision_rules": {"gate": f"G2' at A_half 30 passes in >= {GATE_MIN_PASSES}/5 seeds on both layers",
                                 "ranking": f"seed-mean v3 (averaged over layers) within {V3_MARGIN} of the best gated candidate"},
              "gate_ok": gate_ok.to_dict(), "v3_mean_over_layers": best.to_dict(), "retained": retained,
              "correction_effect_seed1_W": correction_effect,
              "onset_definition": f"first year in the wave window with adults >= {ONSET}/ha", "waves": WAVES,
              "summary": json.loads(summary.to_json(orient="records")),
              "per_arm": json.loads(result.to_json(orient="records")),
              "inputs": {f"{ds}/{p.name}/{n}": sha256(p / n) for ds, (run, rescore, gate) in runs.items()
                         for p, n in ((run, "metadata.toml"), (run, "trajectories.csv"),
                                      (rescore, "detectability_scores.csv"), (gate, "gate_g2prime.csv"))}
                        | {"targets_reef_year.csv": sha256(TARGETS)},
              "script_sha256": sha256(Path(__file__)),
              "outputs": {n: sha256(p) for n, p in paths.items() if n != "recheck_audit.json"}}
    paths["recheck_audit.json"].write_text(json.dumps(report, indent=2, sort_keys=True, default=str) + "\n")
    pd.set_option("display.width", 320)
    pd.set_option("display.max_columns", 80)
    print(result[["dataset", "candidate", "seed", "loss_v3", "loss_v2_manta", "diver_error", "matched_v2",
                  "g2_pass_ahalf30", "g2_pass_ahalf10"]].round(2).to_string(index=False))
    print(summary.round(2).T.to_string())
    print("Correction effect (W, seed 1):", json.dumps(correction_effect, indent=1))
    print("Gate ok:", gate_ok.to_dict(), "| v3 mean:", best.round(2).to_dict(), "| retained:", retained)
    print("Replay max:", max(replay.values()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("run_original", "rescore_original", "gate_original", "run_corrected", "rescore_corrected", "gate_corrected"):
        parser.add_argument(name, type=Path)
    a = parser.parse_args()
    main({"original": (a.run_original.resolve(), a.rescore_original.resolve(), a.gate_original.resolve()),
          "corrected": (a.run_corrected.resolve(), a.rescore_corrected.resolve(), a.gate_corrected.resolve())})
