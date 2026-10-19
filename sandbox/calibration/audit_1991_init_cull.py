"""Audit the starting density x observed culling test (crash_mechanism_protocol.md, "Preregistered
starting density x observed culling", 2026-10-15, fifteenth round): W with the start inversion
{lin, h10, h30} x culling {off, sigma 210 m} x five seed pairs x both outbreak layers.

    python audit_1991_init_cull.py <run_original> <rescore_original> <gate_original> \\
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
PRIOR = {"original": "20261015T_28_recheck_original_1991_v1", "corrected": "20261015T_29_recheck_corrected_1991_v1"}
WAVES = {"w3": (1993, 2003), "w4": (2008, 2015), "w5": (2019, 2026)}
SHORT = {"Lizard Island Reef": "liz", "MacGillivray Reef": "macg", "North Direction Reef": "nd", "Eyrie Reef": "eyr"}
GATED = ("Lizard Island Reef", "MacGillivray Reef", "North Direction Reef")
ONSET = 3.0
GATE_MIN_PASSES, V3_TOLERANCE = 4, 0.5
CULL_WINDOWS = {"2012_16": (2012, 2016), "2022_26": (2022, 2026)}


def parts(arm: str) -> dict:
    return {"start": arm.split("_")[4], "cull": "cull210" if "_cull210_" in arm else "nocull",
            "seed": int(re.search(r"_s(\d+)_", arm).group(1))}


def wave(adults: pd.Series, window: tuple[int, int]) -> dict:
    a = adults.loc[window[0]:window[1]]
    above = a[a >= ONSET]
    return {"onset": int(above.index.min()) if len(above) else np.nan, "peak": float(a.max()),
            "peak_year": int(a.idxmax())}


def main(runs: dict[str, tuple[Path, Path, Path]]) -> None:
    out = runs["corrected"][0]
    paths = {n: out / n for n in ("init_cull_audit.csv", "init_cull_summary.csv", "init_cull_removals.csv",
                                  "init_cull_audit.json", "init_cull.png")}
    if any(p.exists() for p in paths.values()):
        raise FileExistsError("Refusing to overwrite init_cull audit")
    targets = pd.read_csv(TARGETS)
    rows, replay, series, removals = [], {}, {}, []
    observed_first_peak = None
    for ds, (run, rescore, gate) in runs.items():
        tables = {n: pd.read_csv(run / f"{n}.csv") for n in ("trajectories", "flows", "mechanisms", "failures", "initial_reef", "culls")}
        if len(tables["failures"]):
            raise ValueError(f"Failed evaluations in {run.name}")
        controls = sorted(a for a in set(tables["trajectories"].treatment) if parts(a)["start"] == "lin" and parts(a)["cull"] == "nocull")
        for name in ("trajectories", "flows", "mechanisms"):
            prior = pd.read_csv(RUNS / PRIOR[ds] / f"{name}.csv")
            columns = [c for c in prior if c not in ("treatment", "reef_name", "year")]
            for arm in controls:
                replay[f"{ds}/{arm}/{name}"] = exact_replay(tables[name][tables[name].treatment == arm],
                                                            prior[prior.treatment == arm], columns)
        gate_json = json.loads((gate / "gate_g2prime.json").read_text())
        observed_first_peak = {r: v["observed_peaks"][0] for r, v in gate_json["observed_troughs"].items()}
        scores = pd.read_csv(rescore / "detectability_scores.csv")
        g2 = pd.read_csv(gate / "gate_g2prime.csv").set_index(["treatment", "a_half"])
        traj, init = tables["trajectories"], tables["initial_reef"].set_index(["treatment", "reef_name"])
        culls = tables["culls"]
        for arm in sorted(set(traj.treatment)):
            s = scores[(scores.treatment == arm) & (scores.a_half == 0.0)]
            primary = s[s.score_v2_primary.astype(str).str.lower() == "true"].set_index("reef_name")
            row = {"dataset": ds, **parts(arm), "treatment": arm,
                   "loss_v3": float(primary.loss_v3.sum()), "loss_v2_manta": float(primary.loss_v2_manta.sum()),
                   "diver_error": float(primary.diver_error.sum()), "matched_v2": int(primary.matched_peaks.sum()),
                   "g2_pass_ahalf30": bool(g2.loc[(arm, 30.0)].passes_g2prime),
                   "g2_pass_ahalf10": bool(g2.loc[(arm, 10.0)].passes_g2prime)}
            peak_errors = []
            for reef, short in SHORT.items():
                a = traj[(traj.treatment == arm) & (traj.reef_name == reef)].set_index("year").adults_ha
                series[(ds, arm, reef)] = a
                row[f"{short}_start_adults_ha"] = float(init.loc[(arm, reef), "seed_adults_ha"])
                row[f"{short}_diver_error"] = float(primary.loc[reef, "diver_error"])
                for w, window in WAVES.items():
                    for k, v in wave(a, window).items():
                        row[f"{short}_{w}_{k}"] = v
                if reef in GATED:
                    peak_errors.append(abs(row[f"{short}_w3_peak_year"] - observed_first_peak[reef]))
                tg = targets[targets.target_reef == reef].set_index("report_year").target_density_ha
                for y in (2025, 2026):
                    row[f"{short}_{y}_model"] = float(a.get(y, np.nan))
                    row[f"{short}_{y}_target"] = float(tg.get(y, np.nan))
            row["first_peak_year_abs_error"] = float(np.mean(peak_errors))
            rows.append(row)
            c = culls[culls.treatment == arm]
            for scope, g in c.groupby("scope"):
                for wname, (y0, y1) in CULL_WINDOWS.items():
                    w = g[g.year.between(y0, y1)]
                    req = float(w.filter(like="requested_").to_numpy().sum())
                    rem = float(w.filter(like="removed_").to_numpy().sum())
                    removals.append({"dataset": ds, **parts(arm), "scope": scope, "window": wname, "requested": req,
                                     "removed": rem, "unmet_fraction": 1 - rem / req if req > 0 else np.nan})
    result = pd.DataFrame(rows).sort_values(["dataset", "start", "cull", "seed"]).reset_index(drop=True)
    result.to_csv(paths["init_cull_audit.csv"], index=False)
    removals = pd.DataFrame(removals)
    removals_summary = removals[removals.cull == "cull210"].groupby(["dataset", "start", "scope", "window"])[
        ["requested", "removed", "unmet_fraction"]].mean().reset_index()
    removals_summary.to_csv(paths["init_cull_removals.csv"], index=False)

    keys = ["dataset", "start", "cull"]
    summary = result.groupby(keys).agg(
        v3_mean=("loss_v3", "mean"), v3_min=("loss_v3", "min"), v3_max=("loss_v3", "max"),
        manta_mean=("loss_v2_manta", "mean"), diver_mean=("diver_error", "mean"),
        first_peak_year_abs_error=("first_peak_year_abs_error", "mean"),
        g2_passes_ahalf30=("g2_pass_ahalf30", "sum"), g2_passes_ahalf10=("g2_pass_ahalf10", "sum"),
        n_seeds=("seed", "nunique")).reset_index()
    grouped = result.groupby(keys)
    for short in SHORT.values():
        for col in ("start_adults_ha", "w3_onset", "w3_peak_year", "w3_peak", "w4_peak", "w5_peak"):
            summary[f"{short}_{col}_mean"] = grouped[f"{short}_{col}"].mean().values
        for y in (2025, 2026):
            ratio = np.log((result[f"{short}_{y}_model"] + 0.5) / (result[f"{short}_{y}_target"] + 0.5))
            summary[f"{short}_{y}_model_over_target_gmean"] = np.exp(ratio.groupby([result[k] for k in keys]).mean()).values
    summary.to_csv(paths["init_cull_summary.csv"], index=False)

    # Decision rule 1 (start), within each culling level.
    start_rule = {}
    for cull in ("nocull", "cull210"):
        base = summary[(summary.start == "lin") & (summary.cull == cull)].set_index("dataset")
        for start in ("h10", "h30"):
            arm = summary[(summary.start == start) & (summary.cull == cull)].set_index("dataset")
            start_rule[f"{start}/{cull}"] = {
                "g2_ok": bool((arm.g2_passes_ahalf30 >= GATE_MIN_PASSES).all()),
                "d_v3_mean": float((arm.v3_mean - base.v3_mean).mean()),
                "v3_ok": bool((arm.v3_mean - base.v3_mean).mean() <= V3_TOLERANCE),
                "d_first_peak_year_abs_error": float((arm.first_peak_year_abs_error - base.first_peak_year_abs_error).mean()),
                "timing_ok": bool((arm.first_peak_year_abs_error - base.first_peak_year_abs_error).mean() < 0)}
            start_rule[f"{start}/{cull}"]["recommend"] = all(start_rule[f"{start}/{cull}"][k] for k in ("g2_ok", "v3_ok", "timing_ok"))
    cull_effect = {}
    for start in ("lin", "h10", "h30"):
        for ds in runs:
            off = summary[(summary.start == start) & (summary.cull == "nocull") & (summary.dataset == ds)].iloc[0]
            on = summary[(summary.start == start) & (summary.cull == "cull210") & (summary.dataset == ds)].iloc[0]
            cull_effect[f"{start}/{ds}"] = {"d_v3_mean": float(on.v3_mean - off.v3_mean),
                                            "d_diver_mean": float(on.diver_mean - off.diver_mean),
                                            "d_manta_mean": float(on.manta_mean - off.manta_mean),
                                            **{f"{s}_2026_ratio_{k}": float(r[f"{s}_2026_model_over_target_gmean"])
                                               for s in ("liz", "nd", "eyr") for k, r in (("off", off), ("on", on))}}

    colors = {"lin": "#0d366b", "h10": "#3987e5", "h30": "#eb6834"}
    fig, axes = plt.subplots(4, 2, figsize=(13, 13), sharex=True)
    for col, ds in enumerate(("original", "corrected")):
        for r, reef in enumerate(REEFS):
            ax = axes[r, col]
            for start in ("lin", "h10", "h30"):
                for cull, style in (("nocull", "-"), ("cull210", ":")):
                    arms = result[(result.dataset == ds) & (result.start == start) & (result.cull == cull)].treatment
                    m = pd.concat([series[(ds, a, reef)] for a in arms], axis=1)
                    if cull == "nocull":
                        ax.fill_between(m.index, m.min(axis=1), m.max(axis=1), color=colors[start], alpha=0.12, linewidth=0)
                    ax.plot(m.index, m.mean(axis=1), color=colors[start], linestyle=style, linewidth=1.7,
                            label=f"start {start}, {cull}" if r == 0 else None)
            tg = targets[targets.target_reef == reef]
            ax.scatter(tg.report_year, tg.target_density_ha, color="#0b0b0b", s=22, zorder=5,
                       label="diver/cull target" if r == 0 else None)
            if reef in observed_first_peak:
                ax.axvline(observed_first_peak[reef], color="#52514e", linewidth=0.8, linestyle="--")
            ax.axhline(ONSET, color="#52514e", linewidth=0.6, linestyle=":")
            ax.set_yscale("symlog", linthresh=1.0)
            ax.grid(color="#e4e3df", linewidth=.8)
            ax.spines[["top", "right"]].set_visible(False)
            if col == 0:
                ax.set_ylabel(f"{reef}\nadults / ha")
        axes[0, col].set_title(f"supply: {ds} outbreak layer (seed means; dashed = observed first peak)", loc="left", fontsize=10)
        axes[0, col].legend(fontsize=7, frameon=False, ncol=2)
    fig.tight_layout()
    fig.savefig(paths["init_cull.png"], dpi=160, facecolor="#fcfcfb")
    plt.close(fig)

    report = {"status": "bounded_diagnostic_not_promoted", "reference_replay_max_rel_delta": replay,
              "observed_first_peak": observed_first_peak, "start_rule": start_rule, "cull_effect": cull_effect,
              "decision_rules": {"start": f"G2' A_half 30 >= {GATE_MIN_PASSES}/5 seeds on both layers; seed-mean v3 no more than "
                                          f"{V3_TOLERANCE} worse than lin; lower mean |first-peak year error| (Lizard, MacGillivray, North Direction)",
                                 "cull": "report; unmet fraction > 0.25 on a calibration reef limits interpretation"},
              "onset_definition": f"first year in the wave window with adults >= {ONSET}/ha", "waves": WAVES,
              "summary": json.loads(summary.to_json(orient="records")),
              "removals": json.loads(removals_summary.to_json(orient="records")),
              "inputs": {f"{ds}/{p.name}/{n}": sha256(p / n) for ds, (run, rescore, gate) in runs.items()
                         for p, n in ((run, "metadata.toml"), (run, "trajectories.csv"), (run, "culls.csv"),
                                      (rescore, "detectability_scores.csv"), (gate, "gate_g2prime.csv"))}
                        | {"targets_reef_year.csv": sha256(TARGETS)},
              "script_sha256": sha256(Path(__file__)),
              "outputs": {n: sha256(p) for n, p in paths.items() if n != "init_cull_audit.json"}}
    paths["init_cull_audit.json"].write_text(json.dumps(report, indent=2, sort_keys=True, default=str) + "\n")
    pd.set_option("display.width", 320)
    pd.set_option("display.max_columns", 80)
    print(summary.round(2).T.to_string())
    print(removals_summary[removals_summary.scope.isin(["domain", *SHORT])].round(2).to_string(index=False))
    print("Start rule:", json.dumps(start_rule, indent=1))
    print("Cull effect:", json.dumps(cull_effect, indent=1))
    print("Replay max:", max(replay.values()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("run_original", "rescore_original", "gate_original", "run_corrected", "rescore_corrected", "gate_corrected"):
        parser.add_argument(name, type=Path)
    a = parser.parse_args()
    main({"original": (a.run_original.resolve(), a.rescore_original.resolve(), a.gate_original.resolve()),
          "corrected": (a.run_corrected.resolve(), a.rescore_corrected.resolve(), a.gate_corrected.resolve())})
