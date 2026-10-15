"""Audit the external-multiplier x supply x chlorophyll factorial (crash_mechanism_protocol.md,
"Preregistered external-multiplier x supply x chlorophyll factorial", 2026-10-15).

    python audit_1991_ext_chl.py <run_original> <rescore_original> <gate_original> \\
        <run_corrected> <rescore_corrected> <gate_corrected>
Outputs go to <run_corrected>.
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
from audit_1991_initiation import exact_replay  # noqa: E402
from audit_1991_switching_repair import REEFS, sha256  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / "sandbox/calibration/runs"
TARGETS = ROOT / "sandbox/data/derived/density_targets_2020s/targets_reef_year.csv"
PRIOR = {"original": "20261015T_08_chl_survival_owen_global_2018_2023_1991_v1",
         "corrected": "20261015T_09_chl_survival_owen_global_2018_2023_corrected_1991_v1"}
BASE = "wc_b100_q20_eb_lin_obs015_s1_t60"
REPLAYS = {f"{BASE}_e50_chloff_cv2_pnd": f"{BASE}_chloff_cv2_pnd", f"{BASE}_e50_pb30_cv2_pnd": f"{BASE}_pb30_cv2_pnd"}
# ext_chl2 (full chl series): the chl-off arms replay the same-named ext_chl arms (runs 14 / 15).
PRIOR2 = {"original": "20261015T_14_ext_chl_original_1991_v1", "corrected": "20261015T_15_ext_chl_corrected_1991_v1"}
REPLAYS2 = {f"{BASE}_e{e}_chloff_cv2_pnd": f"{BASE}_e{e}_chloff_cv2_pnd" for e in ("10", "25", "50")}
WAVES = {"w3": (1993, 2001), "w4": (2008, 2015), "w5": (2019, 2026)}
CHL_WINDOWS = {"1996_97": (1996, 1997), "2011_12": (2011, 2012), "2019_20": (2019, 2020), "2024_25": (2024, 2025)}
SHORT = {"Lizard Island Reef": "liz", "MacGillivray Reef": "macg", "North Direction Reef": "nd", "Eyrie Reef": "eyr"}
ONSET = 3.0


def parts(arm: str) -> tuple[float, str]:
    p = arm.split("_")
    return int(p[8][1:]) / 10, p[9]


def wave(adults: pd.Series, window: tuple[int, int]) -> dict:
    a = adults.loc[window[0]:window[1]]
    above = a[a >= ONSET]
    onset = int(above.index.min()) if len(above) else np.nan
    return {"onset": onset, "onset_at_window_start": bool(len(above) and above.index.min() == window[0]),
            "peak": float(a.max()), "peak_year": int(a.idxmax())}


def main(runs: dict[str, tuple[Path, Path, Path]]) -> None:
    out = runs["corrected"][0]
    paths = {n: out / n for n in ("ext_chl_audit.csv", "ext_chl_audit.json", "ext_chl.png")}
    if any(p.exists() for p in paths.values()):
        raise FileExistsError("Refusing to overwrite ext-chl audit")
    rows, replay, series = [], {}, {}
    for ds, (run, rescore, gate) in runs.items():
        tables = {n: pd.read_csv(run / f"{n}.csv") for n in ("trajectories", "flows", "mechanisms", "failures")}
        if len(tables["failures"]):
            raise ValueError(f"Failed evaluations in {run.name}")
        second = any(a.split("_")[9] in ("amols", "aprma", "abrma") for a in set(tables["trajectories"].treatment))
        for current, prior_arm in (REPLAYS2 if second else REPLAYS).items():
            for name in ("trajectories", "flows", "mechanisms"):
                prior = pd.read_csv(RUNS / (PRIOR2 if second else PRIOR)[ds] / f"{name}.csv")
                columns = [c for c in prior if c not in ("treatment", "reef_name", "year")]
                replay[f"{ds}/{current}/{name}"] = exact_replay(tables[name][tables[name].treatment == current],
                                                                prior[prior.treatment == prior_arm], columns)
        scores = pd.read_csv(rescore / "detectability_scores.csv")
        g2 = pd.read_csv(gate / "gate_g2prime.csv").set_index(["treatment", "a_half"])
        traj, flows = tables["trajectories"], tables["flows"]
        for arm in sorted(set(traj.treatment)):
            ext, chl = parts(arm)
            s = scores[(scores.treatment == arm) & (scores.a_half == 0.0)]
            primary = s[s.score_v2_primary.astype(str).str.lower() == "true"].set_index("reef_name")
            row = {"dataset": ds, "ext": ext, "chl": chl, "treatment": arm,
                   "loss_v3": float(primary.loss_v3.sum()), "loss_v2_manta": float(primary.loss_v2_manta.sum()),
                   "diver_error": float(primary.diver_error.sum()), "loss_v2_1992_2024": float(primary.loss_v2.sum()),
                   "matched_v2": int(primary.matched_peaks.sum()),
                   "g2_pass_ahalf30": bool(g2.loc[(arm, 30.0)].passes_g2prime),
                   "g2_pass_ahalf10": bool(g2.loc[(arm, 10.0)].passes_g2prime)}
            for reef, short in SHORT.items():
                a = traj[(traj.treatment == arm) & (traj.reef_name == reef)].set_index("year").adults_ha
                f = flows[(flows.treatment == arm) & (flows.reef_name == reef)].set_index("year")
                series[(ds, arm, reef)] = a
                row[f"{short}_diver_error"] = float(primary.loc[reef, "diver_error"])
                for w, window in WAVES.items():
                    for k, v in wave(a, window).items():
                        row[f"{short}_{w}_{k}"] = v
                for k, (y0, y1) in CHL_WINDOWS.items():
                    row[f"{short}_recruits_{k}"] = float(f.loc[y0:y1].settled_recruits.sum())
                for k, (y0, y1) in (("2010_14", (2010, 2014)), ("2019_26", (2019, 2026))):
                    settled = f.loc[y0:y1].settled_recruits.sum()
                    row[f"{short}_local_share_{k}"] = float(1 - f.loc[y0:y1].external_immigration.sum() / settled) if settled > 0 else np.nan
            rows.append(row)
    result = pd.DataFrame(rows).sort_values(["dataset", "ext", "chl"]).reset_index(drop=True)
    # Chlorophyll effect: pb30 minus chloff at the same dataset and Ext.
    effects = []
    for (ds, ext), g in result.groupby(["dataset", "ext"]):
      off = g[g.chl == "chloff"].iloc[0]
      for _, on in g[g.chl != "chloff"].iterrows():
        e = {"dataset": ds, "ext": ext, "chl": on.chl, "d_loss_v3": on.loss_v3 - off.loss_v3,
             "d_loss_v2_manta": on.loss_v2_manta - off.loss_v2_manta, "d_diver_error": on.diver_error - off.diver_error}
        for short in SHORT.values():
            for w in ("w3",):
                e[f"{short}_{w}_onset_shift"] = on[f"{short}_{w}_onset"] - off[f"{short}_{w}_onset"]
                e[f"{short}_{w}_peak_pct"] = float(100 * (on[f"{short}_{w}_peak"] / off[f"{short}_{w}_peak"] - 1))
            for k in CHL_WINDOWS:
                b = off[f"{short}_recruits_{k}"]
                e[f"{short}_recruits_{k}_pct"] = float(100 * (on[f"{short}_recruits_{k}"] / b - 1)) if b > 0 else np.nan
            for w in ("w4", "w5"):
                e[f"{short}_{w}_onset_shift"] = on[f"{short}_{w}_onset"] - off[f"{short}_{w}_onset"]
                e[f"{short}_{w}_peak_pct"] = float(100 * (on[f"{short}_{w}_peak"] / off[f"{short}_{w}_peak"] - 1))
        effects.append(e)
    effects = pd.DataFrame(effects)
    result.to_csv(paths["ext_chl_audit.csv"], index=False)

    targets = pd.read_csv(TARGETS)
    fig, axes = plt.subplots(4, 2, figsize=(13, 13), sharex=True)
    colors = {1.0: "#86b6ef", 2.5: "#3987e5", 5.0: "#0d366b"}
    for col, ds in enumerate(("original", "corrected")):
        for r, reef in enumerate(REEFS):
            ax = axes[r, col]
            for (d, arm, rf), a in series.items():
                if d != ds or rf != reef:
                    continue
                ext, chl = parts(arm)
                if chl not in ("pb30", "aprma", "chloff"):
                    continue
                ax.plot(a.index, a.values, color=colors[ext], linestyle=":" if chl == "chloff" else "-",
                        linewidth=1.3 if chl == "chloff" else 1.8, label=f"Ext x{ext:g}, chl {chl}" if r == 0 else None)
            tg = targets[targets.target_reef == reef]
            ax.scatter(tg.report_year, tg.target_density_ha, color="#eb6834", s=24, zorder=5, label="diver/cull target" if r == 0 else None)
            ax.set_yscale("symlog", linthresh=1.0)
            ax.grid(color="#e4e3df", linewidth=.8)
            if col == 0:
                ax.set_ylabel(f"{reef}\nadults / ha")
        axes[0, col].set_title(f"supply: {ds} outbreak layer (solid chl on, dotted off)", loc="left", fontsize=10)
        axes[0, col].legend(fontsize=7, frameon=False, ncol=2)
    fig.tight_layout()
    fig.savefig(paths["ext_chl.png"], dpi=160, facecolor="#fcfcfb")
    plt.close(fig)

    report = {"status": "bounded_diagnostic_not_promoted", "reference_replay_max_rel_delta": replay,
              "onset_definition": f"first year in the wave window with adults >= {ONSET}/ha", "waves": WAVES,
              "summary": json.loads(result.to_json(orient="records")),
              "chl_effects": json.loads(effects.to_json(orient="records")),
              "inputs": {f"{ds}/{p.name}/{n}": sha256(p / n) for ds, (run, rescore, gate) in runs.items()
                         for p, n in ((run, "metadata.toml"), (run, "trajectories.csv"),
                                      (rescore, "detectability_scores.csv"), (gate, "gate_g2prime.csv"))}
                        | {"targets_reef_year.csv": sha256(TARGETS)},
              "script_sha256": sha256(Path(__file__)),
              "outputs": {n: sha256(p) for n, p in paths.items() if n != "ext_chl_audit.json"}}
    paths["ext_chl_audit.json"].write_text(json.dumps(report, indent=2, sort_keys=True, default=str) + "\n")
    pd.set_option("display.width", 320)
    pd.set_option("display.max_columns", 80)
    print(result[["dataset", "ext", "chl", "loss_v3", "loss_v2_manta", "diver_error", "loss_v2_1992_2024", "matched_v2",
                  "g2_pass_ahalf30", "g2_pass_ahalf10"]].round(2).to_string(index=False))
    for short in SHORT.values():
        cols = [c for c in result.columns if c.startswith(f"{short}_w") and not c.endswith("at_window_start")] + \
               [f"{short}_local_share_2010_14", f"{short}_local_share_2019_26"]
        print(result[["dataset", "ext", "chl"] + cols].round(2).to_string(index=False))
    print(effects.round(1).to_string(index=False))
    print("Replay max:", max(replay.values()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("run_original", "rescore_original", "gate_original", "run_corrected", "rescore_corrected", "gate_corrected"):
        parser.add_argument(name, type=Path)
    a = parser.parse_args()
    main({"original": (a.run_original.resolve(), a.rescore_original.resolve(), a.gate_original.resolve()),
          "corrected": (a.run_corrected.resolve(), a.rescore_corrected.resolve(), a.gate_corrected.resolve())})
