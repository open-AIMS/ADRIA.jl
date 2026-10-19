"""Audit the local vs external production diagnostic (crash_mechanism_protocol.md, "Preregistered
local vs external production diagnostic", 2026-10-15, seventeenth round): W + h10 on the face-value
layer x Ext {5, 2.5} x local larval-output multiplier {1, 2, 3} x five seed pairs.

    python audit_1991_local_ext.py <run> <rescore> <gate>
Outputs go to <run>.
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
PRIOR = "20261015T_42_facevalue_h10_1991_v1"
ARM = re.compile(r"^wc_b100_q20_eb_h10_obs015_s(\d)_t60_e(\d\d)(?:_f(\d{3}))?_aprmak_cv2_pnd$")
SHORT = {"Lizard Island Reef": "liz", "MacGillivray Reef": "macg", "North Direction Reef": "nd", "Eyrie Reef": "eyr"}
WAVES = {"w3": (1993, 2003), "w4": (2008, 2015), "w5": (2019, 2026)}
PERIODS = {"1992_2000": (1992, 2000), "2001_2011": (2001, 2011), "2012_2016": (2012, 2016), "2019_2026": (2019, 2026)}
ONSET = 3.0


def parts(arm: str) -> dict:
    m = ARM.match(arm)
    return {"seed": int(m.group(1)), "ext": int(m.group(2)) / 10, "fec": int(m.group(3) or 100) / 100}


def wave(adults: pd.Series, window: tuple[int, int]) -> dict:
    a = adults.loc[window[0]:window[1]]
    above = a[a >= ONSET]
    return {"onset": int(above.index.min()) if len(above) else np.nan, "peak": float(a.max()),
            "peak_year": int(a.idxmax())}


def main(run: Path, rescore: Path, gate: Path) -> None:
    paths = {n: run / n for n in ("local_ext_audit.csv", "local_ext_summary.csv", "local_ext_audit.json", "local_ext.png")}
    if any(p.exists() for p in paths.values()):
        raise FileExistsError("Refusing to overwrite local_ext audit")
    tables = {n: pd.read_csv(run / f"{n}.csv") for n in ("trajectories", "flows", "mechanisms", "failures")}
    if len(tables["failures"]):
        raise ValueError(f"Failed evaluations in {run.name}")
    arms = sorted(set(tables["trajectories"].treatment))
    if len(arms) != 30 or not all(ARM.match(a) for a in arms):
        raise ValueError("Unexpected arm set")
    replay = {}
    controls = [a for a in arms if parts(a)["ext"] == 5.0 and parts(a)["fec"] == 1.0]
    for name in ("trajectories", "flows", "mechanisms"):
        prior = pd.read_csv(RUNS / PRIOR / f"{name}.csv")
        columns = [c for c in prior if c not in ("treatment", "reef_name", "year")]
        for arm in controls:
            replay[f"{arm}/{name}"] = exact_replay(tables[name][tables[name].treatment == arm],
                                                   prior[prior.treatment == arm], columns)
    targets = pd.read_csv(TARGETS)
    scores = pd.read_csv(rescore / "detectability_scores.csv")
    g2 = pd.read_csv(gate / "gate_g2prime.csv").set_index(["treatment", "a_half"])
    observed_first_peak = {r: v["observed_peaks"][0]
                           for r, v in json.loads((gate / "gate_g2prime.json").read_text())["observed_troughs"].items()}
    traj, flows = tables["trajectories"], tables["flows"]
    rows, series = [], {}
    for arm in arms:
        s = scores[(scores.treatment == arm) & (scores.a_half == 0.0)]
        primary = s[s.score_v2_primary.astype(str).str.lower() == "true"].set_index("reef_name")
        row = {**parts(arm), "treatment": arm, "loss_v3": float(primary.loss_v3.sum()),
               "loss_v2_manta": float(primary.loss_v2_manta.sum()), "diver_error": float(primary.diver_error.sum()),
               "matched_v2": int(primary.matched_peaks.sum()),
               "g2_pass_ahalf30": bool(g2.loc[(arm, 30.0)].passes_g2prime)}
        errors = []
        for reef, short in SHORT.items():
            a = traj[(traj.treatment == arm) & (traj.reef_name == reef)].set_index("year").adults_ha
            f = flows[(flows.treatment == arm) & (flows.reef_name == reef)].set_index("year")
            series[(arm, reef)] = a
            row[f"{short}_diver_error"] = float(primary.loc[reef, "diver_error"])
            for w, window in WAVES.items():
                for k, v in wave(a, window).items():
                    row[f"{short}_{w}_{k}"] = v
            if reef in observed_first_peak:
                errors.append(abs(row[f"{short}_w3_peak_year"] - observed_first_peak[reef]))
            for pname, (y0, y1) in PERIODS.items():
                g = f.loc[y0:y1]
                arrivals = (g.external_immigration + g.internal_immigration + g.local_retention + g.background_immigration).sum()
                row[f"{short}_domain_share_{pname}"] = float((g.internal_immigration + g.local_retention).sum() / arrivals) \
                    if arrivals > 0 else np.nan
            tg = targets[targets.target_reef == reef].set_index("report_year").target_density_ha
            for y in range(2020, 2027):
                row[f"{short}_{y}_model"] = float(a.get(y, np.nan))
                row[f"{short}_{y}_target"] = float(tg.get(y, np.nan))
        row["first_peak_year_abs_error"] = float(np.mean(errors))
        rows.append(row)
    result = pd.DataFrame(rows).sort_values(["ext", "fec", "seed"], ascending=[False, True, True]).reset_index(drop=True)
    result.to_csv(paths["local_ext_audit.csv"], index=False)
    g = result.groupby(["ext", "fec"])
    summary = g.agg(v3_mean=("loss_v3", "mean"), v3_min=("loss_v3", "min"), v3_max=("loss_v3", "max"),
                    manta_mean=("loss_v2_manta", "mean"), diver_mean=("diver_error", "mean"),
                    g2_passes_ahalf30=("g2_pass_ahalf30", "sum"),
                    first_peak_year_abs_error=("first_peak_year_abs_error", "mean"))
    for short in SHORT.values():
        for col in ("w3_onset", "w3_peak_year", "w3_peak", "w4_onset", "w4_peak", "w5_onset", "w5_peak",
                    "domain_share_1992_2000", "domain_share_2019_2026", "diver_error"):
            summary[f"{short}_{col}_mean"] = g[f"{short}_{col}"].mean()
        for y in (2022, 2023, 2024, 2025, 2026):
            ratio = np.log((result[f"{short}_{y}_model"] + 0.5) / (result[f"{short}_{y}_target"] + 0.5))
            summary[f"{short}_{y}_model_mean"] = g[f"{short}_{y}_model"].mean()
            summary[f"{short}_{y}_model_over_target_gmean"] = np.exp(ratio.groupby([result.ext, result.fec]).mean())
    summary = summary.reset_index().sort_values(["ext", "fec"], ascending=[False, True])
    summary.to_csv(paths["local_ext_summary.csv"], index=False)

    colors = {1.0: "#86b6ef", 2.0: "#3987e5", 3.0: "#0d366b"}
    fig, axes = plt.subplots(4, 2, figsize=(13, 13), sharex=True)
    for col, ext in enumerate((5.0, 2.5)):
        for r, reef in enumerate(REEFS):
            ax = axes[r, col]
            for fec in (1.0, 2.0, 3.0):
                sel = result[(result.ext == ext) & (result.fec == fec)].treatment
                m = pd.concat([series[(a, reef)] for a in sel], axis=1)
                ax.fill_between(m.index, m.min(axis=1), m.max(axis=1), color=colors[fec], alpha=0.12, linewidth=0)
                ax.plot(m.index, m.mean(axis=1), color=colors[fec], linewidth=1.8, label=f"local output x{fec:g}" if r == 0 else None)
            tg = targets[targets.target_reef == reef]
            ax.scatter(tg.report_year, tg.target_density_ha, color="#eb6834", s=22, zorder=5,
                       label="diver/cull target" if r == 0 else None)
            if reef in observed_first_peak:
                ax.axvline(observed_first_peak[reef], color="#52514e", linewidth=0.8, linestyle="--")
            ax.set_yscale("symlog", linthresh=1.0)
            ax.set_ylim(bottom=0)
            ax.grid(color="#e4e3df", linewidth=.8)
            ax.spines[["top", "right"]].set_visible(False)
            if col == 0:
                ax.set_ylabel(f"{reef}\nadults / ha")
        axes[0, col].set_title(f"external boost x{ext:g} (seed mean, range; dashed = observed first peak)", loc="left", fontsize=10)
        axes[0, col].legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(paths["local_ext.png"], dpi=160, facecolor="#fcfcfb")
    plt.close(fig)

    report = {"status": "labelled_diagnostic_not_adopted", "reference_replay_max_rel_delta": replay,
              "observed_first_peak": observed_first_peak, "waves": WAVES,
              "onset_definition": f"first year in the wave window with adults >= {ONSET}/ha",
              "summary": json.loads(summary.to_json(orient="records")),
              "inputs": {f"{p.name}/{n}": sha256(p / n) for p, n in ((run, "metadata.toml"), (run, "trajectories.csv"),
                                                                     (run, "flows.csv"), (rescore, "detectability_scores.csv"),
                                                                     (gate, "gate_g2prime.csv"))}
                        | {"targets_reef_year.csv": sha256(TARGETS)},
              "script_sha256": sha256(Path(__file__)),
              "outputs": {n: sha256(p) for n, p in paths.items() if n != "local_ext_audit.json"}}
    paths["local_ext_audit.json"].write_text(json.dumps(report, indent=2, sort_keys=True, default=str) + "\n")
    pd.set_option("display.width", 250)
    print(summary.set_index(["ext", "fec"]).T.round(2).to_string())
    print(result.pivot_table(index="seed", columns=["ext", "fec"], values="loss_v3").round(2).to_string())
    print("Replay max:", max(replay.values()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("run", "rescore", "gate"):
        parser.add_argument(name, type=Path)
    a = parser.parse_args()
    main(a.run.resolve(), a.rescore.resolve(), a.gate.resolve())
