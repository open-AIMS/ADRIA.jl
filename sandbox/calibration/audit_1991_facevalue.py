"""Audit option (b), the face-value correction (crash_mechanism_protocol.md, "Preregistered
face-value correction, option (b)", 2026-10-15): W + h10 start, no cull, x 5 seeds on the
face-value layer, against the same arms on the original and corrected layers (runs 36 / 37).

    python audit_1991_facevalue.py <run_fv> <rescore_fv> <gate_fv> \\
        <run_orig> <rescore_orig> <gate_orig> <run_corr> <rescore_corr> <gate_corr>
Outputs go to <run_fv>.
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
from audit_1991_switching_repair import REEFS, sha256  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
TARGETS = ROOT / "sandbox/data/derived/density_targets_2020s/targets_reef_year.csv"
ARM = re.compile(r"^wc_b100_q20_eb_h10_obs015_s(\d)_t60_e50_aprmak_cv2_pnd$")
SHORT = {"Lizard Island Reef": "liz", "MacGillivray Reef": "macg", "North Direction Reef": "nd", "Eyrie Reef": "eyr"}
YEARS = range(2020, 2027)


def main(sets: dict[str, tuple[Path, Path, Path]]) -> None:
    out = sets["face_value"][0]
    paths = {n: out / n for n in ("facevalue_audit.csv", "facevalue_summary.csv", "facevalue_audit.json", "facevalue.png")}
    if any(p.exists() for p in paths.values()):
        raise FileExistsError("Refusing to overwrite face-value audit")
    targets = pd.read_csv(TARGETS)
    rows, series = [], {}
    for layer, (run, rescore, gate) in sets.items():
        traj, flows = pd.read_csv(run / "trajectories.csv"), pd.read_csv(run / "flows.csv")
        if len(pd.read_csv(run / "failures.csv")):
            raise ValueError(f"Failed evaluations in {run.name}")
        scores = pd.read_csv(rescore / "detectability_scores.csv")
        g2 = pd.read_csv(gate / "gate_g2prime.csv").set_index(["treatment", "a_half"])
        arms = sorted(a for a in set(traj.treatment) if ARM.match(a))
        if len(arms) != 5:
            raise ValueError(f"Expected 5 W + h10 no-cull arms in {run.name}, found {len(arms)}")
        for arm in arms:
            s = scores[(scores.treatment == arm) & (scores.a_half == 0.0)]
            primary = s[s.score_v2_primary.astype(str).str.lower() == "true"].set_index("reef_name")
            row = {"layer": layer, "seed": int(ARM.match(arm).group(1)), "loss_v3": float(primary.loss_v3.sum()),
                   "loss_v2_manta": float(primary.loss_v2_manta.sum()), "diver_error": float(primary.diver_error.sum()),
                   "g2_pass_ahalf30": bool(g2.loc[(arm, 30.0)].passes_g2prime)}
            for reef, short in SHORT.items():
                a = traj[(traj.treatment == arm) & (traj.reef_name == reef)].set_index("year").adults_ha
                f = flows[(flows.treatment == arm) & (flows.reef_name == reef)].set_index("year")
                series[(layer, arm, reef)] = a
                row[f"{short}_diver_error"] = float(primary.loc[reef, "diver_error"])
                tg = targets[targets.target_reef == reef].set_index("report_year").target_density_ha
                for y in YEARS:
                    row[f"{short}_{y}_model"] = float(a.get(y, np.nan))
                    row[f"{short}_{y}_target"] = float(tg.get(y, np.nan))
                    row[f"{short}_{y}_external"] = float(f.loc[y, "external_immigration"])
            rows.append(row)
    result = pd.DataFrame(rows).sort_values(["layer", "seed"]).reset_index(drop=True)
    result.to_csv(paths["facevalue_audit.csv"], index=False)
    g = result.groupby("layer")
    summary = g.agg(v3_mean=("loss_v3", "mean"), v3_min=("loss_v3", "min"), v3_max=("loss_v3", "max"),
                    manta_mean=("loss_v2_manta", "mean"), diver_mean=("diver_error", "mean"),
                    g2_passes_ahalf30=("g2_pass_ahalf30", "sum"))
    for short in SHORT.values():
        summary[f"{short}_diver_error_mean"] = g[f"{short}_diver_error"].mean()
        for y in YEARS:
            summary[f"{short}_{y}_model_mean"] = g[f"{short}_{y}_model"].mean()
            summary[f"{short}_{y}_target"] = g[f"{short}_{y}_target"].first()
            summary[f"{short}_{y}_external_mean"] = g[f"{short}_{y}_external"].mean()
    summary = summary.reset_index()
    summary.to_csv(paths["facevalue_summary.csv"], index=False)
    paired = result.pivot(index="seed", columns="layer", values="loss_v3")

    colors = {"original": "#86b6ef", "face_value": "#0d366b", "corrected": "#eb6834"}
    fig, axes = plt.subplots(2, 2, figsize=(13, 8.5), sharex=True)
    for ax, reef in zip(axes.flat, REEFS):
        for layer in ("original", "face_value", "corrected"):
            arms = sorted(a for (l, a, r) in series if l == layer and r == reef)
            m = pd.concat([series[(layer, a, reef)] for a in arms], axis=1).loc[2015:]
            ax.fill_between(m.index, m.min(axis=1), m.max(axis=1), color=colors[layer], alpha=0.15, linewidth=0)
            ax.plot(m.index, m.mean(axis=1), color=colors[layer], linewidth=1.9, label=f"{layer.replace('_', ' ')} layer")
        tg = targets[(targets.target_reef == reef) & (targets.report_year >= 2015)]
        ax.scatter(tg.report_year, tg.target_density_ha, color="#0b0b0b", s=24, zorder=5, label="diver/cull target")
        ax.set_yscale("symlog", linthresh=1.0)
        ax.set_ylim(bottom=0)
        ax.set_title(reef, loc="left", fontsize=10)
        ax.set_ylabel("adults / ha")
        ax.grid(color="#e4e3df", linewidth=.8)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0, 0].legend(frameon=False, fontsize=8)
    fig.suptitle("W + moderate start, 2015-2026: seed mean and range by outbreak layer", x=0.01, ha="left", fontsize=11)
    fig.tight_layout()
    fig.savefig(paths["facevalue.png"], dpi=160, facecolor="#fcfcfb")
    plt.close(fig)

    report = {"status": "diagnostic_not_adopted", "summary": json.loads(summary.to_json(orient="records")),
              "paired_v3_by_seed": json.loads(paired.to_json(orient="index")),
              "inputs": {f"{layer}/{p.name}/{n}": sha256(p / n) for layer, (run, rescore, gate) in sets.items()
                         for p, n in ((run, "metadata.toml"), (run, "trajectories.csv"), (run, "flows.csv"),
                                      (rescore, "detectability_scores.csv"), (gate, "gate_g2prime.csv"))}
                        | {"targets_reef_year.csv": sha256(TARGETS)},
              "script_sha256": sha256(Path(__file__)),
              "outputs": {n: sha256(p) for n, p in paths.items() if n != "facevalue_audit.json"}}
    paths["facevalue_audit.json"].write_text(json.dumps(report, indent=2, sort_keys=True, default=str) + "\n")
    pd.set_option("display.width", 250)
    print(summary.set_index("layer").T.round(2).to_string())
    print(paired.round(2).to_string())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    names = [f"{k}_{l}" for l in ("fv", "orig", "corr") for k in ("run", "rescore", "gate")]
    for name in names:
        parser.add_argument(name, type=Path)
    a = vars(parser.parse_args())
    main({layer: tuple(a[f"{k}_{s}"].resolve() for k in ("run", "rescore", "gate"))
          for layer, s in (("face_value", "fv"), ("original", "orig"), ("corrected", "corr"))})
