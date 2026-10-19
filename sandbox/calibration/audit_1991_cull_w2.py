"""Audit the cull re-check on the w 2 layer (crash_mechanism_protocol.md, "Preregistered cull
re-check on the w 2 layer", 2026-10-16): W + h10, Q 20 x culling {off, sigma 210 m} x five seeds.

    python audit_1991_cull_w2.py <run> <rescore> <gate>
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
PRIOR = "20261016T_02_pit_q_w2_1991_v1"
ARM = re.compile(r"^wc_b100_q20_eb_h10_obs015_s(\d)_t60_e50_(cull210_)?aprmak_cv2_pnd$")
SHORT = {"Lizard Island Reef": "liz", "MacGillivray Reef": "macg", "North Direction Reef": "nd", "Eyrie Reef": "eyr"}
YEARS = (2022, 2023, 2024, 2025, 2026)


def main(run: Path, rescore: Path, gate: Path) -> None:
    paths = {n: run / n for n in ("cull_w2_audit.csv", "cull_w2_summary.csv", "cull_w2_removals.csv",
                                  "cull_w2_audit.json", "cull_w2.png")}
    if any(p.exists() for p in paths.values()):
        raise FileExistsError("Refusing to overwrite cull_w2 audit")
    t = {n: pd.read_csv(run / f"{n}.csv") for n in ("trajectories", "flows", "mechanisms", "failures", "culls")}
    if len(t["failures"]):
        raise ValueError("Failed evaluations")
    arms = sorted(set(t["trajectories"].treatment))
    if len(arms) != 10 or not all(ARM.match(a) for a in arms):
        raise ValueError("Unexpected arm set")
    replay = {}
    for name in ("trajectories", "flows", "mechanisms"):
        prior = pd.read_csv(RUNS / PRIOR / f"{name}.csv")
        cols = [c for c in prior if c not in ("treatment", "reef_name", "year")]
        for arm in (a for a in arms if not ARM.match(a).group(2)):
            replay[f"{arm}/{name}"] = exact_replay(t[name][t[name].treatment == arm], prior[prior.treatment == arm], cols)
    targets = pd.read_csv(TARGETS)
    scores = pd.read_csv(rescore / "detectability_scores.csv")
    g2 = pd.read_csv(gate / "gate_g2prime.csv").set_index(["treatment", "a_half"])
    rows, series, removals = [], {}, []
    for arm in arms:
        m = ARM.match(arm)
        cull = "cull210" if m.group(2) else "nocull"
        s = scores[(scores.treatment == arm) & (scores.a_half == 0.0)]
        primary = s[s.score_v2_primary.astype(str).str.lower() == "true"].set_index("reef_name")
        row = {"cull": cull, "seed": int(m.group(1)), "loss_v3": float(primary.loss_v3.sum()),
               "loss_v2_manta": float(primary.loss_v2_manta.sum()), "diver_error": float(primary.diver_error.sum()),
               "g2_pass_ahalf30": bool(g2.loc[(arm, 30.0)].passes_g2prime)}
        for reef, short in SHORT.items():
            a = t["trajectories"][(t["trajectories"].treatment == arm) & (t["trajectories"].reef_name == reef)].set_index("year").adults_ha
            series[(cull, int(m.group(1)), reef)] = a
            row[f"{short}_diver_error"] = float(primary.loc[reef, "diver_error"])
            tg = targets[targets.target_reef == reef].set_index("report_year").target_density_ha
            for y in YEARS:
                row[f"{short}_{y}_model"] = float(a.get(y, np.nan))
                row[f"{short}_{y}_target"] = float(tg.get(y, np.nan))
        rows.append(row)
        c = t["culls"][t["culls"].treatment == arm]
        for scope, g in c.groupby("scope"):
            for wname, (y0, y1) in {"2012_16": (2012, 2016), "2022_26": (2022, 2026)}.items():
                w = g[g.year.between(y0, y1)]
                req, rem = float(w.filter(like="requested_").to_numpy().sum()), float(w.filter(like="removed_").to_numpy().sum())
                removals.append({"cull": cull, "seed": int(m.group(1)), "scope": scope, "window": wname,
                                 "requested": req, "removed": rem, "unmet_fraction": 1 - rem / req if req > 0 else np.nan})
    result = pd.DataFrame(rows).sort_values(["cull", "seed"])
    result.to_csv(paths["cull_w2_audit.csv"], index=False)
    rem = pd.DataFrame(removals)
    rem_summary = rem[rem.cull == "cull210"].groupby(["scope", "window"])[["requested", "removed", "unmet_fraction"]].mean().reset_index()
    rem_summary.to_csv(paths["cull_w2_removals.csv"], index=False)
    g = result.groupby("cull")
    summary = g.agg(v3_mean=("loss_v3", "mean"), v3_min=("loss_v3", "min"), v3_max=("loss_v3", "max"),
                    manta_mean=("loss_v2_manta", "mean"), diver_mean=("diver_error", "mean"),
                    g2_passes_ahalf30=("g2_pass_ahalf30", "sum"))
    for short in SHORT.values():
        summary[f"{short}_diver_error_mean"] = g[f"{short}_diver_error"].mean()
        for y in YEARS:
            summary[f"{short}_{y}_model_mean"] = g[f"{short}_{y}_model"].mean()
            summary[f"{short}_{y}_target"] = g[f"{short}_{y}_target"].first()
    summary = summary.reset_index()
    summary.to_csv(paths["cull_w2_summary.csv"], index=False)
    paired = result.pivot(index="seed", columns="cull", values="loss_v3")

    fig, axes = plt.subplots(2, 2, figsize=(13, 8.5), sharex=True)
    for ax, reef in zip(axes.flat, REEFS):
        for cull, color in (("nocull", "#86b6ef"), ("cull210", "#0d366b")):
            m = pd.concat([series[(cull, s, reef)] for s in range(1, 6)], axis=1).loc[2012:]
            ax.fill_between(m.index, m.min(axis=1), m.max(axis=1), color=color, alpha=0.15, linewidth=0)
            ax.plot(m.index, m.mean(axis=1), color=color, linewidth=1.9, label="observed culling" if cull == "cull210" else "no culling")
        tg = targets[(targets.target_reef == reef) & (targets.report_year >= 2012)]
        ax.scatter(tg.report_year, tg.target_density_ha, color="#eb6834", s=24, zorder=5, label="diver/cull target")
        ax.set_yscale("symlog", linthresh=1.0)
        ax.set_ylim(bottom=0)
        ax.set_title(reef, loc="left", fontsize=10)
        ax.set_ylabel("adults / ha")
        ax.grid(color="#e4e3df", linewidth=.8)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0, 0].legend(frameon=False, fontsize=8)
    fig.suptitle("W + moderate start, correction weight 2: culling off vs on (seed mean, range)", x=0.01, ha="left", fontsize=11)
    fig.tight_layout()
    fig.savefig(paths["cull_w2.png"], dpi=150, facecolor="#fcfcfb")
    plt.close(fig)

    report = {"status": "diagnostic_not_adopted", "reference_replay_max_rel_delta": replay,
              "summary": json.loads(summary.to_json(orient="records")),
              "removals": json.loads(rem_summary.to_json(orient="records")),
              "paired_v3_by_seed": json.loads(paired.to_json(orient="index")),
              "inputs": {f"{p.name}/{n}": sha256(p / n) for p, n in ((run, "metadata.toml"), (run, "trajectories.csv"),
                                                                     (run, "culls.csv"), (rescore, "detectability_scores.csv"),
                                                                     (gate, "gate_g2prime.csv"))}
                        | {"targets_reef_year.csv": sha256(TARGETS)},
              "script_sha256": sha256(Path(__file__)),
              "outputs": {n: sha256(p) for n, p in paths.items() if n != "cull_w2_audit.json"}}
    paths["cull_w2_audit.json"].write_text(json.dumps(report, indent=2, sort_keys=True, default=str) + "\n")
    pd.set_option("display.width", 250)
    print(summary.set_index("cull").T.round(3).to_string())
    print(rem_summary[rem_summary.scope.isin(["domain", *SHORT])].round(2).to_string(index=False))
    print(paired.round(2).to_string())
    print("Replay max:", max(replay.values()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("run", "rescore", "gate"):
        parser.add_argument(name, type=Path)
    a = parser.parse_args()
    main(a.run.resolve(), a.rescore.resolve(), a.gate.resolve())
