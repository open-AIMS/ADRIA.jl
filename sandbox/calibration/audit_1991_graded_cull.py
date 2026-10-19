"""Audit the graded-correction weight x culling factorial (crash_mechanism_protocol.md, "Preregistered:
graded correction weight x culling", 2026-10-19, twenty-sixth round).

Sources (W-wave, five seeds; found by run-name pattern under runs/):
- binary layer (corr_w3): run 20261019T_24 no-cull and _cullreef_ arms (re-score 25, G2' 26);
- graded w3: culling off from 20261019T_30 (31, 32); culling on from *_gradedw3_cull_1991_v1;
- graded w0.5 / w1 / w2: *_gradedw{05,1,2}_cull_1991_v1 (both culling settings).
Identity: every no-cull arm equals run 15 for 1991-2017; every cull arm equals its same-layer no-cull arm to 2001.
Selection (preregistered): per culling setting, lowest seed-mean diver error among weights with G2' >= 4/5 and
dv3 <= +0.5 against the binary layer at the same culling setting.

    python audit_1991_graded_cull.py --output <dir>
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
from audit_1991_switching_repair import sha256  # noqa: E402

RUNS = Path(__file__).resolve().parent / "runs"
CONTROL = "20261019T_15_wave_w3_q25_1991_v1"
TARGETS = Path(__file__).resolve().parents[2] / "sandbox/data/derived/density_targets_2020s/targets_reef_year.csv"
REEFS = ("Lizard Island Reef", "MacGillivray Reef", "North Direction Reef", "Eyrie Reef")
RIBBONS = ("Yonge Reef", "Carter Reef", "No Name Reef")
Q, A_HALF = 0.015, 30.0
LAYER_ORDER = ["binary", "g0.5", "g1", "g2", "g3"]
COLORS = {"binary": "#0d366b", "g0.5": "#86b6ef", "g1": "#3987e5", "g2": "#eb6834", "g3": "#b6b5b0"}


def one(pattern: str) -> Path:
    hits = sorted(RUNS.glob(pattern))
    if len(hits) != 1:
        raise FileNotFoundError(f"Expected one run for {pattern}, found {[h.name for h in hits]}")
    return hits[0]


def sources() -> list[tuple[str, Path, Path, Path]]:
    out = [("binary", one("*_24_cull_place_wave_1991_v1"), one("*_25_rescore_cull_place"), one("*_26_gate_cull_place_lizmacg")),
           ("g3", one("*_30_graded_wave_1991_v1"), one("*_31_rescore_graded"), one("*_32_gate_graded_lizmacg"))]
    for w, label in (("05", "g0.5"), ("1", "g1"), ("2", "g2"), ("3", "g3")):
        out.append((label, one(f"*_gradedw{w}_cull_1991_v1"), one(f"*_rescore_gradedw{w}_cull"), one(f"*_gate_gradedw{w}_cull_lizmacg")))
    return out


def keep(arm: str) -> bool:
    return arm.endswith("_wv025") and "_cull210_" not in arm


def main(out: Path) -> None:
    paths = {n: out / n for n in ("graded_cull_audit.csv", "graded_cull_summary.csv", "graded_cull_audit.json", "graded_cull.png")}
    if any(p.exists() for p in paths.values()):
        raise FileExistsError("Refusing to overwrite audit")
    targets = pd.read_csv(TARGETS)
    control = {n: pd.read_csv(RUNS / CONTROL / f"{n}.csv") for n in ("trajectories", "flows", "mechanisms")}
    rows, identity, series, inputs = [], {}, {}, {}
    for label, run, rescore, gate in sources():
        if len(pd.read_csv(run / "failures.csv")):
            raise ValueError(f"Failed evaluations in {run.name}")
        inputs[run.name] = sha256(run / "trajectories.csv")
        tabs = {n: pd.read_csv(run / f"{n}.csv") for n in ("trajectories", "flows", "mechanisms")}
        extra = pd.read_csv(run / "trajectories_extra.csv")
        culls = pd.read_csv(run / "culls.csv") if (run / "culls.csv").exists() else None
        scores = pd.read_csv(rescore / "detectability_scores.csv")
        g2 = pd.read_csv(gate / "gate_g2prime.csv").set_index(["treatment", "a_half"])
        arms = [a for a in sorted(set(tabs["trajectories"].treatment)) if keep(a)]
        for arm in arms:
            cull = "_cullreef_" in arm
            nocull_arm = arm.replace("_cullreef", "")
            for n, cur in tabs.items():
                cols = [c for c in control[n] if c not in ("treatment", "reef_name", "year")]
                if not cull:
                    identity[f"{label}/{arm}/{n}<=2017"] = exact_replay(cur[(cur.treatment == arm) & (cur.year <= 2017)],
                                                                      control[n][(control[n].treatment == arm) & (control[n].year <= 2017)], cols)
                elif nocull_arm in set(cur.treatment):
                    identity[f"{label}/{arm}/{n}<=2001"] = exact_replay(cur[(cur.treatment == arm) & (cur.year <= 2001)],
                                                                      cur[(cur.treatment == nocull_arm) & (cur.year <= 2001)], cols)
            s = scores[(scores.treatment == arm) & (scores.a_half == 0.0)]
            p = s[s.score_v2_primary.astype(str).str.lower() == "true"].set_index("reef_name")
            row = {"layer": label, "cull": cull, "treatment": arm, "seed": int(re.search(r"_s(\d)_", arm).group(1)),
                   "loss_v3": float(p.loss_v3.sum()), "loss_v2_manta": float(p.loss_v2_manta.sum()),
                   "diver_error": float(p.diver_error.sum()), "g2_pass_ahalf30": bool(g2.loc[(arm, 30.0)].passes_g2prime)}
            traj = tabs["trajectories"]
            for reef in REEFS:
                key = reef.split()[0].lower()
                row[f"{key}_diver"] = float(p.diver_error.get(reef, np.nan))
                a = traj[(traj.treatment == arm) & (traj.reef_name == reef)].set_index("year").adults_ha
                series[(label, cull, reef, arm)] = a
                hit = a.loc[2019:][a.loc[2019:] >= 5.0]
                row[f"{key}_first_ge5"] = int(hit.index.min()) if len(hit) else 2027
                for y in (2021, 2022, 2023, 2024, 2026):
                    row[f"{key}_{y}"] = float(a.get(y, np.nan))
            for reef in RIBBONS:
                a = extra[(extra.treatment == arm) & (extra.reef_name == reef)].set_index("year").adults_ha.loc[2018:]
                row[f"{reef.split()[0].lower()}_max_cpue_2018_26"] = float((Q * a ** 2 / (a + A_HALF)).max())
            if cull and culls is not None:
                c = culls[(culls.treatment == arm) & culls.year.between(2019, 2026)]
                for scope in ("domain", "Eyrie Reef", "Lizard Island Reef"):
                    x = c[c.scope == scope]
                    req, got = x.filter(like="requested").sum().sum(), x.filter(like="removed").sum().sum()
                    row[f"unmet_2019_26_{scope.split()[0].lower()}"] = float(1 - got / req) if req > 0 else np.nan
            rows.append(row)
    d = pd.DataFrame(rows).drop_duplicates("treatment" if False else ["layer", "treatment"])
    out.mkdir(parents=True, exist_ok=True)
    d.to_csv(paths["graded_cull_audit.csv"], index=False)
    num = [c for c in d if c not in ("layer", "cull", "treatment", "seed")]
    summary = d.groupby(["cull", "layer"])[num].mean()
    summary["g2_passes_ahalf30"] = d.groupby(["cull", "layer"]).g2_pass_ahalf30.sum()
    summary["n_seeds"] = d.groupby(["cull", "layer"]).seed.nunique()
    summary = summary.reset_index()
    summary["d_v3_vs_binary"] = summary.loss_v3 - summary.groupby("cull").apply(
        lambda g: g.set_index("layer").loss_v3.get("binary"), include_groups=False).reindex(summary.cull).to_numpy()
    summary["eligible"] = (summary.layer != "binary") & (summary.g2_passes_ahalf30 >= 4) & (summary.d_v3_vs_binary <= 0.5)
    summary["layer"] = pd.Categorical(summary.layer, LAYER_ORDER, ordered=True)
    summary = summary.sort_values(["cull", "layer"])
    summary.to_csv(paths["graded_cull_summary.csv"], index=False)
    choice = {}
    for cull, g in summary.groupby("cull"):
        e = g[g.eligible]
        choice["cull_on" if cull else "cull_off"] = str(e.loc[e.diver_error.idxmin(), "layer"]) if len(e) else None

    fig, axes = plt.subplots(2, 4, figsize=(17, 7.5), sharex=True, sharey="col")
    for r, cull in enumerate((False, True)):
        for c, reef in enumerate(REEFS):
            ax = axes[r, c]
            for label in LAYER_ORDER:
                ss = [v for (lb, cu, rf, _), v in series.items() if lb == label and cu == cull and rf == reef]
                if not ss:
                    continue
                m = pd.concat(ss, axis=1).mean(axis=1)
                ax.plot(m.index, m, color=COLORS[label], linewidth=1.7,
                        label="binary (W-wave)" if label == "binary" else f"graded w{label[1:]}")
            t = targets[(targets.target_reef == reef) & (targets.target_density_ha > 0)]
            ax.scatter(t.report_year, t.target_density_ha, s=16, color="#52514e", zorder=5, label="2020s diver/cull target")
            ax.set_yscale("log")
            ax.set_xlim(2010, 2027)
            ax.set_title(f"{reef} — culling {'on' if cull else 'off'}", loc="left", fontsize=9)
            ax.grid(color="#e4e3df", linewidth=.8)
            ax.spines[["top", "right"]].set_visible(False)
    axes[0, 0].set_ylabel("Adults ha$^{-1}$ (log)")
    axes[1, 0].set_ylabel("Adults ha$^{-1}$ (log)")
    axes[0, 0].legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(paths["graded_cull.png"], dpi=150, facecolor="#fcfcfb")
    plt.close(fig)

    report = {"status": "layer_weight_screen", "choice": choice, "identity_max_rel_delta": max(identity.values()),
              "summary": json.loads(summary.astype({"layer": str}).to_json(orient="records")), "run_trajectories_sha256": inputs,
              "script_sha256": sha256(Path(__file__)),
              "outputs": {n: sha256(p) for n, p in paths.items() if n != "graded_cull_audit.json"}}
    paths["graded_cull_audit.json"].write_text(json.dumps(report, indent=2, default=str) + "\n")
    pd.set_option("display.width", 260)
    print(summary[["cull", "layer", "n_seeds", "loss_v3", "d_v3_vs_binary", "loss_v2_manta", "diver_error", "g2_passes_ahalf30", "eligible"] +
                  [f"{r}_diver" for r in ("lizard", "macgillivray", "north", "eyrie")]].round(2).to_string(index=False))
    print(summary[["cull", "layer"] + [f"{r}_first_ge5" for r in ("lizard", "macgillivray", "north", "eyrie")] +
                  [f"{r}_{y}" for r in ("lizard", "eyrie") for y in (2022, 2024, 2026)] + ["macgillivray_2024", "north_2024"]].round(1).to_string(index=False))
    extra_cols = [c for c in summary if c.startswith("unmet_") or c.endswith("_max_cpue_2018_26")]
    print(summary[["cull", "layer"] + extra_cols].round(3).to_string(index=False))
    print("Choice:", choice, "| identity max:", report["identity_max_rel_delta"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    main(parser.parse_args().output.resolve())
