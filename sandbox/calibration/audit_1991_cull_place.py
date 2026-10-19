"""Audit the cull-placement screen (crash_mechanism_protocol.md, "Preregistered: reef-level cull allocation",
2026-10-19, twenty-third round): W-wave x {no cull, site-matched sigma 210 (_cull210_), reef-proportional
(_cullreef_)}, five seeds.

- Replay: no-cull arms must equal the _wv025 arms of 20261019T_15.
- Removals: requested / removed / unmet by reef scope and window (culls.csv).
- Scores: v3, manta part, diver error, G2' (A_half 30); peaks; adults 2020-26 against the diver/cull targets.
Flags (reported, not tuned): reef-proportional unmet > 25% on a calibration reef with culls; G2' < 4/5;
dv3 > +1.0 against no cull.

    python audit_1991_cull_place.py <run> <rescore> <gate>
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
from audit_1991_switching_repair import sha256  # noqa: E402

RUNS = Path(__file__).resolve().parent / "runs"
PRIOR = "20261019T_15_wave_w3_q25_1991_v1"
TARGETS = Path(__file__).resolve().parents[2] / "sandbox/data/derived/density_targets_2020s/targets_reef_year.csv"
REEFS = ("Lizard Island Reef", "MacGillivray Reef", "North Direction Reef", "Eyrie Reef")
WINDOWS = {"2012_16": (2012, 2016), "2019_26": (2019, 2026)}
WAVES = {"w3": (1993, 2003), "w4": (2008, 2016)}
LABEL = {"none": "no cull", "site": "site-matched (σ 210 m)", "reef": "reef-proportional"}
COLOR = {"none": "#0d366b", "site": "#86b6ef", "reef": "#eb6834"}


def placement(arm: str) -> str:
    return "reef" if "_cullreef_" in arm else "site" if "_cull210_" in arm else "none"


def main(run: Path, rescore: Path, gate: Path) -> None:
    paths = {n: run / n for n in ("cull_place_audit.csv", "cull_place_summary.csv", "cull_place_removals.csv",
                                  "cull_place_audit.json", "cull_place.png")}
    if any(p.exists() for p in paths.values()):
        raise FileExistsError("Refusing to overwrite cull-placement audit")
    if len(pd.read_csv(run / "failures.csv")):
        raise ValueError("Failed evaluations")
    traj = pd.read_csv(run / "trajectories.csv")
    extra = pd.read_csv(run / "trajectories_extra.csv")
    culls = pd.read_csv(run / "culls.csv")
    scores = pd.read_csv(rescore / "detectability_scores.csv")
    g2 = pd.read_csv(gate / "gate_g2prime.csv").set_index(["treatment", "a_half"])
    targets = pd.read_csv(TARGETS)
    arms = sorted(set(traj.treatment))

    replay = {}
    prior_arms = {}
    for name in ("trajectories", "flows", "mechanisms"):
        cur, prior = pd.read_csv(run / f"{name}.csv"), pd.read_csv(RUNS / PRIOR / f"{name}.csv")
        cols = [c for c in prior if c not in ("treatment", "reef_name", "year")]
        for arm in (a for a in arms if placement(a) == "none"):
            replay[f"{arm}/{name}"] = exact_replay(cur[cur.treatment == arm], prior[prior.treatment == arm], cols)

    # Removals by scope and window (seed means).
    rem = []
    for arm in arms:
        if placement(arm) == "none":
            continue
        c = culls[culls.treatment == arm]
        for scope in ["domain", *REEFS]:
            for w, (y0, y1) in WINDOWS.items():
                x = c[(c.scope == scope) & c.year.between(y0, y1)]
                req = float(x.filter(like="requested").sum().sum())
                got = float(x.filter(like="removed").sum().sum())
                rem.append({"placement": placement(arm), "seed": int(re.search(r"_s(\d)_", arm).group(1)), "scope": scope,
                            "window": w, "requested": req, "removed": got, "unmet_fraction": 1 - got / req if req > 0 else np.nan})
    rem = pd.DataFrame(rem).groupby(["placement", "scope", "window"])[["requested", "removed", "unmet_fraction"]].mean().reset_index()
    rem.to_csv(paths["cull_place_removals.csv"], index=False)

    rows = []
    for arm in arms:
        s = scores[(scores.treatment == arm) & (scores.a_half == 0.0)]
        p = s[s.score_v2_primary.astype(str).str.lower() == "true"].set_index("reef_name")
        row = {"treatment": arm, "placement": placement(arm), "seed": int(re.search(r"_s(\d)_", arm).group(1)),
               "loss_v3": float(p.loss_v3.sum()), "loss_v2_manta": float(p.loss_v2_manta.sum()),
               "diver_error": float(p.diver_error.sum()), "g2_pass_ahalf30": bool(g2.loc[(arm, 30.0)].passes_g2prime)}
        for reef in REEFS:
            key = reef.split()[0].lower()
            row[f"{key}_v3"] = float(p.loss_v3.get(reef, np.nan))
            a = traj[(traj.treatment == arm) & (traj.reef_name == reef)].set_index("year").adults_ha
            for w, (y0, y1) in WAVES.items():
                row[f"{key}_{w}_peak_year"] = int(a.loc[y0:y1].idxmax())
                row[f"{key}_{w}_peak_adults"] = float(a.loc[y0:y1].max())
            for y in range(2020, 2027):
                row[f"{key}_{y}"] = float(a.get(y, np.nan))
        for reef in ("Martin Reef", "Linnet Reef"):
            a = extra[(extra.treatment == arm) & (extra.reef_name == reef)].set_index("year").adults_ha
            row[f"{reef.split()[0].lower()}_2026"] = float(a.get(2026, np.nan))
            row[f"{reef.split()[0].lower()}_w4_peak_adults"] = float(a.loc[2008:2016].max())
        rows.append(row)
    d = pd.DataFrame(rows)
    d.to_csv(paths["cull_place_audit.csv"], index=False)
    num = [c for c in d if c not in ("treatment", "placement", "seed")]
    summary = d.groupby("placement")[num].mean()
    summary["g2_passes_ahalf30"] = d.groupby("placement").g2_pass_ahalf30.sum()
    summary["d_v3"] = summary.loss_v3 - summary.loc["none", "loss_v3"]
    summary = summary.reindex(["none", "site", "reef"]).reset_index()
    summary.to_csv(paths["cull_place_summary.csv"], index=False)

    unmet_reef = rem[(rem.placement == "reef") & rem.scope.isin(REEFS) & (rem.requested > 0)]
    flags = {"reef_unmet_gt_25pct": unmet_reef[unmet_reef.unmet_fraction > 0.25][["scope", "window", "unmet_fraction"]].to_dict("records"),
             "g2_lt_4": summary[summary.g2_passes_ahalf30 < 4].placement.tolist(),
             "dv3_gt_1": summary[summary.d_v3 > 1.0].placement.tolist()}

    fig, axes = plt.subplots(1, 4, figsize=(17, 4.2), sharex=True)
    for ax, reef in zip(axes, REEFS):
        for pl in ("none", "site", "reef"):
            sel = [a for a in arms if placement(a) == pl]
            a = traj[traj.treatment.isin(sel) & (traj.reef_name == reef)].pivot(index="year", columns="treatment", values="adults_ha")
            ax.plot(a.index, a.mean(axis=1), color=COLOR[pl], linewidth=1.8, label=LABEL[pl])
        t = targets[(targets.target_reef == reef) & (targets.target_density_ha > 0)]
        ax.scatter(t.report_year, t.target_density_ha, s=16, color="#52514e", zorder=5, label="2020s diver/cull target")
        ax.set_yscale("log")
        ax.set_title(reef, loc="left", fontsize=10)
        ax.set_xlim(2005, 2027)
        ax.grid(color="#e4e3df", linewidth=.8)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Adults ha$^{-1}$ (log)")
    axes[0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(paths["cull_place.png"], dpi=150, facecolor="#fcfcfb")
    plt.close(fig)

    report = {"status": "report_only", "flags": flags, "replay_max_rel_delta": max(replay.values()),
              "summary": json.loads(summary.to_json(orient="records")), "removals": json.loads(rem.to_json(orient="records")),
              "inputs": {f"{p.name}/{n}": sha256(p / n) for p, n in ((run, "metadata.toml"), (run, "trajectories.csv"),
                                                                       (run, "culls.csv"), (rescore, "detectability_scores.csv"),
                                                                       (gate, "gate_g2prime.csv"))},
              "script_sha256": sha256(Path(__file__)),
              "outputs": {n: sha256(p) for n, p in paths.items() if n != "cull_place_audit.json"}}
    paths["cull_place_audit.json"].write_text(json.dumps(report, indent=2, default=str) + "\n")
    pd.set_option("display.width", 250)
    print(rem.round(2).to_string(index=False))
    print(summary[["placement", "loss_v3", "d_v3", "loss_v2_manta", "diver_error", "g2_passes_ahalf30"] +
                  [f"{r}_v3" for r in ("lizard", "macgillivray", "north", "eyrie")]].round(2).to_string(index=False))
    print(summary[["placement"] + [c for c in summary if re.search(r"_w[34]_peak_", c)]].round(1).to_string(index=False))
    print(summary[["placement"] + [c for c in summary if re.search(r"_(2022|2024|2026)$", c)] +
                  ["martin_w4_peak_adults", "linnet_w4_peak_adults"]].round(1).to_string(index=False))
    print(json.dumps(flags, indent=1, default=str))
    print("Replay max:", report["replay_max_rel_delta"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("run", "rescore", "gate"):
        parser.add_argument(name, type=Path)
    a = parser.parse_args()
    main(a.run.resolve(), a.rescore.resolve(), a.gate.resolve())
