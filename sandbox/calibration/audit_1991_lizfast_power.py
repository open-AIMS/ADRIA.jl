"""Audit the density-dependent manta observation test and the Lizard share re-selection
(crash_mechanism_protocol.md, 2026-10-19, twenty-ninth round, tests 1 and 3).

W-wave x Lizard 1991 fast increment {0, +0.04, +0.06} x 5 seeds, re-scored under the frozen linear
observation (a_half 0) and the power forms P1 (a_half -1, LTMP pairs) and P2 (a_half -2, all pairs).

Replay: 0 arms = the _wv025 arms of 20261019T_15; _lf04_ arms = those of 20261019T_48 (exact).
Per observation form, the increment with the lowest seed-mean summed v3 is "selected". Test 3: Lizard
total coral 1996 -> 1999 against Pratchett 0.32 -> 0.23; other calibration reefs' coral may not move by
more than 0.01 in any year.

    python audit_1991_lizfast_power.py <run> <rescore> <gate>
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
PRIORS = {0.0: "20261019T_15_wave_w3_q25_1991_v1", 0.04: "20261019T_48_lizfast_wave_1991_v1"}
DATA = Path(__file__).resolve().parents[2] / "sandbox" / "data"
REEFS = ("Lizard Island Reef", "MacGillivray Reef", "North Direction Reef", "Eyrie Reef")
FORMS = {0.0: "linear q 0.015", -1.0: "P1 power (LTMP pairs)", -2.0: "P2 power (all pairs)"}
COLORS = {0.0: "#0d366b", 0.04: "#eb6834", 0.06: "#2a7f62"}
CORAL_TOLERANCE = 0.01


def delta(arm: str) -> float:
    m = re.search(r"_lf(\d{2})_", arm)
    return 0.0 if m is None else int(m.group(1)) / 100


def main(run: Path, rescore: Path, gate: Path) -> None:
    paths = {n: run / n for n in ("lizfast_power_audit.csv", "lizfast_power_summary.csv",
                                  "lizfast_power_audit.json", "lizfast_power.png")}
    if any(p.exists() for p in paths.values()):
        raise FileExistsError("Refusing to overwrite audit")
    if len(pd.read_csv(run / "failures.csv")):
        raise ValueError("Failed evaluations")
    traj, mech = pd.read_csv(run / "trajectories.csv"), pd.read_csv(run / "mechanisms.csv")
    scores = pd.read_csv(rescore / "detectability_scores.csv")
    series = pd.read_csv(rescore / "detectability_series.csv")
    g2 = pd.read_csv(gate / "gate_g2prime.csv").set_index(["treatment", "a_half"])
    arms = sorted(set(traj.treatment))
    replay = {}
    for d0, prior in PRIORS.items():
        for name in ("trajectories", "flows", "mechanisms"):
            cur, ref = pd.read_csv(run / f"{name}.csv"), pd.read_csv(RUNS / prior / f"{name}.csv")
            cols = [c for c in ref if c not in ("treatment", "reef_name", "year")]
            for arm in (a for a in arms if delta(a) == d0):
                replay[f"{arm}/{name}"] = exact_replay(cur[cur.treatment == arm], ref[ref.treatment == arm], cols)
    rows = []
    for arm in arms:
        base = {"treatment": arm, "delta": delta(arm), "seed": int(re.search(r"_s(\d)_", arm).group(1))}
        liz = traj[(traj.treatment == arm) & (traj.reef_name == "Lizard Island Reef")].set_index("year")
        m = mech[(mech.treatment == arm) & (mech.reef_name == "Lizard Island Reef")].set_index("year")
        w3 = liz.adults_ha.loc[1993:2003]
        state = {"lizard_w3_peak_adults": float(w3.max()), "lizard_w3_peak_year": int(w3.idxmax()),
                 "lizard_w4_peak_adults": float(liz.adults_ha.loc[2008:2016].max()),
                 "lizard_coral_1996": float(liz.coral_cover.loc[1996]), "lizard_coral_1999": float(liz.coral_cover.loc[1999]),
                 "lizard_fast_1993_95": float(m.fast_coral.loc[1993:1995].mean())}
        for reef in REEFS[1:]:
            a = traj[(traj.treatment == arm) & (traj.reef_name == reef)].set_index("year")
            state[f"{reef.split()[0].lower()}_w3_peak_adults"] = float(a.adults_ha.loc[1993:2003].max())
        for a_half, form in FORMS.items():
            s = scores[(scores.treatment == arm) & (scores.a_half == a_half)]
            p = s[s.score_v2_primary.astype(str).str.lower() == "true"].set_index("reef_name")
            row = base | state | {"a_half": a_half, "obs_form": form, "loss_v3": float(p.loss_v3.sum()),
                                  "loss_v2_manta": float(p.loss_v2_manta.sum()), "diver_error": float(p.diver_error.sum()),
                                  "g2_pass": bool(g2.loc[(arm, a_half)].passes_g2prime)}
            for reef in REEFS:
                key = reef.split()[0].lower()
                row[f"{key}_v3"] = float(p.loss_v3.get(reef, np.nan))
                x = series[(series.treatment == arm) & (series.a_half == a_half) & (series.reef_name == reef)].set_index("year")
                row[f"{key}_w3_peak_cpue"] = float(x.detected_cpue.loc[1993:2003].max())
                row[f"{key}_trough_cpue_2000_2010"] = float(x.detected_cpue.loc[2000:2010].mean())
            rows.append(row)
    d = pd.DataFrame(rows)
    d.to_csv(paths["lizfast_power_audit.csv"], index=False)
    num = [c for c in d if c not in ("treatment", "delta", "seed", "a_half", "obs_form", "g2_pass")]
    summary = d.groupby(["a_half", "delta"])[num].mean()
    summary["g2_passes"] = d.groupby(["a_half", "delta"]).g2_pass.sum()
    summary = summary.reset_index()
    summary["d_v3"] = summary.loss_v3 - summary.groupby("a_half").loss_v3.transform(
        lambda v: v[summary.loc[v.index, "delta"] == 0.0].iloc[0])
    summary["selected"] = summary.loss_v3 == summary.groupby("a_half").loss_v3.transform("min")
    summary.to_csv(paths["lizfast_power_summary.csv"], index=False)

    # Test 3: other reefs' coral must not move with the Lizard increment.
    coral = traj.assign(delta=traj.treatment.map(delta), seed=traj.treatment.str.extract(r"_s(\d)_")[0])
    piv = coral[coral.reef_name != "Lizard Island Reef"].pivot_table(index=["reef_name", "seed", "year"], columns="delta",
                                                                      values="coral_cover")
    other_coral_max = {str(dl): float((piv[dl] - piv[0.0]).abs().max()) for dl in piv.columns if dl != 0.0}

    obs = pd.read_csv(DATA / "reef_cots.csv")
    o = obs[obs.reef_name == "Lizard Isles"].groupby("year").cotsptow.mean()
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.3), sharey=True)
    for ax, (a_half, form) in zip(axes, FORMS.items()):
        for dl, color in COLORS.items():
            sel = [a for a in arms if delta(a) == dl]
            x = series[series.treatment.isin(sel) & (series.a_half == a_half) & (series.reef_name == "Lizard Island Reef")]
            x = x.pivot(index="year", columns="treatment", values="detected_cpue")
            ax.plot(x.index, x.mean(axis=1), color=color, linewidth=1.9, label="W-wave" if dl == 0 else f"Lizard +{dl:.2f}")
            ax.fill_between(x.index, x.min(axis=1), x.max(axis=1), color=color, alpha=.15, linewidth=0)
        ax.scatter(o.index, o.values, s=14, color="#1f1e1d", zorder=5, label="LTMP manta")
        ax.set_title(f"Lizard manta-equivalent CPUE: {form}", loc="left", fontsize=10)
        ax.set_xlim(1991, 2027)
        ax.grid(color="#e4e3df", linewidth=.8)
        ax.spines[["top", "right"]].set_visible(False)
        ax.legend(frameon=False, fontsize=7.5)
    axes[0].set_ylabel("COTS per tow")
    fig.tight_layout()
    fig.savefig(paths["lizfast_power.png"], dpi=150, facecolor="#fcfcfb")
    plt.close(fig)

    report = {"status": "observation_sensitivity_not_promoted",
              "protocol": "crash_mechanism_protocol.md, 2026-10-19, twenty-ninth round, tests 1 and 3",
              "replay_max_rel_delta": max(replay.values()), "observation_forms": FORMS,
              "other_reef_coral_max_abs_delta": other_coral_max,
              "other_reef_coral_within_tolerance": all(v <= CORAL_TOLERANCE for v in other_coral_max.values()),
              "summary": json.loads(summary.to_json(orient="records")),
              "inputs": {f"{p.name}/{n}": sha256(p / n) for p, n in (
                  (run, "metadata.toml"), (run, "trajectories.csv"), (rescore, "detectability_scores.csv"),
                  (rescore, "metadata.toml"), (gate, "gate_g2prime.csv"))},
              "script_sha256": sha256(Path(__file__)),
              "outputs": {n: sha256(p) for n, p in paths.items() if n != "lizfast_power_audit.json"}}
    paths["lizfast_power_audit.json"].write_text(json.dumps(report, indent=2, default=str) + "\n")
    pd.set_option("display.width", 260)
    print(summary[["a_half", "delta", "loss_v3", "d_v3", "loss_v2_manta", "diver_error", "g2_passes", "selected"] +
                  [f"{r}_v3" for r in ("lizard", "macgillivray", "north", "eyrie")]].round(3).to_string(index=False))
    print(summary[["a_half", "delta"] + [f"{r}_w3_peak_cpue" for r in ("lizard", "macgillivray", "north", "eyrie")] +
                  ["lizard_trough_cpue_2000_2010", "macgillivray_trough_cpue_2000_2010"]].round(3).to_string(index=False))
    print(summary[summary.a_half == 0.0][["delta", "lizard_w3_peak_adults", "lizard_w3_peak_year", "lizard_w4_peak_adults",
                                          "lizard_coral_1996", "lizard_coral_1999", "lizard_fast_1993_95",
                                          "macgillivray_w3_peak_adults", "north_w3_peak_adults", "eyrie_w3_peak_adults"]]
          .round(3).to_string(index=False))
    print("Other-reef coral max |delta|:", other_coral_max, "Replay max:", report["replay_max_rel_delta"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("run", "rescore", "gate"):
        parser.add_argument(name, type=Path)
    a = parser.parse_args()
    main(a.run.resolve(), a.rescore.resolve(), a.gate.resolve())
