"""Audit the juvenile-bank retest (crash_mechanism_protocol.md, "Preregistered: juvenile bank retest on the
new coral", 2026-10-19, twenty-second round): W-next x bank {off, on (_jb_)} x wave {off, _wv025}, five seeds.
Also audits the scarcity-only bank ("Preregistered: scarcity-only juvenile bank"; arms _sb_), same criteria.

- Replay: bank-off arms must equal the same arms of 20261019T_15 (wave screen).
- Paired by seed and wave setting, the bank fails if: seed-mean v3 rises > 0.5; G2' (A_half 30) passes < 4/5;
  a 1990s peak (Lizard, MacGillivray, North Direction) moves later by > 1 y; or the first year with model
  adults >= 5/ha at Lizard and Eyrie (from 2019) does not move earlier. Diver error must fall.

    python audit_1991_bank.py <run> <rescore> <gate>
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
W3_REEFS = ("Lizard Island Reef", "MacGillivray Reef", "North Direction Reef")
V3_TOLERANCE, MIN_SEEDS, PEAK_SHIFT = 0.5, 4, 1.0


def parts(arm: str) -> tuple[bool, bool, int]:
    return ("_jb_" in arm or "_sb_" in arm), arm.endswith("_wv025"), int(re.search(r"_s(\d)_", arm).group(1))


def first_year(a: pd.Series, level: float = 5.0):
    hit = a.loc[2019:][a.loc[2019:] >= level]
    return int(hit.index.min()) if len(hit) else 2027


def main(run: Path, rescore: Path, gate: Path) -> None:
    paths = {n: run / n for n in ("bank_audit.csv", "bank_summary.csv", "bank_audit.json", "bank.png")}
    if any(p.exists() for p in paths.values()):
        raise FileExistsError("Refusing to overwrite bank audit")
    if len(pd.read_csv(run / "failures.csv")):
        raise ValueError("Failed evaluations")
    traj = pd.read_csv(run / "trajectories.csv")
    mech = pd.read_csv(run / "mechanisms.csv")
    scores = pd.read_csv(rescore / "detectability_scores.csv")
    g2 = pd.read_csv(gate / "gate_g2prime.csv").set_index(["treatment", "a_half"])
    arms = sorted(set(traj.treatment))

    replay = {}
    for name in ("trajectories", "flows", "mechanisms"):
        cur, prior = pd.read_csv(run / f"{name}.csv"), pd.read_csv(RUNS / PRIOR / f"{name}.csv")
        cols = [c for c in prior if c not in ("treatment", "reef_name", "year")]
        for arm in (a for a in arms if not parts(a)[0]):
            replay[f"{arm}/{name}"] = exact_replay(cur[cur.treatment == arm], prior[prior.treatment == arm], cols)

    rows = []
    for arm in arms:
        jb, wv, seed = parts(arm)
        s = scores[(scores.treatment == arm) & (scores.a_half == 0.0)]
        p = s[s.score_v2_primary.astype(str).str.lower() == "true"]
        row = {"treatment": arm, "bank": jb, "wave": wv, "seed": seed, "loss_v3": float(p.loss_v3.sum()),
               "loss_v2_manta": float(p.loss_v2_manta.sum()), "diver_error": float(p.diver_error.sum()),
               "g2_pass_ahalf30": bool(g2.loc[(arm, 30.0)].passes_g2prime)}
        for reef in REEFS:
            key = reef.split()[0].lower()
            a = traj[(traj.treatment == arm) & (traj.reef_name == reef)].set_index("year").adults_ha
            m = mech[(mech.treatment == arm) & (mech.reef_name == reef)].set_index("year")
            row[f"{key}_w3_peak_year"] = int(a.loc[1993:2003].idxmax())
            row[f"{key}_w3_peak_adults"] = float(a.loc[1993:2003].max())
            row[f"{key}_w4_peak_year"] = int(a.loc[2008:2016].idxmax())
            row[f"{key}_w4_peak_adults"] = float(a.loc[2008:2016].max())
            row[f"{key}_first_year_ge5"] = first_year(a)
            for y in (2020, 2022, 2024, 2026):
                row[f"{key}_{y}_adults"] = float(a.get(y, np.nan))
            yld = m.maturation_ha / m.subadults_ha.shift(1)
            row[f"{key}_yield_2015_18"] = float(yld.loc[2015:2018].mean())
            row[f"{key}_subadults_2018"] = float(m.subadults_ha.get(2018, np.nan))
            row[f"{key}_pit_2015_20"] = float(m.juvenile_predation_ha.loc[2015:2020].sum())
        rows.append(row)
    d = pd.DataFrame(rows).sort_values(["wave", "bank", "seed"])
    d.to_csv(paths["bank_audit.csv"], index=False)

    num = [c for c in d if c not in ("treatment", "bank", "wave", "seed")]
    summary = d.groupby(["wave", "bank"])[num].mean()
    summary["g2_passes_ahalf30"] = d.groupby(["wave", "bank"]).g2_pass_ahalf30.sum()
    verdict = {}
    for wv in (False, True):
        off, on = summary.loc[(wv, False)], summary.loc[(wv, True)]
        shifts = {r: on[f"{r.split()[0].lower()}_w3_peak_year"] - off[f"{r.split()[0].lower()}_w3_peak_year"] for r in W3_REEFS}
        earlier = {r: off[f"{r.split()[0].lower()}_first_year_ge5"] - on[f"{r.split()[0].lower()}_first_year_ge5"]
                   for r in ("Lizard Island Reef", "Eyrie Reef")}
        checks = {"d_v3": float(on.loss_v3 - off.loss_v3), "d_diver_error": float(on.diver_error - off.diver_error),
                  "g2_passes": int(on.g2_passes_ahalf30), "w3_peak_shift_years": shifts, "years_earlier_ge5": earlier}
        checks["pass"] = bool(checks["d_v3"] <= V3_TOLERANCE and checks["g2_passes"] >= MIN_SEEDS and
                              max(shifts.values()) <= PEAK_SHIFT and min(earlier.values()) > 0 and checks["d_diver_error"] < 0)
        verdict["wave_025" if wv else "wave_off"] = checks
    summary = summary.reset_index()
    summary.to_csv(paths["bank_summary.csv"], index=False)

    targets = pd.read_csv(TARGETS) if TARGETS.exists() else None
    fig, axes = plt.subplots(2, 4, figsize=(17, 7.5), sharex="col")
    styles = {(False, False): ("#0d366b", "-", "bank off"), (False, True): ("#eb6834", "-", "bank on"),
              (True, False): ("#0d366b", "--", "bank off, wave 0.25"), (True, True): ("#eb6834", "--", "bank on, wave 0.25")}
    for col, reef in enumerate(REEFS):
        for (wv, jb), (color, ls, label) in styles.items():
            sel = [a for a in arms if parts(a)[:2] == (jb, wv)]
            a = traj[traj.treatment.isin(sel) & (traj.reef_name == reef)].pivot(index="year", columns="treatment", values="adults_ha")
            axes[0, col].plot(a.index, a.mean(axis=1), color=color, linestyle=ls, linewidth=1.7, label=label)
            m = mech[mech.treatment.isin(sel) & (mech.reef_name == reef)].pivot(index="year", columns="treatment", values="subadults_ha")
            axes[1, col].plot(m.index, m.mean(axis=1), color=color, linestyle=ls, linewidth=1.7, label=label)
        if targets is not None:
            t = targets[(targets.target_reef == reef) & (targets.target_density_ha > 0)]
            if len(t):
                axes[0, col].scatter(t.report_year, t.target_density_ha, s=16, color="#52514e", zorder=5,
                                     label="2020s diver/cull target")
        axes[0, col].set_title(reef, loc="left", fontsize=10)
        axes[0, col].set_yscale("log")
        axes[1, col].set_yscale("log")
        for a in axes[:, col]:
            a.set_xlim(1991, 2027)
            a.grid(color="#e4e3df", linewidth=.8)
            a.spines[["top", "right"]].set_visible(False)
            a.axvspan(2016, 2017.5, color="#f3d9cc", alpha=0.5, linewidth=0)
    axes[0, 0].set_ylabel("Adults ha$^{-1}$ (log)")
    axes[1, 0].set_ylabel("Subadults ha$^{-1}$ (log)")
    axes[0, 0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(paths["bank.png"], dpi=150, facecolor="#fcfcfb")
    plt.close(fig)

    report = {"status": "mechanism_screen", "verdict": verdict, "replay_max_rel_delta": max(replay.values()),
              "summary": json.loads(summary.to_json(orient="records")),
              "inputs": {f"{p.name}/{n}": sha256(p / n) for p, n in ((run, "metadata.toml"), (run, "trajectories.csv"),
                                                                       (run, "mechanisms.csv"),
                                                                       (rescore, "detectability_scores.csv"), (gate, "gate_g2prime.csv"))},
              "script_sha256": sha256(Path(__file__)),
              "outputs": {n: sha256(p) for n, p in paths.items() if n != "bank_audit.json"}}
    paths["bank_audit.json"].write_text(json.dumps(report, indent=2, default=str) + "\n")
    pd.set_option("display.width", 250)
    show = ["wave", "bank", "loss_v3", "loss_v2_manta", "diver_error", "g2_passes_ahalf30"] + \
        [f"{r}_{c}" for r in ("lizard", "eyrie") for c in ("first_year_ge5", "yield_2015_18", "subadults_2018", "2022_adults", "2026_adults")]
    print(summary[show].round(2).to_string(index=False))
    print(summary[["wave", "bank"] + [c for c in summary if re.search(r"_w[34]_peak_(year|adults)$", c)]].round(1).to_string(index=False))
    print(json.dumps(verdict, indent=1))
    print("Replay max:", report["replay_max_rel_delta"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("run", "rescore", "gate"):
        parser.add_argument(name, type=Path)
    a = parser.parse_args()
    main(a.run.resolve(), a.rescore.resolve(), a.gate.resolve())
