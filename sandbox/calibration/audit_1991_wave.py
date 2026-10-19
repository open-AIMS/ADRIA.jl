"""Audit the wave-exposure settlement screen (crash_mechanism_protocol.md, "Preregistered: wave-exposure
settlement mechanism", 2026-10-19, twenty-second round): W-next x {off, floor 0.5, 0.25, 0.1}, five seeds.

- Replay: off arms must equal the q25 arms of 20261016T_03 (trajectories, flows, mechanisms).
- Ribbon reefs (Yonge, Carter, No Name; never >= 0.22 COTS/tow in LTMP): success if detected CPUE
  (q 0.015, A_half 30) stays < 0.22 in every year 1992-2026 in >= 4/5 seeds.
- Calibration: failure if seed-mean v3 rises by > 0.5 over off, or G2' (A_half 30) passes < 4/5 seeds.
- Martin and Linnet must keep a 1990s and a 2014 outbreak (seed-mean detected CPUE >= 0.22).

    python audit_1991_wave.py <run> <rescore> <gate>
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
PRIOR = "20261016T_03_pit_q_w3_1991_v1"
OBS = Path(__file__).resolve().parents[2] / "sandbox/data/reef_cots.csv"
RIBBONS = ("Yonge Reef", "Carter Reef", "No Name Reef")
OUTBREAK_EXTRA = ("Martin Reef", "Linnet Reef")
SCORE_REEFS = ("Lizard Island Reef", "MacGillivray Reef", "North Direction Reef", "Eyrie Reef")
Q, A_HALF, OUTBREAK_CPUE = 0.015, 30.0, 0.22
WAVES = {"w3": (1993, 2003), "w4": (2008, 2016), "w5": (2019, 2026)}
V3_TOLERANCE, MIN_SEEDS = 0.5, 4
COLORS = {1.0: "#0d366b", 0.5: "#3987e5", 0.25: "#eb6834", 0.1: "#b6b5b0"}


def floor_of(arm: str) -> float:
    m = re.search(r"_wv(\d{3})$", arm)
    return 1.0 if m is None else int(m.group(1)) / 100


def detected(adults):
    return Q * adults ** 2 / (adults + A_HALF)


def main(run: Path, rescore: Path, gate: Path) -> None:
    paths = {n: run / n for n in ("wave_audit.csv", "wave_summary.csv", "wave_audit.json", "wave.png")}
    if any(p.exists() for p in paths.values()):
        raise FileExistsError("Refusing to overwrite wave audit")
    if len(pd.read_csv(run / "failures.csv")):
        raise ValueError("Failed evaluations")
    traj = pd.read_csv(run / "trajectories.csv")
    extra = pd.read_csv(run / "trajectories_extra.csv")
    fextra = pd.read_csv(run / "flows_extra.csv")
    scores = pd.read_csv(rescore / "detectability_scores.csv")
    g2 = pd.read_csv(gate / "gate_g2prime.csv").set_index(["treatment", "a_half"])
    arms = sorted(set(traj.treatment))

    replay = {}
    for name in ("trajectories", "flows", "mechanisms"):
        cur, prior = pd.read_csv(run / f"{name}.csv"), pd.read_csv(RUNS / PRIOR / f"{name}.csv")
        cols = [c for c in prior if c not in ("treatment", "reef_name", "year")]
        for arm in (a for a in arms if floor_of(a) == 1.0):
            replay[f"{arm}/{name}"] = exact_replay(cur[cur.treatment == arm], prior[prior.treatment == arm], cols)

    rows = []
    for arm in arms:
        s = scores[(scores.treatment == arm) & (scores.a_half == 0.0)]
        p = s[s.score_v2_primary.astype(str).str.lower() == "true"]
        row = {"treatment": arm, "floor": floor_of(arm), "seed": int(re.search(r"_s(\d)_", arm).group(1)),
               "loss_v3": float(p.loss_v3.sum()), "loss_v2_manta": float(p.loss_v2_manta.sum()),
               "diver_error": float(p.diver_error.sum()), "g2_pass_ahalf30": bool(g2.loc[(arm, 30.0)].passes_g2prime)}
        for reef in RIBBONS + OUTBREAK_EXTRA:
            a = extra[(extra.treatment == arm) & (extra.reef_name == reef)].set_index("year").adults_ha.loc[1992:]
            c = detected(a)
            key = reef.split()[0].lower()
            row[f"{key}_max_adults"] = float(a.max())
            row[f"{key}_max_cpue"] = float(c.max())
            row[f"{key}_years_ge_022"] = int((c >= OUTBREAK_CPUE).sum())
            for w, (y0, y1) in WAVES.items():
                row[f"{key}_{w}_max_cpue"] = float(c.loc[y0:y1].max())
        for reef in SCORE_REEFS:
            a = traj[(traj.treatment == arm) & (traj.reef_name == reef)].set_index("year").adults_ha
            key = reef.split()[0].lower()
            for w, (y0, y1) in WAVES.items():
                row[f"{key}_{w}_peak_year"] = int(a.loc[y0:y1].idxmax())
                row[f"{key}_{w}_peak_adults"] = float(a.loc[y0:y1].max())
        f = fextra[(fextra.treatment == arm) & fextra.reef_name.isin(RIBBONS) & fextra.year.between(1993, 2016)]
        row["ribbon_settled_recruits_mean"] = float(f.settled_recruits.mean())
        row["ribbon_external_share"] = float(f.external_immigration.sum() /
                                             (f.external_immigration + f.internal_immigration + f.local_retention).sum())
        rows.append(row)
    d = pd.DataFrame(rows).sort_values(["floor", "seed"], ascending=[False, True])
    d.to_csv(paths["wave_audit.csv"], index=False)

    num = [c for c in d if c not in ("treatment", "floor", "seed")]
    summary = d.groupby("floor")[num].mean()
    summary["g2_passes_ahalf30"] = d.groupby("floor").g2_pass_ahalf30.sum()
    for reef in RIBBONS:
        key = reef.split()[0].lower()
        summary[f"{key}_seeds_clear"] = d.groupby("floor")[f"{key}_years_ge_022"].apply(lambda s: int((s == 0).sum()))
    base = float(summary.loc[1.0, "loss_v3"])
    summary["d_v3"] = summary.loss_v3 - base
    summary["ribbons_clear"] = np.all([summary[f"{r.split()[0].lower()}_seeds_clear"] >= MIN_SEEDS for r in RIBBONS], axis=0)
    summary["calibration_ok"] = (summary.d_v3 <= V3_TOLERANCE) & (summary.g2_passes_ahalf30 >= MIN_SEEDS)
    summary["outbreak_extra_kept"] = np.all([(summary[f"{r.split()[0].lower()}_{w}_max_cpue"] >= OUTBREAK_CPUE)
                                             for r in OUTBREAK_EXTRA for w in ("w3", "w4")], axis=0)
    summary = summary.reset_index().sort_values("floor", ascending=False)
    summary.to_csv(paths["wave_summary.csv"], index=False)

    obs = pd.read_csv(OBS)
    fig, axes = plt.subplots(2, 4, figsize=(17, 7.5), sharex=True)
    panels = list(RIBBONS) + ["Martin Reef", "Linnet Reef"] + ["Lizard Island Reef", "Eyrie Reef", "MacGillivray Reef"]
    for ax, reef in zip(axes.flat, panels):
        src = extra if reef in RIBBONS + OUTBREAK_EXTRA else traj
        for fl, color in COLORS.items():
            sel = [a for a in arms if floor_of(a) == fl]
            a = src[src.treatment.isin(sel) & (src.reef_name == reef)].pivot(index="year", columns="treatment", values="adults_ha")
            ax.plot(a.index, detected(a).mean(axis=1), color=color, linewidth=1.8,
                    label="off" if fl == 1.0 else f"floor {fl:g}")
        o = obs[obs.reef_name == {"Lizard Island Reef": "Lizard Isles", "MacGillivray Reef": "Macgillivray Reef"}.get(reef, reef)]
        ax.scatter(o.year, o.cotsptow, s=14, color="#52514e", zorder=5, label="LTMP manta")
        ax.axhline(OUTBREAK_CPUE, color="#b6b5b0", linewidth=0.8, linestyle="--")
        ax.set_title(reef, loc="left", fontsize=10)
        ax.set_xlim(1991, 2027)
        ax.grid(color="#e4e3df", linewidth=.8)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0, 0].set_ylabel("COTS/tow (model detected, A½ 30)")
    axes[1, 0].set_ylabel("COTS/tow (model detected, A½ 30)")
    axes[0, 0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(paths["wave.png"], dpi=150, facecolor="#fcfcfb")
    plt.close(fig)

    report = {"status": "mechanism_screen", "base_v3": base, "replay_max_rel_delta": max(replay.values()),
              "replay": replay, "summary": json.loads(summary.to_json(orient="records")),
              "inputs": {**{f"{p.name}/{n}": sha256(p / n) for p, n in ((run, "metadata.toml"), (run, "trajectories.csv"),
                                                                           (run, "trajectories_extra.csv"),
                                                                           (rescore, "detectability_scores.csv"), (gate, "gate_g2prime.csv"))},
                         "reef_cots.csv": sha256(OBS)},
              "script_sha256": sha256(Path(__file__)),
              "outputs": {n: sha256(p) for n, p in paths.items() if n != "wave_audit.json"}}
    paths["wave_audit.json"].write_text(json.dumps(report, indent=2, default=str) + "\n")
    pd.set_option("display.width", 250)
    show = ["floor", "loss_v3", "d_v3", "g2_passes_ahalf30", "yonge_max_cpue", "carter_max_cpue", "no_max_cpue",
            "yonge_seeds_clear", "carter_seeds_clear", "no_seeds_clear", "martin_w3_max_cpue", "martin_w4_max_cpue",
            "linnet_w3_max_cpue", "linnet_w4_max_cpue", "ribbon_settled_recruits_mean", "ribbon_external_share",
            "ribbons_clear", "calibration_ok", "outbreak_extra_kept"]
    print(summary[show].round(3).to_string(index=False))
    peaks = [c for c in summary if re.search(r"_(w3|w4|w5)_peak_(year|adults)$", c)]
    print(summary[["floor"] + peaks].round(1).to_string(index=False))
    print("Replay max:", report["replay_max_rel_delta"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("run", "rescore", "gate"):
        parser.add_argument(name, type=Path)
    a = parser.parse_args()
    main(a.run.resolve(), a.rescore.resolve(), a.gate.resolve())
