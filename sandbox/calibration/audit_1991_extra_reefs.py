"""Audit the extra reported reefs (crash_mechanism_protocol.md, "Preregistered extra reefs",
2026-10-19): Martin, Linnet (2020s increase in manta) and Yonge, Carter, No Name (almost never
COTS), for W + h10 at Q 20 on the corr_w2 and corr_w3 layers, against LTMP manta by reporting year.

Model detected CPUE: q 0.015 x adults (linear) and q A^2 / (A + 30) (A_half 30 low-density
detection, as in G2').

    python audit_1991_extra_reefs.py <run_w2> <run_w3>
Outputs go to <run_w2>.
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
from audit_1991_switching_repair import sha256  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / "sandbox/calibration/runs"
OBS = ROOT / "sandbox/data/reef_cots.csv"
PRIOR = {2: "20261016T_02_pit_q_w2_1991_v1", 3: "20261016T_03_pit_q_w3_1991_v1"}
OUTBREAK_REEFS = ("Martin Reef", "Linnet Reef")
ABSENCE_REEFS = ("Yonge Reef", "Carter Reef", "No Name Reef")
Q, A_HALF, OUTBREAK_CPUE = 0.015, 30.0, 0.22
WAVES = {"w3": (1993, 2003), "w4": (2008, 2016), "w5": (2019, 2026)}


def detected(adults, a_half=0.0):
    return Q * adults if a_half == 0.0 else Q * adults ** 2 / (adults + a_half)


def main(runs: dict[int, Path]) -> None:
    out = runs[2]
    paths = {n: out / n for n in ("extra_reefs_audit.csv", "extra_reefs_audit.json", "extra_reefs.png")}
    if any(p.exists() for p in paths.values()):
        raise FileExistsError("Refusing to overwrite extra-reef audit")
    obs = pd.read_csv(OBS)
    replay, rows, series = {}, [], {}
    for w, run in runs.items():
        for name in ("trajectories", "flows", "mechanisms"):
            cur, prior = pd.read_csv(run / f"{name}.csv"), pd.read_csv(RUNS / PRIOR[w] / f"{name}.csv")
            cols = [c for c in prior if c not in ("treatment", "reef_name", "year")]
            for arm in sorted(set(cur.treatment)):
                replay[f"w{w}/{arm}/{name}"] = exact_replay(cur[cur.treatment == arm], prior[prior.treatment == arm], cols)
        traj = pd.read_csv(run / "trajectories_extra.csv")
        for reef in OUTBREAK_REEFS + ABSENCE_REEFS:
            a = traj[traj.reef_name == reef].pivot(index="year", columns="treatment", values="adults_ha")
            series[(w, reef)] = a
            o = obs[obs.reef_name == reef].set_index("year").cotsptow
            row = {"w": w, "reef": reef, "group": "outbreak" if reef in OUTBREAK_REEFS else "absence",
                   "obs_max_cpue": float(o.max()), "obs_years_ge_022": int((o >= OUTBREAK_CPUE).sum()),
                   "model_mean_adults_1992_2026": float(a.loc[1992:].mean(axis=1).mean()),
                   "model_max_adults": float(a.loc[1992:].mean(axis=1).max())}
            for label, ah in (("lin", 0.0), ("ah30", A_HALF)):
                c = detected(a.loc[1992:], ah)
                row[f"model_max_cpue_{label}"] = float(c.mean(axis=1).max())
                row[f"model_years_ge_022_{label}"] = float((c >= OUTBREAK_CPUE).sum().mean())
                obs_years = [y for y in o.index if 1992 <= y <= 2026]
                if obs_years:
                    mc = c.mean(axis=1).reindex(obs_years)
                    row[f"rmse_log_cpue_{label}"] = float(np.sqrt(np.mean((np.log(mc + 0.01) - np.log(o.loc[obs_years] + 0.01)) ** 2)))
            for wv, (y0, y1) in WAVES.items():
                m = a.loc[y0:y1].mean(axis=1)
                row[f"model_{wv}_peak_year"], row[f"model_{wv}_peak_adults"] = int(m.idxmax()), float(m.max())
                ow = o.loc[y0:y1]
                row[f"obs_{wv}_peak_year"] = int(ow.idxmax()) if len(ow) and ow.max() > 0 else None
                row[f"obs_{wv}_peak_cpue"] = float(ow.max()) if len(ow) else np.nan
            for y in (2022, 2023, 2024, 2025, 2026):
                row[f"model_{y}_adults"] = float(a.loc[y].mean()) if y in a.index else np.nan
                row[f"model_{y}_cpue_ah30"] = float(detected(a.loc[y], A_HALF).mean()) if y in a.index else np.nan
                row[f"obs_{y}_cpue"] = float(o.get(y, np.nan))
            row["obs_2027_cpue_late_2026_survey"] = float(o.get(2027, np.nan))
            rows.append(row)
    result = pd.DataFrame(rows)
    result.to_csv(paths["extra_reefs_audit.csv"], index=False)

    reefs = OUTBREAK_REEFS + ABSENCE_REEFS
    fig, axes = plt.subplots(len(reefs), 1, figsize=(11, 2.6 * len(reefs)), sharex=True)
    for ax, reef in zip(axes, reefs):
        for w, color in ((2, "#3987e5"), (3, "#0d366b")):
            c = detected(series[(w, reef)], A_HALF)
            ax.fill_between(c.index, c.min(axis=1), c.max(axis=1), color=color, alpha=0.15, linewidth=0)
            ax.plot(c.index, c.mean(axis=1), color=color, linewidth=1.8,
                    label=f"model, correction weight {w} (detected, A_half 30)")
        o = obs[obs.reef_name == reef]
        ax.scatter(o.year, o.cotsptow, color="#eb6834", s=20, zorder=5, label="LTMP manta (reporting year)")
        ax.axhline(OUTBREAK_CPUE, color="#52514e", linewidth=0.7, linestyle=":")
        ax.set_title(reef, loc="left", fontsize=10)
        ax.set_ylabel("COTS / tow")
        ax.grid(color="#e4e3df", linewidth=.8)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].legend(frameon=False, fontsize=8)
    axes[-1].set_xlabel("Year")
    fig.tight_layout()
    fig.savefig(paths["extra_reefs.png"], dpi=150, facecolor="#fcfcfb")
    plt.close(fig)

    report = {"status": "report_only_not_scored", "reference_replay_max_rel_delta": replay,
              "detection": f"q {Q}; A_half {A_HALF}: q A^2 / (A + A_half)", "waves": WAVES,
              "summary": json.loads(result.to_json(orient="records")),
              "inputs": {f"w{w}/{n}": sha256(run / n) for w, run in runs.items()
                         for n in ("metadata.toml", "trajectories_extra.csv")} | {"reef_cots.csv": sha256(OBS)},
              "script_sha256": sha256(Path(__file__)),
              "outputs": {n: sha256(p) for n, p in paths.items() if n != "extra_reefs_audit.json"}}
    paths["extra_reefs_audit.json"].write_text(json.dumps(report, indent=2, sort_keys=True, default=str) + "\n")
    pd.set_option("display.width", 250)
    print(result.set_index(["w", "reef"]).T.to_string(float_format=lambda v: f"{v:.3g}"))
    print("Replay max:", max(replay.values()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_w2", type=Path)
    parser.add_argument("run_w3", type=Path)
    a = parser.parse_args()
    main({2: a.run_w2.resolve(), 3: a.run_w3.resolve()})
