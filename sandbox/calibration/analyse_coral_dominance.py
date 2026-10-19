"""Is fast-coral dominance too slow in this CoralBlox calibration? (lead's hypothesis, 2026-10-19,
twentieth round: the 2020s COTS build-up is late because fast-growing coral, the preferred food,
recovers too slowly after 2016-17.)

Compares, for Lizard, MacGillivray and North Direction (reefs with LTMP photo transects):
- observed, from LTMP 9 m photo transects (sandbox/data/reef_photo_transect.csv, refreshed from the
  2026-10-16 snapshot): hard coral = HARD CORAL / GROUP_LEVEL; fast coral = HARD CORAL / COMPOSITION
  rows with reefpage_category Acropora + Pocilloporidae; fast share = fast / hard coral;
  (version 2, 2026-10-19: version 1 used manta hard coral x the old photo share, which mixed habitats
  and inflated North Direction; it is kept as coral_dominance.* and superseded by coral_dominance_photo.*)
- model fast share = fast coral / total coral and total cover (mechanisms.csv / trajectories.csv)
  for W at Q 20 on the corr_w2 layer (run 20261016T_02, 5 seeds; mean and range).
Model cover is a fraction of habitable area; manta cover is a fraction of reef surface (scales are
compared for timing and shape, not level).

    python analyse_coral_dominance.py [--run <run_dir>] [--output <dir>]
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "sandbox" / "data"
PHOTO = DATA / "reef_photo_transect.csv"
FAST_TAXA = ("Acropora", "Pocilloporidae")
DEFAULT_RUN = ROOT / "sandbox/calibration/runs/20261016T_02_pit_q_w2_1991_v1"
ARM = "wc_b100_q20_eb_h10_obs015_s{k}_t60_e50_aprmak_cv2_pnd"
REEFS = {"Lizard Island Reef": "Lizard Isles", "MacGillivray Reef": "Macgillivray Reef",
         "North Direction Reef": "North Direction Island"}
DEFAULT_OUT = ROOT / "sandbox/calibration/runs/20261016T_02_pit_q_w2_1991_v1"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(run: Path, out: Path) -> None:
    paths = {n: out / n for n in ("coral_dominance_photo.csv", "coral_dominance_photo.json", "coral_dominance_photo.png")}
    if any(p.exists() for p in paths.values()):
        raise FileExistsError("Refusing to overwrite coral-dominance outputs")
    photo = pd.read_csv(PHOTO, low_memory=False)
    photo = photo[(photo.variable == "HARD CORAL") & (photo.project_code == "LTMP") & (photo.depth == 9)]
    mech = pd.read_csv(run / "mechanisms.csv")
    traj = pd.read_csv(run / "trajectories.csv")
    arms = [ARM.format(k=k) for k in range(1, 6)]
    rows = []
    for reef, ltmp in REEFS.items():
        m = mech[mech.treatment.isin(arms) & (mech.reef_name == reef)]
        fs = m.assign(fs=m.fast_coral / m.total_coral).pivot(index="year", columns="treatment", values="fs")
        tc = traj[traj.treatment.isin(arms) & (traj.reef_name == reef)].pivot(index="year", columns="treatment", values="coral_cover")
        fc = m.pivot(index="year", columns="treatment", values="fast_coral")
        ph = photo[photo.domain_name == ltmp]
        obs_c = ph[ph.purpose == "GROUP_LEVEL"].groupby("report_year")["median"].first()
        obs_f = ph[(ph.purpose == "COMPOSITION") & ph.reefpage_category.isin(FAST_TAXA)].groupby("report_year")["median"].sum()
        obs_s = (obs_f / obs_c).dropna()
        for y in range(1992, 2027):
            rows.append({"reef": reef, "year": y,
                         "model_fast_share": float(fs.loc[y].mean()) if y in fs.index else np.nan,
                         "model_total_cover": float(tc.loc[y].mean()) if y in tc.index else np.nan,
                         "model_fast_cover": float(fc.loc[y].mean()) if y in fc.index else np.nan,
                         "obs_fast_share": float(obs_s.get(y, np.nan)), "obs_photo_hc": float(obs_c.get(y, np.nan)),
                         "obs_fast_cover": float(obs_f.get(y, np.nan))})
    d = pd.DataFrame(rows)
    # Recovery summary after the 2016-17 bleaching: trough year and years to regain 2x trough cover.
    summary = {}
    for reef, g in d.groupby("reef"):
        g = g.set_index("year")
        out_r = {}
        for label, col in (("model_fast_cover", "model_fast_cover"), ("model_total_cover", "model_total_cover"),
                           ("obs_photo_hc", "obs_photo_hc"), ("obs_fast_cover", "obs_fast_cover")):
            s = g.loc[2015:2026, col].dropna()
            if len(s) < 2:
                continue
            trough_year = int(s.idxmin())
            after = s.loc[trough_year:]
            reach = after[after >= 2 * s.min()]
            out_r[label] = {"trough_year": trough_year, "trough": float(s.min()),
                            "first_year_at_2x_trough": int(reach.index.min()) if len(reach) else None,
                            "value_2022": float(g.loc[2022, col]) if np.isfinite(g.loc[2022, col]) else None}
        for label, col in (("model_fast_share", "model_fast_share"), ("obs_fast_share", "obs_fast_share")):
            s = g.loc[2015:2026, col].dropna()
            out_r[label] = {str(int(y)): round(float(v), 3) for y, v in s.items()}
        summary[reef] = out_r
    d.to_csv(paths["coral_dominance_photo.csv"], index=False)

    fig, axes = plt.subplots(2, 3, figsize=(15, 7.5), sharex=True)
    for col, reef in enumerate(REEFS):
        g = d[d.reef == reef].set_index("year")
        ax = axes[0, col]
        ax.plot(g.index, g.model_fast_share, color="#0d366b", linewidth=1.9, label="model fast share")
        ax.scatter(g.index, g.obs_fast_share, color="#eb6834", s=22, zorder=5, label="LTMP photo: Acropora + Pocilloporidae share")
        ax.set_title(reef, loc="left", fontsize=10)
        ax.set_ylim(0, 1)
        ax = axes[1, col]
        ax.plot(g.index, g.model_total_cover, color="#86b6ef", linewidth=1.6, label="model total cover (habitable)")
        ax.plot(g.index, g.model_fast_cover, color="#0d366b", linewidth=1.9, label="model fast cover")
        ax.scatter(g.index, g.obs_photo_hc, color="#b6b5b0", s=18, zorder=4, label="LTMP photo hard coral (9 m)")
        ax.scatter(g.index, g.obs_fast_cover, color="#eb6834", s=22, zorder=5, label="LTMP photo Acropora + Pocilloporidae")
        for a in axes[:, col]:
            a.grid(color="#e4e3df", linewidth=.8)
            a.spines[["top", "right"]].set_visible(False)
            a.axvspan(2016, 2017.5, color="#f3d9cc", alpha=0.5, linewidth=0)
    axes[0, 0].set_ylabel("Fast-coral share of hard coral")
    axes[1, 0].set_ylabel("Cover (fraction)")
    axes[0, 0].legend(frameon=False, fontsize=8)
    axes[1, 0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(paths["coral_dominance_photo.png"], dpi=150, facecolor="#fcfcfb")
    plt.close(fig)
    report = {"status": "diagnostic", "summary": summary,
              "inputs": {"reef_photo_transect.csv": sha256(PHOTO),
                         "mechanisms.csv": sha256(run / "mechanisms.csv"), "trajectories.csv": sha256(run / "trajectories.csv")},
              "script_sha256": sha256(Path(__file__)),
              "outputs": {n: sha256(p) for n, p in paths.items() if n != "coral_dominance_photo.json"}}
    paths["coral_dominance_photo.json"].write_text(json.dumps(report, indent=2) + "\n")
    pd.set_option("display.width", 250)
    print(d[d.year >= 2012].round(3).to_string(index=False))
    print(json.dumps(summary, indent=1)[:3000])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, default=DEFAULT_RUN)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    a = parser.parse_args()
    main(a.run.resolve(), a.output.resolve())
