"""Audit the preregistered Mellin background-mortality follow-up
(crash_mechanism_protocol.md, 2026-10-13): growth (s, f) x modelled-cyclone strength kappa
with the Mellin et al. 2019 / COTS_Model site x year background mortality, COTS forced
from observed CPUE / q. Each arm is paired with the joint run's mu = 0 arm at the same
(s, f, kappa); the full model is selected across both runs by the frozen objective J.
Run audit_1991_joint_coral.py on the joint run first.
"""

from __future__ import annotations

import argparse
import json
import sys
from itertools import product
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from audit_1991_coral_calibration import observed_photo  # noqa: E402
from audit_1991_growth_starvation import PHOTO, PHOTO_REEFS  # noqa: E402
from audit_1991_joint_coral import FAST, KAPPA, PERIODS, SLOW, arm_name, objective, pct3, reef_errors  # noqa: E402
from audit_1991_switching_repair import REEFS, sha256  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
JOINT_RUNS = [ROOT / "sandbox/calibration/runs/20261013T_02_joint_coral_1991_v1",
              ROOT / "sandbox/calibration/runs/20261013T_04_joint_coral_1991_v1_part2",
              ROOT / "sandbox/calibration/runs/20261013T_05_joint_coral_1991_v1_part3"]
JOINT_AUDIT = JOINT_RUNS[-1] / "joint_coral_audit.csv"
N_YEARS = 34


def mellin_name(s, f, k) -> str:
    return f"mellin_s{pct3(s)}_f{pct3(f)}_k{pct3(k)}"


ARMS = [mellin_name(*v) for v in product(SLOW, FAST, KAPPA)]


def main(run: Path) -> None:
    paths = {name: run / name for name in (
        "mellin_coral_audit.csv", "mellin_coral_audit.json", "mellin_coral_fit.png")}
    if any(path.exists() for path in paths.values()):
        raise FileExistsError("Refusing to overwrite Mellin audit")
    failures = pd.read_csv(run / "failures.csv")
    if len(failures):
        raise ValueError("Failed evaluations in run")
    mechanisms = pd.read_csv(run / "mechanisms.csv")
    if set(mechanisms.treatment) != set(ARMS) or len(mechanisms) != len(ARMS) * len(REEFS) * N_YEARS:
        raise ValueError("Missing, extra or incomplete treatments")
    if (mechanisms[["adults_ha", "subadults_ha", "recruits_ha"]].abs().to_numpy() > 1e-12).any():
        raise ValueError("COTS dynamics active in a forced arm")
    cover = mechanisms[["fast_coral", "total_coral", "total_coral_site_fraction"]].to_numpy()
    if not np.isfinite(cover).all() or (cover < -1e-12).any() or (cover > 1 + 1e-9).any():
        raise ValueError("Nonfinite or out-of-range coral cover")
    joint = pd.read_csv(JOINT_AUDIT)
    joint_mechanisms = pd.concat([pd.read_csv(r / "mechanisms.csv") for r in JOINT_RUNS], ignore_index=True)

    photo = observed_photo()
    manta_raw = pd.read_csv(ROOT / "sandbox/data/reef_manta.csv")
    manta_raw = manta_raw[(manta_raw.data_type == "manta") & (manta_raw.domain_category == "reef") &
                          (manta_raw.variable == "HC") & (manta_raw.purpose == "MANTA") &
                          (manta_raw.project_code == "LTMP") & (manta_raw.depth == 9)]
    rows = []
    for arm in ARMS:
        s, f, k = (float(p[1:]) / 100 for p in arm.split("_")[1:4])
        for reef, observed_name in REEFS.items():
            mech = mechanisms[(mechanisms.treatment == arm) & (mechanisms.reef_name == reef)].set_index("year")
            mech = mech.assign(slow_coral=mech.total_coral - mech.fast_coral)
            manta = manta_raw[(manta_raw.reef_name == observed_name) &
                              manta_raw.report_year.between(1992, 2024)].groupby("report_year")["median"].mean()
            row = {"treatment": arm, "reef_name": reef, "slow_multiplier": s, "fast_multiplier": f,
                   "cyclone_kappa": k, "background_mortality": "mellin"}
            row |= reef_errors(mech, reef, photo, manta, (1992, 2024))
            for label, window in PERIODS.items():
                row |= {f"{key}_{label}": value for key, value in
                        reef_errors(mech, reef, photo, manta, window).items()}
            row |= {"model_total_site_2017_2019": float(mech.total_coral_site_fraction.loc[2017:2019].mean()),
                    "model_fast_2008_2013": float(mech.fast_coral.loc[2008:2013].mean()),
                    "model_slow_2012": float(mech.slow_coral.loc[2012])}
            rows.append(row)
    result = pd.DataFrame(rows)
    result.to_csv(paths["mellin_coral_audit.csv"], index=False)

    J_mellin = objective(result)
    J_joint = objective(joint)
    paired = []
    for s, f, k in product(SLOW, FAST, KAPPA):
        m, c = mellin_name(s, f, k), arm_name(s, f, k, 0.0)
        paired.append({"slow": s, "fast": f, "kappa": k, "J_mellin": float(J_mellin[m]),
                       "J_mu0": float(J_joint[c]), "delta": float(J_mellin[m] - J_joint[c])})
    paired = pd.DataFrame(paired)
    combined = pd.concat([J_mellin.rename("J"), J_joint.rename("J")])
    best_overall = combined.idxmin()
    best_mellin = J_mellin.idxmin()
    best_joint = J_joint.idxmin()

    shown = {best_joint: ("#7c3aed", f"best joint {best_joint}", joint_mechanisms),
             best_mellin: ("#0d9488", f"best Mellin {best_mellin}", mechanisms)}
    fig, axes = plt.subplots(4, 3, figsize=(16, 13), sharex=True)
    for index, (reef, observed_name) in enumerate(REEFS.items()):
        fast_ax, slow_ax, total_ax = axes[index]
        for arm, (color, label, table) in shown.items():
            mech = table[(table.treatment == arm) & (table.reef_name == reef)].sort_values("year")
            kwargs = dict(color=color, linewidth=1.6, label=label if index == 0 else None)
            fast_ax.plot(mech.year, mech.fast_coral, **kwargs)
            slow_ax.plot(mech.year, mech.total_coral - mech.fast_coral, **kwargs)
            total_ax.plot(mech.year, mech.total_coral_site_fraction, **kwargs)
        if reef in PHOTO_REEFS:
            p = photo[photo.model_reef == reef]
            fast_ax.scatter(p.year, p.obs_fast, s=14, color="#111827", zorder=5,
                            label="LTMP observations" if index == 0 else None)
            slow_ax.scatter(p.year, p.obs_slow, s=14, color="#111827", zorder=5)
        m = manta_raw[(manta_raw.reef_name == observed_name) & manta_raw.report_year.between(1992, 2024)]
        total_ax.scatter(m.report_year, m["median"], s=14, color="#111827", zorder=5)
        fast_ax.set_ylabel(f"{reef}\nfast (habitable)")
        slow_ax.set_ylabel("slow (habitable)")
        total_ax.set_ylabel("total (site area)")
        for axis in (fast_ax, slow_ax, total_ax):
            axis.set_xlim(1991, 2024)
            axis.set_ylim(bottom=0)
            axis.grid(alpha=.15)
    axes[0, 0].set_title("Fast coral vs photo transects")
    axes[0, 1].set_title("Slow coral vs photo transects")
    axes[0, 2].set_title("Total coral (site area) vs LTMP manta")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, fontsize=9, bbox_to_anchor=(.5, -.01))
    fig.tight_layout(rect=(0, .04, 1, 1))
    fig.savefig(paths["mellin_coral_fit.png"], dpi=160, bbox_inches="tight")
    plt.close(fig)

    report = {
        "status": "coral_calibration_not_promoted",
        "best_overall_arm": best_overall, "best_overall_J": float(combined.min()),
        "best_mellin_arm": best_mellin, "best_mellin_J": float(J_mellin.min()),
        "best_joint_arm": best_joint, "best_joint_J": float(J_joint.min()),
        "paired_mellin_minus_mu0": json.loads(paired.to_json(orient="records")),
        "mellin_J": json.loads(J_mellin.to_json()),
        "inputs": {name: sha256(run / name) for name in ("metadata.toml", "mechanisms.csv")} |
                  {"joint_coral_audit.csv": sha256(JOINT_AUDIT),
                   "reef_photo_transect.csv": sha256(PHOTO)},
        "script_sha256": sha256(Path(__file__)),
        "outputs": {name: sha256(path) for name, path in paths.items() if name != "mellin_coral_audit.json"},
    }
    paths["mellin_coral_audit.json"].write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    pd.set_option("display.width", 250)
    print(paired.round(4).to_string(index=False))
    print("Best Mellin:", best_mellin, round(J_mellin.min(), 4), "| best joint:", best_joint,
          round(J_joint.min(), 4), "| best overall:", best_overall)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    main(parser.parse_args().run)
