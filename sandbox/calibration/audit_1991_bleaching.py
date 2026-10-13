"""Audit the preregistered bleaching-response sensitivity (crash_mechanism_protocol.md,
2026-10-13): CoralBlox growth untouched (s = f = 1), background {Mellin, Mellin + 3%/yr},
cyclone kappa {2.5, 5}, bleaching scalar beta {1, 1.5, 2.5}; COTS forced from observed CPUE.
Frozen objective J as in the coral calibration; reference replays run 03 Mellin arm."""

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
from audit_1991_growth_starvation import PHOTO, PHOTO_REEFS, slope  # noqa: E402
from audit_1991_initiation import exact_replay  # noqa: E402
from audit_1991_joint_coral import PERIODS, objective, pct3, reef_errors  # noqa: E402
from audit_1991_switching_repair import REEFS, sha256  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
PRIOR = ROOT / "sandbox/calibration/runs/20261013T_03_mellin_coral_1991_v1"
LEVELS, KAPPA, BETA = ("mel", "mel3"), (2.5, 5.0), (1.0, 1.5, 2.5)


def arm_name(level, kappa, beta) -> str:
    return f"bleach_{level}_k{pct3(kappa)}_b{pct3(beta)}"


ARMS = [arm_name(*v) for v in product(LEVELS, KAPPA, BETA)]
REFERENCE, PRIOR_REFERENCE = arm_name("mel", 2.5, 1.0), "mellin_s100_f100_k250"
N_YEARS = 34


def main(run: Path) -> None:
    paths = {name: run / name for name in (
        "bleaching_audit.csv", "bleaching_audit.json", "bleaching_fit.png")}
    if any(path.exists() for path in paths.values()):
        raise FileExistsError("Refusing to overwrite bleaching audit")
    tables = {name: pd.read_csv(run / f"{name}.csv")
              for name in ("trajectories", "flows", "mechanisms", "failures")}
    if len(tables["failures"]):
        raise ValueError("Failed evaluations in run")
    mechanisms = tables["mechanisms"]
    if set(mechanisms.treatment) != set(ARMS) or len(mechanisms) != len(ARMS) * len(REEFS) * N_YEARS:
        raise ValueError("Missing, extra or incomplete treatments")
    replay = {}
    for name in ("trajectories", "flows", "mechanisms"):
        prior = pd.read_csv(PRIOR / f"{name}.csv")
        columns = [c for c in prior if c not in ("treatment", "reef_name", "year")]
        replay[name] = exact_replay(tables[name][tables[name].treatment == REFERENCE],
                                    prior[prior.treatment == PRIOR_REFERENCE], columns)
    if (mechanisms[["adults_ha", "subadults_ha", "recruits_ha"]].abs().to_numpy() > 1e-12).any():
        raise ValueError("COTS dynamics active in a forced arm")
    cover = mechanisms[["fast_coral", "total_coral", "total_coral_site_fraction"]].to_numpy()
    if not np.isfinite(cover).all() or (cover < -1e-12).any() or (cover > 1 + 1e-9).any():
        raise ValueError("Nonfinite or out-of-range coral cover")

    photo = observed_photo()
    manta_raw = pd.read_csv(ROOT / "sandbox/data/reef_manta.csv")
    manta_raw = manta_raw[(manta_raw.data_type == "manta") & (manta_raw.domain_category == "reef") &
                          (manta_raw.variable == "HC") & (manta_raw.purpose == "MANTA") &
                          (manta_raw.project_code == "LTMP") & (manta_raw.depth == 9)]
    rows = []
    for arm in ARMS:
        _, level, k, b = arm.split("_")
        for reef, observed_name in REEFS.items():
            mech = mechanisms[(mechanisms.treatment == arm) & (mechanisms.reef_name == reef)].set_index("year")
            mech = mech.assign(slow_coral=mech.total_coral - mech.fast_coral)
            manta = manta_raw[(manta_raw.reef_name == observed_name) &
                              manta_raw.report_year.between(1992, 2024)].groupby("report_year")["median"].mean()
            row = {"treatment": arm, "reef_name": reef, "background": level,
                   "cyclone_kappa": int(k[1:]) / 100, "bleaching_scalar": int(b[1:]) / 100}
            row |= reef_errors(mech, reef, photo, manta, (1992, 2024))
            for label, window in PERIODS.items():
                row |= {f"{key}_{label}": value for key, value in
                        reef_errors(mech, reef, photo, manta, window).items()}
            row |= {"model_total_site_2017_2019": float(mech.total_coral_site_fraction.loc[2017:2019].mean()),
                    "observed_manta_2017_2019": float(manta.loc[2017:2019].mean()) if len(manta.loc[2017:2019]) else np.nan,
                    "model_fast_2005_2013": float(mech.fast_coral.loc[2005:2013].mean()),
                    "model_fast_slope_2020_2024": slope(range(2020, 2025), mech.fast_coral.loc[2020:2024])}
            if reef in PHOTO_REEFS:
                p = photo[(photo.model_reef == reef)].set_index("year")
                row |= {"observed_fast_2005_2013": float(p.loc[2005:2013].obs_fast.mean()),
                        "observed_fast_slope_2020_2024": slope(p.loc[2020:2024].index, p.loc[2020:2024].obs_fast)}
            rows.append(row)
    result = pd.DataFrame(rows)
    result.to_csv(paths["bleaching_audit.csv"], index=False)
    summary = result.groupby("treatment").agg(
        background=("background", "first"), cyclone_kappa=("cyclone_kappa", "first"),
        bleaching_scalar=("bleaching_scalar", "first"),
        mean_rmse_fast=("rmse_fast_vs_photo", "mean"), mean_rmse_slow=("rmse_slow_vs_photo", "mean"),
        mean_rmse_total_manta=("rmse_total_vs_manta", "mean"),
        model_total_2017_2019=("model_total_site_2017_2019", "mean"),
        observed_manta_2017_2019=("observed_manta_2017_2019", "mean"),
        model_fast_2005_2013=("model_fast_2005_2013", "mean"),
        observed_fast_2005_2013=("observed_fast_2005_2013", "mean"),
        model_fast_slope_2020_2024=("model_fast_slope_2020_2024", "mean"),
        observed_fast_slope_2020_2024=("observed_fast_slope_2020_2024", "mean"),
    ).reindex(ARMS)
    summary["J"] = objective(result).reindex(ARMS)
    for label in PERIODS:
        g = result.groupby("treatment")
        summary[f"J_{label}"] = (g[f"rmse_fast_vs_photo_{label}"].mean() + g[f"rmse_slow_vs_photo_{label}"].mean() +
                                 g[f"rmse_total_vs_manta_{label}"].mean()).reindex(ARMS)
    selected = summary.J.idxmin()
    boundary = {"bleaching_scalar": bool(summary.loc[selected, "bleaching_scalar"] in (min(BETA), max(BETA))),
                "cyclone_kappa": bool(summary.loc[selected, "cyclone_kappa"] in (min(KAPPA), max(KAPPA)))}

    shown = {REFERENCE: ("#d97706", "Mellin, kappa 2.5, beta 1 (reference)"), selected: ("#0d9488", f"selected {selected}")}
    fig, axes = plt.subplots(4, 3, figsize=(16, 13), sharex=True)
    for index, (reef, observed_name) in enumerate(REEFS.items()):
        fast_ax, slow_ax, total_ax = axes[index]
        for arm, (color, label) in shown.items():
            mech = mechanisms[(mechanisms.treatment == arm) & (mechanisms.reef_name == reef)].sort_values("year")
            kwargs = dict(color=color, linewidth=1.8, label=label if index == 0 else None)
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
    fig.savefig(paths["bleaching_fit.png"], dpi=160, bbox_inches="tight")
    plt.close(fig)

    report = {
        "status": "coral_calibration_not_promoted",
        "reference_replay_max_rel_delta": replay,
        "selected_arm": selected, "selected_J": float(summary.J.min()), "selected_at_grid_boundary": boundary,
        "treatments": json.loads(summary.to_json(orient="index")),
        "inputs": {name: sha256(run / name) for name in (
            "metadata.toml", "trajectories.csv", "flows.csv", "mechanisms.csv")} |
                  {"reef_photo_transect.csv": sha256(PHOTO)},
        "script_sha256": sha256(Path(__file__)),
        "outputs": {name: sha256(path) for name, path in paths.items() if name != "bleaching_audit.json"},
    }
    paths["bleaching_audit.json"].write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    pd.set_option("display.width", 250)
    print(summary.drop(columns=["mean_rmse_fast", "mean_rmse_slow", "mean_rmse_total_manta"]).round(4).to_string())
    print("Selected:", selected, "boundary:", boundary)
    print("Reference replay max delta:", max(replay.values()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    main(parser.parse_args().run.resolve())
