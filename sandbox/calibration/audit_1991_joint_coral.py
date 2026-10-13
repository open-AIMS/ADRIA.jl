"""Audit the preregistered joint coral fit (crash_mechanism_protocol.md, 2026-10-13):
growth (s, f) x modelled-cyclone strength kappa x background mortality mu, COTS grazing
forced from observed CPUE / q. Frozen objective (as in the coral calibration):

    J = mean RMSE_fast + mean RMSE_slow (LTMP photo transects, habitable fraction)
        + mean RMSE_total (LTMP manta, site-area fraction), 1992-2024.
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
from audit_1991_coral_calibration import observed_photo, rmse  # noqa: E402
from audit_1991_growth_starvation import PHOTO, PHOTO_REEFS, slope  # noqa: E402
from audit_1991_initiation import exact_replay  # noqa: E402
from audit_1991_switching_repair import REEFS, sha256  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
PRIOR = ROOT / "sandbox/calibration/runs/20261013T_01_coral_calibration_1991_v2"
SLOW, FAST, KAPPA, MU = (0.55, 0.75, 1.00), (0.8, 1.0), (1.0, 2.5, 5.0), (0.0, 0.03, 0.06)


def pct3(x: float) -> str:
    return f"{round(100 * x):03d}"


def arm_name(s, f, k, b) -> str:
    return f"forced_s{pct3(s)}_f{pct3(f)}_k{pct3(k)}_b{pct3(b)}"


ARMS = [arm_name(*v) for v in product(SLOW, FAST, KAPPA, MU)]
REFERENCE, PRIOR_REFERENCE = arm_name(1.0, 1.0, 1.0, 0.0), "forced_s100_f100"
N_YEARS = 34
PLATEAU = 0.10
PERIODS = {"1992_2012": (1992, 2012), "2013_2024": (2013, 2024)}


def reef_errors(mech: pd.DataFrame, reef: str, photo: pd.DataFrame, manta: pd.Series,
                window: tuple[int, int]) -> dict:
    out = {}
    m = manta[(manta.index >= window[0]) & (manta.index <= window[1])]
    out["rmse_total_vs_manta"] = rmse(mech.total_coral_site_fraction.loc[m.index], m) if len(m) else np.nan
    if reef in PHOTO_REEFS:
        p = photo[(photo.model_reef == reef) & photo.year.between(*window)].set_index("year")
        out["rmse_fast_vs_photo"] = rmse(mech.fast_coral.loc[p.index], p.obs_fast) if len(p) else np.nan
        out["rmse_slow_vs_photo"] = rmse(mech.slow_coral.loc[p.index], p.obs_slow) if len(p) else np.nan
    return out


def objective(frame: pd.DataFrame, reefs=None) -> pd.Series:
    data = frame if reefs is None else frame[frame.reef_name.isin(reefs)]
    grouped = data.groupby("treatment")
    return (grouped.rmse_fast_vs_photo.mean() + grouped.rmse_slow_vs_photo.mean() +
            grouped.rmse_total_vs_manta.mean())


def main(runs: list[Path]) -> None:
    # The grid may be split across run directories (an interrupted run plus a run of
    # the remaining arms); outputs go to the last directory.
    run = runs[-1]
    paths = {name: run / name for name in (
        "joint_coral_audit.csv", "joint_coral_audit.json", "joint_coral_fit.png",
        "joint_coral_objective.png")}
    if any(path.exists() for path in paths.values()):
        raise FileExistsError("Refusing to overwrite joint-coral audit")
    tables = {name: pd.concat([pd.read_csv(r / f"{name}.csv") for r in runs], ignore_index=True)
              for name in ("trajectories", "flows", "mechanisms", "failures")}
    per_run = [set(pd.read_csv(r / "mechanisms.csv").treatment) for r in runs]
    if sum(len(arms) for arms in per_run) != len(set().union(*per_run)):
        raise ValueError("An arm appears in more than one run directory")
    if len(tables["failures"]):
        raise ValueError("Failed evaluations in run")
    mechanisms = tables["mechanisms"]
    if set(mechanisms.treatment) != set(ARMS):
        raise ValueError("Missing or extra treatment")
    if any(len(t) != len(ARMS) * len(REEFS) * N_YEARS for t in tables.values() if "treatment" in t and len(t)):
        raise ValueError("Incomplete reef-year outputs")
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
        s, f, k, b = (float(p[1:]) / 100 for p in arm.split("_")[1:5])
        for reef, observed_name in REEFS.items():
            mech = mechanisms[(mechanisms.treatment == arm) &
                              (mechanisms.reef_name == reef)].set_index("year")
            mech = mech.assign(slow_coral=mech.total_coral - mech.fast_coral)
            manta = manta_raw[(manta_raw.reef_name == observed_name) &
                              manta_raw.report_year.between(1992, 2024)].groupby("report_year")["median"].mean()
            row = {"treatment": arm, "reef_name": reef, "slow_multiplier": s, "fast_multiplier": f,
                   "cyclone_kappa": k, "background_mortality": b}
            row |= reef_errors(mech, reef, photo, manta, (1992, 2024))
            for label, window in PERIODS.items():
                row |= {f"{key}_{label}": value for key, value in
                        reef_errors(mech, reef, photo, manta, window).items()}
            row |= {
                "model_total_site_2017_2019": float(mech.total_coral_site_fraction.loc[2017:2019].mean()),
                "observed_manta_2017_2019": float(manta.loc[2017:2019].mean()) if len(manta.loc[2017:2019]) else np.nan,
                "model_fast_2008_2013": float(mech.fast_coral.loc[2008:2013].mean()),
                "model_slow_2012": float(mech.slow_coral.loc[2012]),
                "model_fast_slope_2017_2024": slope(range(2017, 2025), mech.fast_coral.loc[2017:2024]),
            }
            rows.append(row)
    result = pd.DataFrame(rows)
    result.to_csv(paths["joint_coral_audit.csv"], index=False)

    summary = result.groupby("treatment").agg(
        slow_multiplier=("slow_multiplier", "first"), fast_multiplier=("fast_multiplier", "first"),
        cyclone_kappa=("cyclone_kappa", "first"), background_mortality=("background_mortality", "first"),
        mean_rmse_fast=("rmse_fast_vs_photo", "mean"), mean_rmse_slow=("rmse_slow_vs_photo", "mean"),
        mean_rmse_total_manta=("rmse_total_vs_manta", "mean"),
        model_total_2017_2019=("model_total_site_2017_2019", "mean"),
        observed_manta_2017_2019=("observed_manta_2017_2019", "mean"),
        model_fast_2008_2013=("model_fast_2008_2013", "mean"), model_slow_2012=("model_slow_2012", "mean"),
    ).reindex(ARMS)
    summary["J"] = objective(result).reindex(ARMS)
    for label in PERIODS:
        g = result.groupby("treatment")
        summary[f"J_{label}"] = (g[f"rmse_fast_vs_photo_{label}"].mean() + g[f"rmse_slow_vs_photo_{label}"].mean() +
                                 g[f"rmse_total_vs_manta_{label}"].mean()).reindex(ARMS)
    best = summary.J.min()
    selected = summary.J.idxmin()
    chosen = summary.loc[selected, ["slow_multiplier", "fast_multiplier", "cyclone_kappa", "background_mortality"]]
    grids = {"slow_multiplier": SLOW, "fast_multiplier": FAST, "cyclone_kappa": KAPPA, "background_mortality": MU}
    boundary = {name: bool(chosen[name] in (min(grid), max(grid))) for name, grid in grids.items()}
    plateau = summary[summary.J <= (1 + PLATEAU) * best].index.tolist()
    plateau_ranges = {name: [float(summary.loc[plateau, name].min()), float(summary.loc[plateau, name].max())]
                      for name in grids}
    held_out = {}
    for reef in REEFS:
        others = [r for r in REEFS if r != reef]
        pick = objective(result, others).idxmin()
        row = result[(result.treatment == pick) & (result.reef_name == reef)].iloc[0]
        ref = result[(result.treatment == REFERENCE) & (result.reef_name == reef)].iloc[0]
        error = np.nansum([row.get("rmse_fast_vs_photo", np.nan), row.get("rmse_slow_vs_photo", np.nan),
                           row.rmse_total_vs_manta])
        reference_error = np.nansum([ref.get("rmse_fast_vs_photo", np.nan), ref.get("rmse_slow_vs_photo", np.nan),
                                     ref.rmse_total_vs_manta])
        held_out[reef] = {"selected_on_other_reefs": pick, "held_out_error": float(error),
                          "reference_error": float(reference_error)}

    shown = {arm_name(1.0, 1.0, 1.0, 0.0): ("#d97706", "reference s1.00 f1.0 k1 mu0"), selected: ("#0d9488", "selected")}
    fig, axes = plt.subplots(4, 3, figsize=(16, 13), sharex=True)
    for index, (reef, observed_name) in enumerate(REEFS.items()):
        fast_ax, slow_ax, total_ax = axes[index]
        for arm, (color, label) in shown.items():
            mech = mechanisms[(mechanisms.treatment == arm) & (mechanisms.reef_name == reef)].sort_values("year")
            kwargs = dict(color=color, linewidth=2.0 if arm == selected else 1.3,
                          label=(label if arm != selected else f"selected {arm}") if index == 0 else None)
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
    fig.savefig(paths["joint_coral_fit.png"], dpi=160, bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 4, figsize=(16, 3.8))
    for axis, (name, grid) in zip(axes, grids.items()):
        best_by_level = summary.groupby(name).J.min().reindex(grid)
        axis.plot(grid, best_by_level, marker="o")
        axis.set_xlabel(name)
        axis.set_ylabel("best J at this level")
        axis.grid(alpha=.2)
    fig.suptitle("Joint coral fit: profile of the frozen objective")
    fig.tight_layout()
    fig.savefig(paths["joint_coral_objective.png"], dpi=160)
    plt.close(fig)

    report = {
        "status": "coral_calibration_not_promoted",
        "reference_replay_max_rel_delta": replay,
        "selected_arm": selected, "selected": {k: float(v) for k, v in chosen.items()},
        "selected_J": float(best), "selected_at_grid_boundary": boundary,
        "plateau_within_10pct": plateau, "plateau_ranges": plateau_ranges,
        "leave_one_reef_out": held_out,
        "treatments": json.loads(summary.to_json(orient="index")),
        "runs": [str(r.resolve().relative_to(ROOT)).replace("\\", "/") for r in runs],
        "inputs": {f"{r.name}/{name}": sha256(r / name) for r in runs for name in (
            "trajectories.csv", "flows.csv", "mechanisms.csv") + (("metadata.toml",) if (r / "metadata.toml").exists() else ())} |
                  {"reef_photo_transect.csv": sha256(PHOTO)},
        "script_sha256": sha256(Path(__file__)),
        "outputs": {name: sha256(path) for name, path in paths.items() if name != "joint_coral_audit.json"},
    }
    paths["joint_coral_audit.json"].write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    pd.set_option("display.width", 250)
    print(summary.sort_values("J").head(12).round(4).to_string())
    print("Reference:", summary.loc[REFERENCE, ["J", "J_1992_2012", "J_2013_2024"]].round(4).to_dict())
    print("Selected:", selected, "boundary:", boundary)
    print("Plateau ranges:", plateau_ranges)
    print("Leave-one-reef-out:", json.dumps(held_out, indent=1))
    print("Reference replay max delta:", max(replay.values()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("runs", type=Path, nargs="+")
    main(parser.parse_args().runs)
