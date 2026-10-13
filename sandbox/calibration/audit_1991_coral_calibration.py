"""Audit the preregistered coral-side calibration with observed COTS forcing
(crash_mechanism_protocol.md, 2026-10-12): slow (massive) x fast linear-extension
multipliers, COTS grazing forced from observed CPUE / q, frozen coral objective

    J = mean RMSE_fast + mean RMSE_slow (LTMP photo transects, habitable fraction)
        + mean RMSE_total (LTMP manta, site-area fraction).
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
from audit_1991_growth_starvation import FAST_TAXA, PHOTO, PHOTO_REEFS, slope  # noqa: E402
from audit_1991_initiation import exact_replay  # noqa: E402
from audit_1991_switching_repair import REEFS, sha256  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
PRIOR = ROOT / "sandbox/calibration/runs/20261012T_15_food_grace_1991_v1"
REFERENCE, PRIOR_REFERENCE = "dynamic_ref", "acrpoc_g0"
SLOW = (0.40, 0.55, 0.70, 0.85, 1.00)
FAST = (1.0, 1.2, 1.4)


def arm_name(slow: float, fast: float) -> str:
    return f"forced_s{round(100 * slow):03d}_f{round(100 * fast)}"


FORCED = [arm_name(s, f) for s in SLOW for f in FAST]
ARMS = [REFERENCE] + FORCED
N_YEARS = 34
PLATEAU = 0.10


def observed_photo() -> pd.DataFrame:
    photo = pd.read_csv(PHOTO, low_memory=False)
    photo = photo[(photo.purpose == "COMPOSITION") & (photo.variable == "HARD CORAL") &
                  (photo.depth == 9) & photo.reef_name.isin(PHOTO_REEFS.values())]
    rows = []
    for (reef, year), group in photo.groupby(["reef_name", "report_year"]):
        total = group["median"].sum()
        fast = group[group.reefpage_category.isin(FAST_TAXA)]["median"].sum()
        rows.append({"model_reef": {v: k for k, v in PHOTO_REEFS.items()}[reef], "year": int(year),
                     "obs_fast": fast, "obs_slow": total - fast, "obs_total": total})
    return pd.DataFrame(rows)


def rmse(model: pd.Series, observed: pd.Series) -> float:
    return float(np.sqrt(np.mean((model.to_numpy() - observed.to_numpy()) ** 2)))


def main(run: Path) -> None:
    paths = {name: run / name for name in (
        "coral_calibration_audit.csv", "coral_calibration_audit.json",
        "coral_calibration_fit.png", "coral_calibration_objective.png")}
    if any(path.exists() for path in paths.values()):
        raise FileExistsError("Refusing to overwrite coral-calibration audit")
    tables = {name: pd.read_csv(run / f"{name}.csv")
              for name in ("trajectories", "flows", "mechanisms", "failures")}
    if len(tables["failures"]):
        raise ValueError("Failed evaluations in run")
    trajectories, flows, mechanisms = (tables[k] for k in ("trajectories", "flows", "mechanisms"))
    if set(mechanisms.treatment) != set(ARMS):
        raise ValueError("Missing or extra treatment")
    if any(len(t) != len(ARMS) * len(REEFS) * N_YEARS for t in (trajectories, flows, mechanisms)):
        raise ValueError("Incomplete reef-year outputs")
    replay = {}
    for name in ("trajectories", "flows", "mechanisms"):
        prior = pd.read_csv(PRIOR / f"{name}.csv")
        columns = [c for c in prior if c not in ("treatment", "reef_name", "year")]
        replay[name] = exact_replay(tables[name][tables[name].treatment == REFERENCE],
                                    prior[prior.treatment == PRIOR_REFERENCE], columns)
    forced = mechanisms[mechanisms.treatment.isin(FORCED)]
    if (forced[["adults_ha", "subadults_ha", "recruits_ha"]].abs().to_numpy() > 1e-12).any():
        raise ValueError("COTS dynamics active in a forced arm")
    cover = mechanisms[["fast_coral", "total_coral", "total_coral_site_fraction"]].to_numpy()
    if not np.isfinite(cover).all() or (cover < -1e-12).any() or (cover > 1 + 1e-9).any():
        raise ValueError("Nonfinite or out-of-range coral cover")

    photo = observed_photo()
    manta = pd.read_csv(ROOT / "sandbox/data/reef_manta.csv")
    manta = manta[(manta.data_type == "manta") & (manta.domain_category == "reef") &
                  (manta.variable == "HC") & (manta.purpose == "MANTA") &
                  (manta.project_code == "LTMP") & (manta.depth == 9)]
    rows = []
    for arm in ARMS:
        for reef, observed_name in REEFS.items():
            mech = mechanisms[(mechanisms.treatment == arm) &
                              (mechanisms.reef_name == reef)].set_index("year")
            mech = mech.assign(slow_coral=mech.total_coral - mech.fast_coral)
            m = manta[(manta.reef_name == observed_name) & manta.report_year.between(1992, 2024)]
            m = m.groupby("report_year")["median"].mean()
            row = {"treatment": arm, "reef_name": reef,
                   "slow_multiplier": float(arm[8:11]) / 100 if arm != REFERENCE else 1.0,
                   "fast_multiplier": float(arm.split("_f")[1]) / 100 if arm != REFERENCE else 1.0,
                   "rmse_total_vs_manta": rmse(mech.total_coral_site_fraction.loc[m.index], m),
                   "peak_slow_coral": float(mech.slow_coral.loc[1992:2024].max()),
                   "slow_coral_2012": float(mech.slow_coral.loc[2012]),
                   "mean_fast_2008_2013": float(mech.fast_coral.loc[2008:2013].mean()),
                   "model_fast_slope_2017_2024": slope(range(2017, 2025), mech.fast_coral.loc[2017:2024]),
                   "model_fast_slope_2001_2013": slope(range(2001, 2014), mech.fast_coral.loc[2001:2013])}
            if reef in PHOTO_REEFS:
                p = photo[(photo.model_reef == reef) & photo.year.between(1992, 2024)].set_index("year")
                late, early = p.loc[2017:2024], p.loc[2001:2013]
                row.update({
                    "rmse_fast_vs_photo": rmse(mech.fast_coral.loc[p.index], p.obs_fast),
                    "rmse_slow_vs_photo": rmse(mech.slow_coral.loc[p.index], p.obs_slow),
                    "obs_fast_slope_2017_2024": slope(late.index, late.obs_fast),
                    "obs_fast_slope_2001_2013": slope(early.index, early.obs_fast),
                    "obs_mean_fast_2008_2013": float(p.loc[2008:2013].obs_fast.mean()),
                    "obs_peak_slow": float(p.obs_slow.max()),
                })
            rows.append(row)
    result = pd.DataFrame(rows)
    result.to_csv(paths["coral_calibration_audit.csv"], index=False)
    summary = result.groupby("treatment").agg(
        slow_multiplier=("slow_multiplier", "first"), fast_multiplier=("fast_multiplier", "first"),
        mean_rmse_fast=("rmse_fast_vs_photo", "mean"), mean_rmse_slow=("rmse_slow_vs_photo", "mean"),
        mean_rmse_total_manta=("rmse_total_vs_manta", "mean"),
        mean_slow_2012=("slow_coral_2012", "mean"), mean_fast_2008_2013=("mean_fast_2008_2013", "mean"),
        mean_fast_slope_2017_2024=("model_fast_slope_2017_2024", "mean"),
    ).reindex(ARMS)
    summary["J"] = summary.mean_rmse_fast + summary.mean_rmse_slow + summary.mean_rmse_total_manta
    forced_summary = summary.loc[FORCED]
    best = forced_summary.J.min()
    selected = forced_summary.J.idxmin()
    s_star, f_star = forced_summary.loc[selected, ["slow_multiplier", "fast_multiplier"]]
    plateau = forced_summary[forced_summary.J <= (1 + PLATEAU) * best].index.tolist()
    boundary = {"slow": bool(s_star in (min(SLOW), max(SLOW))),
                "fast": bool(f_star in (min(FAST), max(FAST)))}
    photo_obs = result[result.reef_name.isin(PHOTO_REEFS)].groupby("reef_name").first()[[
        "obs_fast_slope_2017_2024", "obs_fast_slope_2001_2013", "obs_mean_fast_2008_2013", "obs_peak_slow"]]

    # Plots ---------------------------------------------------------------------------
    shown = {REFERENCE: ("#475569", "COTS dynamics (leading candidate)"),
             arm_name(1.0, 1.0): ("#d97706", "forced, s1.00 f1.0"),
             selected: ("#0d9488", f"forced, selected s{s_star:.2f} f{f_star:.1f}")}
    fig, axes = plt.subplots(4, 3, figsize=(16, 13), sharex=True)
    for index, (reef, observed_name) in enumerate(REEFS.items()):
        fast_ax, slow_ax, total_ax = axes[index]
        for arm, (color, label) in shown.items():
            mech = mechanisms[(mechanisms.treatment == arm) & (mechanisms.reef_name == reef)].sort_values("year")
            kwargs = dict(color=color, linewidth=2.0 if arm == selected else 1.3,
                          label=label if index == 0 else None)
            fast_ax.plot(mech.year, mech.fast_coral, **kwargs)
            slow_ax.plot(mech.year, mech.total_coral - mech.fast_coral, **kwargs)
            total_ax.plot(mech.year, mech.total_coral_site_fraction, **kwargs)
        if reef in PHOTO_REEFS:
            p = photo[photo.model_reef == reef]
            fast_ax.scatter(p.year, p.obs_fast, s=14, color="#111827", zorder=5,
                            label="LTMP photo transect" if index == 0 else None)
            slow_ax.scatter(p.year, p.obs_slow, s=14, color="#111827", zorder=5)
        m = manta[(manta.reef_name == observed_name) & manta.report_year.between(1992, 2024)]
        total_ax.scatter(m.report_year, m["median"], s=14, color="#111827", zorder=5)
        fast_ax.set_ylabel(f"{reef}\nfast (habitable)")
        slow_ax.set_ylabel("slow (habitable)")
        total_ax.set_ylabel("total (site area)")
        for axis in (fast_ax, slow_ax, total_ax):
            axis.set_xlim(1991, 2024)
            axis.set_ylim(bottom=0)
            axis.grid(alpha=.15)
    axes[0, 0].set_title("Fast coral vs photo transects (Acropora + Pocilloporidae)")
    axes[0, 1].set_title("Slow coral vs photo transects (other hard coral)")
    axes[0, 2].set_title("Total coral (site area) vs LTMP manta")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, fontsize=9, bbox_to_anchor=(.5, -.01))
    fig.tight_layout(rect=(0, .04, 1, 1))
    fig.savefig(paths["coral_calibration_fit.png"], dpi=160, bbox_inches="tight")
    plt.close(fig)

    grid = forced_summary.pivot(index="slow_multiplier", columns="fast_multiplier", values="J")
    fig, axis = plt.subplots(figsize=(6, 4.5))
    image = axis.imshow(grid.to_numpy(), origin="lower", cmap="viridis_r", aspect="auto")
    axis.set_xticks(range(len(grid.columns)), [f"{v:.1f}" for v in grid.columns])
    axis.set_yticks(range(len(grid.index)), [f"{v:.2f}" for v in grid.index])
    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            axis.text(j, i, f"{grid.iat[i, j]:.3f}", ha="center", va="center", color="white", fontsize=8)
    axis.set_xlabel("fast linear-extension multiplier f")
    axis.set_ylabel("slow (massive) multiplier s")
    axis.set_title("Frozen coral objective J (lower is better)")
    fig.colorbar(image, ax=axis)
    fig.tight_layout()
    fig.savefig(paths["coral_calibration_objective.png"], dpi=160)
    plt.close(fig)

    report = {
        "status": "coral_calibration_not_promoted",
        "reference_replay_max_rel_delta": replay,
        "objective": "J = mean RMSE_fast + mean RMSE_slow (photo transects, habitable) + mean RMSE_total (manta, site area)",
        "selected_arm": selected, "selected_slow_multiplier": float(s_star),
        "selected_fast_multiplier": float(f_star), "selected_J": float(best),
        "selected_at_grid_boundary": boundary, "plateau_within_10pct": plateau,
        "observed_targets": json.loads(photo_obs.to_json(orient="index")),
        "treatments": json.loads(summary.to_json(orient="index")),
        "inputs": {name: sha256(run / name) for name in (
            "metadata.toml", "trajectories.csv", "flows.csv", "mechanisms.csv")} |
                  {"reef_photo_transect.csv": sha256(PHOTO)},
        "script_sha256": sha256(Path(__file__)),
        "outputs": {name: sha256(path) for name, path in paths.items()
                    if name != "coral_calibration_audit.json"},
    }
    paths["coral_calibration_audit.json"].write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    pd.set_option("display.width", 250)
    print(summary.round(4).to_string())
    print("Selected:", selected, "boundary:", boundary, "plateau:", plateau)
    print(photo_obs.round(4).to_string())
    print("Reference replay max delta:", max(replay.values()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    main(parser.parse_args().run)
