"""Audit the preregistered coral re-fit v2 (crash_mechanism_protocol.md, 2026-10-14): COTS
forced from observed CPUE / 0.015, low-DHW bleaching scale s {1, 0.25, 0} below 4 DHW x chronic
mu {0.06..0.14} x forcing alignment {rep, srv}; frozen coral J; reference replays run 09."""

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
from audit_1991_joint_coral import PERIODS, objective, reef_errors  # noqa: E402
from audit_1991_switching_repair import REEFS, sha256  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
PRIOR = ROOT / "sandbox/calibration/runs/20261013T_09_bleaching_sensitivity_1991_v1"
SCALES, MUS, ALIGNS = (100, 25, 0), (60, 80, 100, 120, 140), ("rep", "srv")
REFERENCE, PRIOR_REFERENCE = "refit2_s100_mu030_rep", "bleach_mel3_k500_b150"
ARMS = [REFERENCE] + [f"refit2_s{s:03d}_mu{m:03d}_{a}" for s, m, a in product(SCALES, MUS, ALIGNS)]
N_YEARS = 34
FLOOR_TOLERANCE = 0.05


def main(run: Path) -> None:
    paths = {name: run / name for name in ("coral_refit_v2_audit.csv", "coral_refit_v2_audit.json", "coral_refit_v2.png")}
    if any(path.exists() for path in paths.values()):
        raise FileExistsError("Refusing to overwrite coral re-fit audit")
    tables = {name: pd.read_csv(run / f"{name}.csv") for name in ("trajectories", "flows", "mechanisms", "failures")}
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
        _, s, m, a = arm.split("_")
        for reef, observed_name in REEFS.items():
            mech = mechanisms[(mechanisms.treatment == arm) & (mechanisms.reef_name == reef)].set_index("year")
            mech = mech.assign(slow_coral=mech.total_coral - mech.fast_coral)
            manta = manta_raw[(manta_raw.reef_name == observed_name) &
                              manta_raw.report_year.between(1992, 2024)].groupby("report_year")["median"].mean()
            row = {"treatment": arm, "reef_name": reef, "low_dhw_scale": int(s[1:]) / 100, "mu": int(m[2:]) / 1000,
                   "alignment": a}
            row |= reef_errors(mech, reef, photo, manta, (1992, 2024))
            for label, window in PERIODS.items():
                row |= {f"{key}_{label}": value for key, value in reef_errors(mech, reef, photo, manta, window).items()}
            row |= {"fast_1992_1996": float(mech.loc[1992:1996].fast_coral.mean()),
                    "fast_2005_2013": float(mech.loc[2005:2013].fast_coral.mean()),
                    "total_site_2017_2019": float(mech.total_coral_site_fraction.loc[2017:2019].mean()),
                    "fast_slope_2020_2024": slope(range(2020, 2025), mech.fast_coral.loc[2020:2024])}
            if reef in PHOTO_REEFS:
                p = photo[photo.model_reef == reef].set_index("year")
                row |= {"obs_fast_1992_1996": float(p.loc[1992:1996].obs_fast.mean()),
                        "obs_fast_2005_2013": float(p.loc[2005:2013].obs_fast.mean())}
            row["obs_manta_2017_2019"] = float(manta.loc[2017:2019].mean()) if len(manta.loc[2017:2019]) else np.nan
            rows.append(row)
    result = pd.DataFrame(rows)
    result.to_csv(paths["coral_refit_v2_audit.csv"], index=False)
    summary = result.groupby("treatment").agg(
        low_dhw_scale=("low_dhw_scale", "first"), mu=("mu", "first"), alignment=("alignment", "first"),
        rmse_fast=("rmse_fast_vs_photo", "mean"), rmse_slow=("rmse_slow_vs_photo", "mean"),
        rmse_manta=("rmse_total_vs_manta", "mean"),
        fast_1992_1996=("fast_1992_1996", "mean"), obs_fast_1992_1996=("obs_fast_1992_1996", "mean"),
        fast_2005_2013=("fast_2005_2013", "mean"), obs_fast_2005_2013=("obs_fast_2005_2013", "mean"),
        total_2017_2019=("total_site_2017_2019", "mean"), obs_manta_2017_2019=("obs_manta_2017_2019", "mean"),
        fast_slope_2020_2024=("fast_slope_2020_2024", "mean")).reindex(ARMS)
    summary["J"] = objective(result).reindex(ARMS)
    for label in PERIODS:
        g = result.groupby("treatment")
        summary[f"J_{label}"] = (g[f"rmse_fast_vs_photo_{label}"].mean() + g[f"rmse_slow_vs_photo_{label}"].mean() +
                                 g[f"rmse_total_vs_manta_{label}"].mean()).reindex(ARMS)
    liz = result[result.reef_name == "Lizard Island Reef"].set_index("treatment")
    summary["liz_fast_1992_1996"] = liz.fast_1992_1996

    selection = {}
    grid = summary.drop(index=REFERENCE)
    for s, a in product(SCALES, ALIGNS):
        # Column access by key: `grid.align` would resolve to DataFrame.align (a method).
        sub = grid[(grid["low_dhw_scale"] == s / 100) & (grid["alignment"] == a)]
        if sub.empty:
            raise ValueError(f"No arms for s={s}, alignment={a}")
        best = sub.J.idxmin()
        selection[f"s{s:03d}_{a}"] = {"mu": float(sub.loc[best, "mu"]), "J": float(sub.loc[best, "J"]), "arm": best,
                                      "mu_at_grid_bound": bool(sub.loc[best, "mu"] in (min(MUS) / 1000, max(MUS) / 1000))}
    floor_choice = {}
    for a in ALIGNS:
        base = selection[f"s100_{a}"]["J"]
        floor_choice[a] = {f"s{s:03d}": {"best_J": selection[f"s{s:03d}_{a}"]["J"],
                                         "within_5pct_of_default": bool(selection[f"s{s:03d}_{a}"]["J"] <= (1 + FLOOR_TOLERANCE) * base)}
                           for s in (25, 0)}

    fig, axes = plt.subplots(3, 3, figsize=(16, 11), sharex=True)
    shown = {REFERENCE: ("#52514e", "--")}
    for key, color in zip(("s100_srv", "s025_srv", "s000_srv"), ("#2a78d6", "#0d366b", "#eb6834")):
        shown[selection[key]["arm"]] = (color, "-")
    for index, reef in enumerate(list(REEFS)[:3]):
        fast_ax, slow_ax, total_ax = axes[index]
        for arm, (color, style) in shown.items():
            mech = mechanisms[(mechanisms.treatment == arm) & (mechanisms.reef_name == reef)].sort_values("year")
            kw = dict(color=color, linestyle=style, linewidth=1.8, label=arm if index == 0 else None)
            fast_ax.plot(mech.year, mech.fast_coral, **kw)
            slow_ax.plot(mech.year, mech.total_coral - mech.fast_coral, **kw)
            total_ax.plot(mech.year, mech.total_coral_site_fraction, **kw)
        p = photo[photo.model_reef == reef]
        fast_ax.scatter(p.year, p.obs_fast, s=14, color="#0b0b0b", zorder=5)
        slow_ax.scatter(p.year, p.obs_slow, s=14, color="#0b0b0b", zorder=5)
        mm = manta_raw[(manta_raw.reef_name == REEFS[reef]) & manta_raw.report_year.between(1992, 2024)]
        total_ax.scatter(mm.report_year, mm["median"], s=14, color="#0b0b0b", zorder=5)
        fast_ax.set_ylabel(f"{reef}\nfast (habitable)")
        slow_ax.set_ylabel("slow (habitable)")
        total_ax.set_ylabel("total (site area)")
        for axis in (fast_ax, slow_ax, total_ax):
            axis.set_xlim(1991, 2024)
            axis.set_ylim(bottom=0)
            axis.grid(color="#e4e3df", linewidth=.8)
    axes[0, 0].set_title("Fast coral vs photo transects", loc="left", fontsize=10)
    axes[0, 1].set_title("Slow coral vs photo transects", loc="left", fontsize=10)
    axes[0, 2].set_title("Total (site area) vs LTMP manta", loc="left", fontsize=10)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, fontsize=8, frameon=False, bbox_to_anchor=(.5, -.01))
    fig.tight_layout(rect=(0, .04, 1, 1))
    fig.savefig(paths["coral_refit_v2.png"], dpi=160, bbox_inches="tight", facecolor="#fcfcfb")
    plt.close(fig)

    report = {"status": "coral_calibration_not_promoted", "reference_replay_max_rel_delta": replay,
              "selection_mu_by_floor_q": selection, "floor_choice": floor_choice,
              "summary": json.loads(summary.reset_index().to_json(orient="records")),
              "inputs": {name: sha256(run / name) for name in ("metadata.toml", "trajectories.csv", "flows.csv", "mechanisms.csv")}
                        | {"reef_photo_transect.csv": sha256(PHOTO)},
              "script_sha256": sha256(Path(__file__)),
              "outputs": {name: sha256(path) for name, path in paths.items() if name != "coral_refit_v2_audit.json"}}
    paths["coral_refit_v2_audit.json"].write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    pd.set_option("display.width", 260)
    pd.set_option("display.max_columns", 30)
    print(summary.round(3).to_string())
    print(json.dumps(selection, indent=1))
    print(json.dumps(floor_choice, indent=1))
    print("Reference replay max delta:", max(replay.values()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    main(parser.parse_args().run.resolve())
