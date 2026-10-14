"""Audit the coral heat-tolerance recalibration (crash_mechanism_protocol.md, "Preregistered
coral heat-tolerance recalibration", 2026-10-15): forced COTS; bleaching scalar removed; mean
tolerance shift Delta x heritability x chronic mu x post-2017 mu multiplier.

    python audit_1991_coral_thermal.py <run>
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
from audit_1991_coral_calibration import observed_photo  # noqa: E402
from audit_1991_growth_starvation import PHOTO, PHOTO_REEFS  # noqa: E402
from audit_1991_initiation import exact_replay  # noqa: E402
from audit_1991_joint_coral import objective, reef_errors  # noqa: E402
from audit_1991_switching_repair import REEFS, sha256  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
PRIOR = ROOT / "sandbox/calibration/runs/20261014T_20_coral_refit_v2_1991_v1"
REFERENCE, PRIOR_REFERENCE = "therm_b150_d000_h030_m120_p100", "refit2_s025_mu120_rep"
RECOVERY_REEFS = ("Lizard Island Reef", "MacGillivray Reef", "North Direction Reef")
J_SLACK = 0.10


def parts(name: str) -> dict:
    p = name.split("_")
    return {"beta": int(p[1][1:]) / 100, "delta": int(p[2][1:]) / 100, "h": int(p[3][1:]) / 100,
            "mu": int(p[4][1:]) / 1000, "post": int(p[5][1:]) / 100}


def main(run: Path) -> None:
    paths = {n: run / n for n in ("coral_thermal_audit.csv", "coral_thermal_audit.json", "coral_thermal.png")}
    if any(p.exists() for p in paths.values()):
        raise FileExistsError("Refusing to overwrite coral-thermal audit")
    tables = {n: pd.read_csv(run / f"{n}.csv") for n in ("trajectories", "mechanisms", "failures")}
    if len(tables["failures"]):
        raise ValueError("Failed evaluations in run")
    mechanisms = tables["mechanisms"]
    arms = sorted(set(mechanisms.treatment))
    if len(arms) != 49 or REFERENCE not in arms:
        raise ValueError("Expected 49 arms including the reference")
    replay = {}
    for name in ("trajectories", "mechanisms"):
        prior = pd.read_csv(PRIOR / f"{name}.csv")
        columns = [c for c in prior if c not in ("treatment", "reef_name", "year")]
        replay[name] = exact_replay(tables[name][tables[name].treatment == REFERENCE],
                                    prior[prior.treatment == PRIOR_REFERENCE], columns)

    photo = observed_photo()
    manta_raw = pd.read_csv(ROOT / "sandbox/data/reef_manta.csv")
    manta_raw = manta_raw[(manta_raw.data_type == "manta") & (manta_raw.domain_category == "reef") &
                          (manta_raw.variable == "HC") & (manta_raw.purpose == "MANTA") &
                          (manta_raw.project_code == "LTMP") & (manta_raw.depth == 9)]
    manta = {reef: manta_raw[(manta_raw.reef_name == REEFS[reef]) & manta_raw.report_year.between(1992, 2024)]
             .groupby("report_year")["median"].mean() for reef in REEFS}
    windows = {"all": (1992, 2024), "pre": (1992, 2016), "post": (2018, 2024)}
    errs = {k: [] for k in windows}
    rows = []
    for (name, reef), mech in mechanisms.groupby(["treatment", "reef_name"]):
        mech = mech.set_index("year").assign(slow_coral=lambda m: m.total_coral - m.fast_coral)
        for k, w in windows.items():
            errs[k].append({"treatment": name, "reef_name": reef} | reef_errors(mech, reef, photo, manta[reef], w))
        if reef in RECOVERY_REEFS:
            m = manta[reef]
            rows.append({"treatment": name, "reef_name": reef,
                         "model_change_2023_2024": float(mech.loc[2024, "total_coral_site_fraction"] - mech.loc[2023, "total_coral_site_fraction"]),
                         "obs_change_2023_2024": float(m.get(2024, np.nan) - m.get(2023, np.nan)),
                         "model_change_2015_2017": float(mech.loc[2017, "total_coral_site_fraction"] - mech.loc[2015, "total_coral_site_fraction"]),
                         "obs_change_2015_2017": float(m.get(2017, np.nan) - m.get(2015, np.nan)),
                         "model_hc_2024": float(mech.loc[2024, "total_coral_site_fraction"]),
                         "obs_hc_2024": float(m.get(2024, np.nan)),
                         "model_fast_2023": float(mech.loc[2023, "fast_coral"])})
    j = {k: objective(pd.DataFrame(v)) for k, v in errs.items()}
    post = pd.DataFrame(errs["post"])
    post = post[post.reef_name.isin(RECOVERY_REEFS)].groupby("treatment")[["rmse_total_vs_manta", "rmse_fast_vs_photo"]].mean()
    reef_rows = pd.DataFrame(rows)
    per_arm = reef_rows.groupby("treatment").agg(
        liz_macg_change_2024=("model_change_2023_2024", lambda v: float(v.iloc[:2].mean())),
        obs_liz_macg_change_2024=("obs_change_2023_2024", lambda v: float(v.iloc[:2].mean())),
        change_2015_2017=("model_change_2015_2017", "mean"), obs_change_2015_2017=("obs_change_2015_2017", "mean"),
        model_hc_2024=("model_hc_2024", "mean"), obs_hc_2024=("obs_hc_2024", "mean"))
    result = pd.DataFrame([parts(a) | {"treatment": a} for a in arms]).set_index("treatment")
    result["J_all"], result["J_pre2017"] = j["all"], j["pre"]
    result["post_rmse_manta_hc"] = post.rmse_total_vs_manta
    result["post_rmse_photo_fast"] = post.rmse_fast_vs_photo
    result["post_recovery_error"] = post.mean(axis=1)
    result = result.join(per_arm).reset_index()

    ref = result.set_index("treatment").loc[REFERENCE]
    eligible = result[(result.J_pre2017 <= ref.J_pre2017 * (1 + J_SLACK)) & (result.treatment != REFERENCE)]
    front = [r.treatment for r in result.itertuples() if not any(
        o.J_all <= r.J_all and o.post_recovery_error <= r.post_recovery_error and
        (o.J_all < r.J_all or o.post_recovery_error < r.post_recovery_error) for o in result.itertuples())]
    result["on_front"] = result.treatment.isin(front)
    result["eligible_pre2017_J"] = result.treatment.isin(eligible.treatment)
    result.sort_values(["delta", "h", "mu", "post"]).to_csv(paths["coral_thermal_audit.csv"], index=False)

    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    colors = {0.0: "#cde2fb", 1.0: "#86b6ef", 2.0: "#3987e5", 3.0: "#0d366b"}
    for r in result.itertuples():
        marker = "s" if r.post == 0.5 else "o"
        color = "#eb6834" if r.treatment == REFERENCE else colors[r.delta]
        ax.scatter(r.J_all, r.post_recovery_error, s=20 + 300 * r.mu, marker=marker, color=color,
                   edgecolors="#0b0b0b" if r.h == 0.5 else "#fcfcfb", linewidths=1.0, zorder=3)
    ax.set_xlabel("Coral fit J 1992-2024 (lower is better)")
    ax.set_ylabel("Post-2017 recovery error (RMSE 2018-24, Lizard/MacGillivray/North Direction)")
    ax.set_title("Heat-tolerance shift (colour: +0/+1/+2/+3 DHW), mu (size), post-2017 halving (square),\n"
                 "heritability 0.5 (dark edge); orange: beta 1.5 reference", loc="left", fontsize=9)
    ax.grid(color="#e4e3df", linewidth=.8)
    fig.tight_layout()
    fig.savefig(paths["coral_thermal.png"], dpi=160, facecolor="#fcfcfb")
    plt.close(fig)

    report = {"status": "bounded_diagnostic_not_promoted", "reference_replay_max_rel_delta": replay,
              "selection": {"rule": f"J 1992-2016 within {J_SLACK:.0%} of the reference; front of J_all vs post-2017 recovery error; lead chooses",
                            "eligible": eligible.treatment.tolist(), "front": front},
              "summary": json.loads(result.to_json(orient="records")),
              "inputs": {n: sha256(run / n) for n in ("metadata.toml", "trajectories.csv", "mechanisms.csv")}
                        | {"reef_photo_transect.csv": sha256(PHOTO)},
              "script_sha256": sha256(Path(__file__)),
              "outputs": {n: sha256(p) for n, p in paths.items() if n != "coral_thermal_audit.json"}}
    paths["coral_thermal_audit.json"].write_text(json.dumps(report, indent=2, sort_keys=True, default=str) + "\n")
    pd.set_option("display.width", 300)
    pd.set_option("display.max_columns", 40)
    show = result.sort_values("post_recovery_error")
    print(show[["treatment", "J_all", "J_pre2017", "post_recovery_error", "post_rmse_manta_hc", "post_rmse_photo_fast",
                "liz_macg_change_2024", "obs_liz_macg_change_2024", "change_2015_2017", "obs_change_2015_2017",
                "model_hc_2024", "obs_hc_2024", "on_front", "eligible_pre2017_J"]].round(3).to_string(index=False))
    print("Replay:", replay)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    main(parser.parse_args().run.resolve())
