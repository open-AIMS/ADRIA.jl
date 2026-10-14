"""Audit the observed-culling screen (crash_mechanism_protocol.md, "Preregistered
observed-culling switch", 2026-10-14; close-out C2): working configuration x {no cull,
observed cull with fuzzy site matching sigma 0 / 210 / 420 m}.

    python audit_1991_culling.py <run> <rescore_q015> <gate_q015>
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
from audit_1991_growth_starvation import PHOTO  # noqa: E402
from audit_1991_initiation import exact_replay  # noqa: E402
from audit_1991_joint_coral import objective, reef_errors  # noqa: E402
from audit_1991_switching_repair import REEFS, sha256  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
PRIOR = ROOT / "sandbox/calibration/runs/20261014T_32_timing_1991_v1"
BASE = "wc_b150_q20_eb_lin_obs015_s1"
ARMS = {c: f"{BASE}_{c}_cv2_pnd" for c in ("nocull", "cull000", "cull210", "cull420")}
REFERENCE, PRIOR_REFERENCE = ARMS["nocull"], "wc_b150_q20_eb_lin_obs015_cv2_pnd"
UNMET_LIMIT = 0.25
SHORT = {"Lizard Island Reef": "liz", "MacGillivray Reef": "macg", "North Direction Reef": "nd", "Eyrie Reef": "eyr"}


def main(run: Path, rescore: Path, gate: Path) -> None:
    paths = {n: run / n for n in ("culling_audit.csv", "culling_removals.csv", "culling_audit.json", "culling.png")}
    if any(p.exists() for p in paths.values()):
        raise FileExistsError("Refusing to overwrite culling audit")
    tables = {n: pd.read_csv(run / f"{n}.csv") for n in ("trajectories", "flows", "mechanisms", "failures", "culls")}
    if len(tables["failures"]):
        raise ValueError("Failed evaluations in run")
    trajectories, flows, mechanisms, culls = (tables[k] for k in ("trajectories", "flows", "mechanisms", "culls"))
    if set(trajectories.treatment) != set(ARMS.values()):
        raise ValueError("Unexpected treatments")
    replay = {}
    for name in ("trajectories", "flows", "mechanisms"):
        prior = pd.read_csv(PRIOR / f"{name}.csv")
        columns = [c for c in prior if c not in ("treatment", "reef_name", "year")]
        replay[name] = exact_replay(tables[name][tables[name].treatment == REFERENCE],
                                    prior[prior.treatment == PRIOR_REFERENCE], columns)
    if culls[culls.treatment == REFERENCE].filter(like="removed").to_numpy().any():
        raise ValueError("No-cull arm removed COTS")

    # Removal diagnostics: requested vs applied, by scope and period.
    req = [c for c in culls if c.startswith("requested")]
    rem = [c for c in culls if c.startswith("removed")]
    c = culls.assign(requested=culls[req].sum(axis=1), removed=culls[rem].sum(axis=1),
                     period=pd.cut(culls.year, [1990, 2011, 2016, 2020, 2024], labels=["1991-2011", "2012-16", "2017-20", "2021-24"]))
    removals = (c[c.treatment != REFERENCE].groupby(["treatment", "scope", "period"], observed=True)
                [["requested", "removed"]].sum().reset_index())
    removals["unmet_share"] = np.where(removals.requested > 0, 1 - removals.removed / removals.requested, np.nan)
    removals.to_csv(paths["culling_removals.csv"], index=False)
    totals = c[c.treatment != REFERENCE].groupby(["treatment", "scope"])[["requested", "removed"]].sum()
    totals["unmet_share"] = 1 - totals.removed / totals.requested.where(totals.requested > 0)
    flagged = totals[(totals.index.get_level_values("scope") != "domain") & (totals.unmet_share > UNMET_LIMIT)]

    photo = observed_photo()
    manta_raw = pd.read_csv(ROOT / "sandbox/data/reef_manta.csv")
    manta_raw = manta_raw[(manta_raw.data_type == "manta") & (manta_raw.domain_category == "reef") &
                          (manta_raw.variable == "HC") & (manta_raw.purpose == "MANTA") &
                          (manta_raw.project_code == "LTMP") & (manta_raw.depth == 9)]
    coral_rows = []
    for (name, reef), mech in mechanisms.groupby(["treatment", "reef_name"]):
        mech = mech.set_index("year").assign(slow_coral=lambda m: m.total_coral - m.fast_coral)
        manta = manta_raw[(manta_raw.reef_name == REEFS[reef]) &
                          manta_raw.report_year.between(1992, 2024)].groupby("report_year")["median"].mean()
        coral_rows.append({"treatment": name, "reef_name": reef} | reef_errors(mech, reef, photo, manta, (1992, 2024)))
    coral_j = objective(pd.DataFrame(coral_rows))

    scores = pd.read_csv(rescore / "detectability_scores.csv")
    series = pd.read_csv(rescore / "detectability_series.csv")
    g2 = pd.read_csv(gate / "gate_g2prime.csv").set_index(["treatment", "a_half"])
    rows = []
    for label, name in ARMS.items():
        s = scores[(scores.run == run.name) & (scores.treatment == name) & (scores.a_half == 0.0)]
        primary = s[s.score_v2_primary.astype(str).str.lower() == "true"].set_index("reef_name")
        row = {"arm": label, "treatment": name, "summed_loss_v2": float(primary.loss_v2.sum()),
               "summed_loss_v1": float(s[s.observation_treatment == "smoothed_3y"].loss.sum()),
               "matched_v2": int(primary.matched_peaks.sum()), "coral_J": float(coral_j[name]),
               "g2_pass_ahalf10": bool(g2.loc[(name, 10.0)].passes_g2prime),
               "liz_trough_ahalf10": float(g2.loc[(name, 10.0)].liz_trough_model_mean_cpue),
               "macg_trough_ahalf10": float(g2.loc[(name, 10.0)].macg_trough_model_mean_cpue)}
        for reef, short in SHORT.items():
            cpue = series[(series.run == run.name) & (series.treatment == name) & (series.a_half == 0.0) &
                          (series.reef_name == reef)].set_index("year").detected_cpue
            f = flows[(flows.treatment == name) & (flows.reef_name == reef)].set_index("year")
            m = mechanisms[(mechanisms.treatment == name) & (mechanisms.reef_name == reef)].set_index("year")
            row |= {f"{short}_loss_v2": float(primary.loc[reef, "loss_v2"]),
                    f"{short}_peaks": primary.loc[reef, "sim_peak_years"],
                    f"{short}_second_wave_cpue": float(cpue.loc[2008:2018].max()),
                    f"{short}_second_wave_year": int(cpue.loc[2008:2018].idxmax()),
                    f"{short}_late_mean_cpue_2016_2024": float(cpue.loc[2016:2024].mean()),
                    f"{short}_internal_imm_2012_2024": float(f.loc[2012:2024].internal_immigration.sum()),
                    f"{short}_total_coral_2016_2024": float(m.loc[2016:2024].total_coral.mean())}
        rows.append(row)
    result = pd.DataFrame(rows)
    result.to_csv(paths["culling_audit.csv"], index=False)

    fig, axes = plt.subplots(4, 1, figsize=(11, 12), sharex=True)
    colors = {"nocull": "#52514e", "cull000": "#86b6ef", "cull210": "#1baf7a", "cull420": "#0d366b"}
    cots_obs = pd.read_csv(ROOT / "sandbox/data/reef_cots.csv")
    for i, (reef, observed_name) in enumerate(REEFS.items()):
        axis = axes[i]
        for label, name in ARMS.items():
            x = series[(series.run == run.name) & (series.treatment == name) & (series.a_half == 0.0) &
                       (series.reef_name == reef) & (series.year >= 1992)]
            axis.plot(x.year, x.detected_cpue, color=colors[label], linewidth=1.8 if label in ("nocull", "cull210") else 1.2,
                      linestyle="--" if label == "nocull" else "-", label=label if i == 0 else None)
        obs = cots_obs[(cots_obs.reef_name == observed_name) & cots_obs.year.between(1992, 2024)]
        axis.scatter(obs.year, obs.cotsptow, s=16, color="#0b0b0b", zorder=5, label="LTMP" if i == 0 else None)
        axis.set_ylabel(f"{reef}\nCOTS/tow (q 0.015)")
        axis.set_xlim(1991, 2024)
        axis.set_ylim(bottom=0)
        axis.grid(color="#e4e3df", linewidth=.8)
    axes[0].set_title("Working configuration with and without observed culling (fuzzy site matching, sigma m)",
                      loc="left", fontsize=10)
    axes[0].legend(fontsize=8, frameon=False, ncol=5)
    fig.tight_layout()
    fig.savefig(paths["culling.png"], dpi=160, bbox_inches="tight", facecolor="#fcfcfb")
    plt.close(fig)

    report = {"status": "bounded_diagnostic_not_promoted", "reference_replay_max_rel_delta": replay,
              "removal_totals": {f"{a}|{b}": v for (a, b), v in totals.to_dict(orient="index").items()},
              "unmet_over_limit_on_calibration_reefs": [f"{a}|{b}" for a, b in flagged.index],
              "unmet_limit": UNMET_LIMIT,
              "summary": json.loads(result.to_json(orient="records")),
              "inputs": {n: sha256(run / n) for n in ("metadata.toml", "trajectories.csv", "flows.csv",
                                                       "mechanisms.csv", "culls.csv")}
                        | {"detectability_scores.csv": sha256(rescore / "detectability_scores.csv"),
                           "gate_g2prime.csv": sha256(gate / "gate_g2prime.csv"), "reef_photo_transect.csv": sha256(PHOTO)},
              "script_sha256": sha256(Path(__file__)),
              "outputs": {n: sha256(p) for n, p in paths.items() if n != "culling_audit.json"}}
    paths["culling_audit.json"].write_text(json.dumps(report, indent=2, sort_keys=True, default=str) + "\n")
    pd.set_option("display.width", 320)
    pd.set_option("display.max_columns", 60)
    print(totals.round(1).to_string())
    print(result.round(3).T.to_string())
    print("Replay:", replay, "| unmet > 25% on calibration reefs:", report["unmet_over_limit_on_calibration_reefs"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("rescore", type=Path)
    parser.add_argument("gate", type=Path)
    a = parser.parse_args()
    main(a.run.resolve(), a.rescore.resolve(), a.gate.resolve())
