"""Audit the timing factorial (crash_mechanism_protocol.md, "Preregistered timing factorial",
2026-10-14; close-out B2 with the A2 consistent-q arms): cv2 x predator-pit Q {10, 15, 20} x
share set {pick, pnd} x observation q {0.015, 0.030}.

Each arm is read at its own q: obs015 arms from the q = 0.015 re-score / G2' evaluation,
obs030 arms from the --q=0.030 ones. Primary score: frozen score v2 (score_v2.jl).

    python audit_1991_timing.py <run> <rescore_q015> <rescore_q030> <gate_q015> <gate_q030>
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
PRIOR = ROOT / "sandbox/calibration/runs/20261014T_25_initial_share_1991_v2"
REFERENCE, PRIOR_REFERENCE = "wc_b150_q20_eb_lin_obs015_cv2_pick", "wc_b150_q20_eb_lin_cv2_pick"
TIMING_REEFS = ("Lizard Island Reef", "MacGillivray Reef", "North Direction Reef")
SHORT = {"Lizard Island Reef": "liz", "MacGillivray Reef": "macg", "North Direction Reef": "nd", "Eyrie Reef": "eyr"}
N_YEARS, TIE = 34, 0.10


def parts(arm: str) -> dict:
    p = arm.split("_")
    return {"Q": int(p[2][1:]), "q": int(p[5][3:]) / 1000, "share_set": p[-1]}


def years(text) -> list[int]:
    return [int(v) for v in str(text).split(";") if v and v != "nan"]


def main(run: Path, rescores: dict[float, Path], gates: dict[float, Path]) -> None:
    paths = {n: run / n for n in ("timing_audit.csv", "timing_audit.json", "timing.png")}
    if any(p.exists() for p in paths.values()):
        raise FileExistsError("Refusing to overwrite timing audit")
    tables = {n: pd.read_csv(run / f"{n}.csv") for n in ("trajectories", "flows", "mechanisms", "failures")}
    if len(tables["failures"]):
        raise ValueError("Failed evaluations in run")
    trajectories, flows, mechanisms = (tables[k] for k in ("trajectories", "flows", "mechanisms"))
    arms = sorted(set(trajectories.treatment))
    if len(arms) != 12 or REFERENCE not in arms or len(trajectories) != len(arms) * 4 * N_YEARS:
        raise ValueError("Missing or incomplete treatments")

    # Controls: exact replay of run 25's cv2_pick, identical external supply and 1991 total cover.
    replay = {}
    for name in ("trajectories", "flows", "mechanisms"):
        prior = pd.read_csv(PRIOR / f"{name}.csv")
        columns = [c for c in prior if c not in ("treatment", "reef_name", "year")]
        replay[name] = exact_replay(tables[name][tables[name].treatment == REFERENCE],
                                    prior[prior.treatment == PRIOR_REFERENCE], columns)
    ext = flows.pivot_table(index=["reef_name", "year"], columns="treatment", values="external_immigration")
    external = float(ext.sub(ext[REFERENCE], axis=0).abs().to_numpy().max())
    start = mechanisms[mechanisms.year == 1991].pivot(index="reef_name", columns="treatment", values="total_coral")
    total_delta = float(start.sub(start[REFERENCE], axis=0).abs().to_numpy().max())
    if external > 1e-9 or total_delta > 1e-9:
        raise ValueError(f"Control inputs differ between arms: external {external}, cover {total_delta}")

    photo = observed_photo()
    manta_raw = pd.read_csv(ROOT / "sandbox/data/reef_manta.csv")
    manta_raw = manta_raw[(manta_raw.data_type == "manta") & (manta_raw.domain_category == "reef") &
                          (manta_raw.variable == "HC") & (manta_raw.purpose == "MANTA") &
                          (manta_raw.project_code == "LTMP") & (manta_raw.depth == 9)]
    coral_rows = []
    for (arm, reef), mech in mechanisms.groupby(["treatment", "reef_name"]):
        mech = mech.set_index("year").assign(slow_coral=lambda m: m.total_coral - m.fast_coral)
        manta = manta_raw[(manta_raw.reef_name == REEFS[reef]) &
                          manta_raw.report_year.between(1992, 2024)].groupby("report_year")["median"].mean()
        coral_rows.append({"treatment": arm, "reef_name": reef} | reef_errors(mech, reef, photo, manta, (1992, 2024)))
    coral_j = objective(pd.DataFrame(coral_rows))

    scores = {q: pd.read_csv(r / "detectability_scores.csv") for q, r in rescores.items()}
    series = {q: pd.read_csv(r / "detectability_series.csv") for q, r in rescores.items()}
    g2 = {q: pd.read_csv(g / "gate_g2prime.csv").set_index(["treatment", "a_half"]) for q, g in gates.items()}
    rows = []
    for arm in arms:
        p = parts(arm)
        s = scores[p["q"]]
        s = s[(s.run == run.name) & (s.treatment == arm) & (s.a_half == 0.0)]
        primary = s[s.score_v2_primary.astype(str).str.lower() == "true"].set_index("reef_name")
        if len(primary) != 4 or len(s) != 8:
            raise ValueError(f"Missing re-score rows for {arm} at q {p['q']}")
        x = series[p["q"]]
        x = x[(x.run == run.name) & (x.treatment == arm) & (x.a_half == 0.0)]
        row = {"treatment": arm, **p,
               "summed_loss_v1": float(s[s.observation_treatment == "smoothed_3y"].loss.sum()),
               "summed_loss_v2": float(primary.loss_v2.sum()),
               "matched_v2": int(primary.matched_peaks.sum()),
               "observed_v2": int(sum(len(years(v)) for v in primary.obs_peak_years)),
               "coral_J": float(coral_j[arm])}
        errors = []
        for reef, short in SHORT.items():
            sim, obs = years(primary.loc[reef, "sim_peak_years"]), years(primary.loc[reef, "obs_peak_years"])
            reef_x = x[x.reef_name == reef].set_index("year")
            row |= {f"{short}_loss_v2": float(primary.loc[reef, "loss_v2"]),
                    f"{short}_first_peak": sim[0] if sim else np.nan,
                    f"{short}_obs_first_peak": obs[0] if obs else np.nan,
                    f"{short}_first_wave_cpue": float(reef_x.loc[1993:2001].detected_cpue.max()),
                    f"{short}_second_wave_cpue": float(reef_x.loc[2009:2017].detected_cpue.max())}
            if reef in TIMING_REEFS and sim and obs:
                errors.append(abs(sim[0] - obs[0]))
                row[f"{short}_first_peak_error"] = sim[0] - obs[0]
        row["mean_abs_first_peak_error"] = float(np.mean(errors)) if len(errors) == 3 else np.nan
        row["liz_peak_adults_1995_2000"] = float(trajectories[(trajectories.treatment == arm) &
            (trajectories.reef_name == "Lizard Island Reef") & trajectories.year.between(1995, 2000)].adults_ha.max())
        for a_half in (0.0, 10.0):
            gate = g2[p["q"]].loc[(arm, a_half)]
            tag = f"ahalf{a_half:g}"
            row |= {f"g2_c1_c3_{tag}": bool(gate.c1_two_waves and gate.c2_peaks_matched_raw_and_smoothed and
                                            gate.c3_distinct_waves),
                    f"g2_pass_{tag}": bool(gate.passes_g2prime),
                    f"liz_trough_{tag}": float(gate.liz_trough_model_mean_cpue),
                    f"macg_trough_{tag}": float(gate.macg_trough_model_mean_cpue)}
        rows.append(row)
    result = pd.DataFrame(rows).sort_values("summed_loss_v2").reset_index(drop=True)
    result.to_csv(paths["timing_audit.csv"], index=False)

    # Preregistered selection: eligible = 7/7 v2 peaks and G2' criteria 1-3 (A_half 0); lowest v2;
    # within 0.10 -> mean |first-peak error| -> Lizard + MacGillivray trough at A_half 10.
    def select(frame: pd.DataFrame):
        ok = frame[(frame.matched_v2 == frame.observed_v2) & frame.g2_c1_c3_ahalf0]
        if ok.empty:
            return None
        near = ok[ok.summed_loss_v2 <= ok.summed_loss_v2.min() + TIE].assign(
            trough=lambda f: f.liz_trough_ahalf10 + f.macg_trough_ahalf10)
        return near.sort_values(["mean_abs_first_peak_error", "trough", "summed_loss_v2"]).iloc[0].treatment
    c4 = select(result)
    c2 = c4 if c4 is not None and parts(c4)["q"] == 0.030 else select(result[result.q == 0.030])
    reference = result.set_index("treatment").loc[REFERENCE]
    improved = result[result.mean_abs_first_peak_error < reference.mean_abs_first_peak_error].treatment.tolist()

    fig, axes = plt.subplots(4, 2, figsize=(13, 12), sharex=True)
    colors = {10: "#eb6834", 15: "#3987e5", 20: "#1baf7a"}
    cots_obs = pd.read_csv(ROOT / "sandbox/data/reef_cots.csv")
    for col, q in enumerate((0.015, 0.030)):
        x = series[q][(series[q].run == run.name) & (series[q].a_half == 0.0)]
        for r, (reef, observed_name) in enumerate(REEFS.items()):
            axis = axes[r, col]
            for arm in [a for a in arms if parts(a)["q"] == q]:
                p = parts(arm)
                y = x[(x.treatment == arm) & (x.reef_name == reef) & (x.year >= 1992)]
                axis.plot(y.year, y.detected_cpue, color=colors[p["Q"]], linewidth=2.0 if arm == REFERENCE else 1.5,
                          linestyle="-" if p["share_set"] == "pick" else "--",
                          label=f"Q {p['Q']} {p['share_set']}" if r == 0 else None)
            obs = cots_obs[(cots_obs.reef_name == observed_name) & cots_obs.year.between(1992, 2024)]
            axis.scatter(obs.year, obs.cotsptow, s=16, color="#0b0b0b", zorder=5, label="LTMP" if r == 0 else None)
            axis.set_xlim(1991, 2024)
            axis.set_ylim(bottom=0)
            axis.grid(color="#e4e3df", linewidth=.8)
            if col == 0:
                axis.set_ylabel(f"{reef}\nCOTS/tow")
        axes[0, col].set_title(f"observation q {q} (seeds inverted at the same q)", loc="left", fontsize=10)
        axes[0, col].legend(fontsize=7, frameon=False, ncol=2)
    fig.suptitle("Timing factorial: cv2 x Q x share set (solid pick: North Direction q25; dashed pnd: min)",
                 x=.01, ha="left", fontsize=11)
    fig.tight_layout()
    fig.savefig(paths["timing.png"], dpi=160, bbox_inches="tight", facecolor="#fcfcfb")
    plt.close(fig)

    report = {"status": "bounded_diagnostic_not_promoted",
              "reference_replay_max_rel_delta": replay, "external_immigration_max_abs_delta": external,
              "initial_total_cover_max_abs_delta": total_delta,
              "selection": {"C4_timing_fix": c4, "C2_pratchett_q": c2,
                            "arms_improving_first_peak_timing_over_reference": improved,
                            "rule": "7/7 v2 peaks and G2' c1-c3 at A_half 0; min summed v2; within 0.10 -> "
                                    "mean |first-peak error| (Lizard, MacGillivray, North Direction) -> "
                                    "Lizard + MacGillivray trough at A_half 10"},
              "summary": json.loads(result.to_json(orient="records")),
              "inputs": {n: sha256(run / n) for n in ("metadata.toml", "trajectories.csv", "flows.csv", "mechanisms.csv")}
                        | {f"rescore_q{q}/detectability_scores.csv": sha256(r / "detectability_scores.csv") for q, r in rescores.items()}
                        | {f"gate_q{q}/gate_g2prime.csv": sha256(g / "gate_g2prime.csv") for q, g in gates.items()}
                        | {"reef_photo_transect.csv": sha256(PHOTO)},
              "script_sha256": sha256(Path(__file__)),
              "outputs": {n: sha256(p) for n, p in paths.items() if n != "timing_audit.json"}}
    paths["timing_audit.json"].write_text(json.dumps(report, indent=2, sort_keys=True, default=str) + "\n")
    pd.set_option("display.width", 320)
    pd.set_option("display.max_columns", 60)
    show = ["treatment", "summed_loss_v1", "summed_loss_v2", "matched_v2", "coral_J",
            "liz_first_peak", "macg_first_peak", "nd_first_peak", "eyr_first_peak", "mean_abs_first_peak_error",
            "liz_first_wave_cpue", "macg_first_wave_cpue", "nd_first_wave_cpue", "liz_peak_adults_1995_2000",
            "liz_loss_v2", "macg_loss_v2", "nd_loss_v2", "eyr_loss_v2",
            "g2_c1_c3_ahalf0", "g2_pass_ahalf10", "liz_trough_ahalf10", "macg_trough_ahalf10"]
    print(result[show].round(3).to_string(index=False))
    print("Replay:", replay, "| C4:", c4, "| C2:", c2, "| improving timing:", improved)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("run", "rescore_q015", "rescore_q030", "gate_q015", "gate_q030"):
        parser.add_argument(name, type=Path)
    a = parser.parse_args()
    main(a.run.resolve(), {0.015: a.rescore_q015.resolve(), 0.030: a.rescore_q030.resolve()},
         {0.015: a.gate_q015.resolve(), 0.030: a.gate_q030.resolve()})
