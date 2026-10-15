"""Audit the chlorophyll larval-survival x new coral x supply-dataset experiment
(crash_mechanism_protocol.md, "Preregistered chlorophyll larval-survival switch and 2020s
experiment", 2026-10-15). Domain v0.4, 1991-2026, W at seed pair 1.

    python audit_1991_chl_survival.py <run_original> <rescore_original> <gate_original> \\
        <run_corrected> <rescore_corrected> <gate_corrected>
Outputs go to <run_corrected>.
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
from audit_1991_switching_repair import REEFS, sha256  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
RUNS = ROOT / "sandbox/calibration/runs"
TARGETS = ROOT / "sandbox/data/derived/density_targets_2020s/targets_reef_year.csv"
SMOKE = {"original": "20261015T_06_v04_smoke_owen_global_2018_2023_1991_v1",
         "corrected": "20261015T_07_v04_smoke_owen_global_2018_2023_corrected_1991_v1"}
REPLAY_ARM = "wc_b150_q20_eb_lin_obs015_s1_t40_chloff_cv2_pnd"
SMOKE_ARM = "wc_b150_q20_eb_lin_obs015_s1_cv2_pnd"
SHORT = {"Lizard Island Reef": "liz", "MacGillivray Reef": "macg", "North Direction Reef": "nd", "Eyrie Reef": "eyr"}


def label(arm: str) -> str:
    return "old coral (W)" if "_t40_" in arm else "new coral, " + arm.split("_t60_")[1].split("_")[0]


def main(runs: dict[str, tuple[Path, Path, Path]]) -> None:
    out = runs["corrected"][0]
    paths = {n: out / n for n in ("chl_survival_audit.csv", "chl_survival_audit.json", "chl_survival.png")}
    if any(p.exists() for p in paths.values()):
        raise FileExistsError("Refusing to overwrite chl-survival audit")
    targets = pd.read_csv(TARGETS)
    rows, replay, series = [], {}, {}
    for ds, (run, rescore, gate) in runs.items():
        tables = {n: pd.read_csv(run / f"{n}.csv") for n in ("trajectories", "flows", "mechanisms", "failures")}
        if len(tables["failures"]):
            raise ValueError(f"Failed evaluations in {run.name}")
        for name in ("trajectories", "flows", "mechanisms"):
            prior = pd.read_csv(RUNS / SMOKE[ds] / f"{name}.csv")
            columns = [c for c in prior if c not in ("treatment", "reef_name", "year")]
            replay[f"{ds}/{name}"] = exact_replay(tables[name][tables[name].treatment == REPLAY_ARM],
                                                  prior[prior.treatment == SMOKE_ARM], columns)
        scores = pd.read_csv(rescore / "detectability_scores.csv")
        g2 = pd.read_csv(gate / "gate_g2prime.csv").set_index(["treatment", "a_half"])
        traj, flows, mech = tables["trajectories"], tables["flows"], tables["mechanisms"]
        for arm in sorted(set(traj.treatment)):
            s = scores[(scores.treatment == arm) & (scores.a_half == 0.0)]
            primary = s[s.score_v2_primary.astype(str).str.lower() == "true"].set_index("reef_name")
            row = {"dataset": ds, "treatment": arm, "label": label(arm),
                   "loss_v3": float(primary.loss_v3.sum()), "loss_v2_manta": float(primary.loss_v2_manta.sum()),
                   "diver_error": float(primary.diver_error.sum()), "loss_v2_1992_2024": float(primary.loss_v2.sum()),
                   "g2_pass_ahalf30": bool(g2.loc[(arm, 30.0)].passes_g2prime),
                   "g2_pass_ahalf10": bool(g2.loc[(arm, 10.0)].passes_g2prime)}
            for reef, short in SHORT.items():
                t = traj[(traj.treatment == arm) & (traj.reef_name == reef)].set_index("year")
                f = flows[(flows.treatment == arm) & (flows.reef_name == reef)].set_index("year")
                m = mech[(mech.treatment == arm) & (mech.reef_name == reef)].set_index("year")
                tg = targets[targets.target_reef == reef].set_index("report_year").target_density_ha
                row |= {f"{short}_adults_2019_2026": ";".join(f"{v:.1f}" for v in t.loc[2019:2026].adults_ha),
                        f"{short}_target_2020_2026": ";".join(f"{y}:{v:.1f}" for y, v in tg.loc[2020:2026].items()),
                        f"{short}_diver_error": float(primary.loc[reef, "diver_error"]),
                        f"{short}_recruits_2019_2020": float(f.loc[2019:2020].settled_recruits.sum()),
                        f"{short}_local_production_2019": float(f.loc[2019].local_fecundity),
                        f"{short}_second_wave_max_adults_2008_2018": float(t.loc[2008:2018].adults_ha.max()),
                        f"{short}_fast_coral_2023": float(m.loc[2023].fast_coral),
                        f"{short}_fast_coral_2026": float(m.loc[2026].fast_coral),
                        f"{short}_total_site_2024": float(m.loc[2024].total_coral_site_fraction),
                        f"{short}_total_site_2026": float(m.loc[2026].total_coral_site_fraction)}
                series[(ds, arm, reef)] = t.adults_ha
            rows.append(row)
    result = pd.DataFrame(rows)
    result.to_csv(paths["chl_survival_audit.csv"], index=False)

    fig, axes = plt.subplots(4, 2, figsize=(13, 13), sharex=True)
    colors = {"old coral (W)": "#9a9893", "new coral, chloff": "#52514e", "new coral, mb15": "#86b6ef",
              "new coral, mb30": "#3987e5", "new coral, pb15": "#cba0e6", "new coral, pb30": "#8b3fc4"}
    for col, ds in enumerate(("original", "corrected")):
        for r, reef in enumerate(REEFS):
            ax = axes[r, col]
            for (d, arm, rf), s in series.items():
                if d != ds or rf != reef:
                    continue
                s = s.loc[2008:2026]
                ax.plot(s.index, s.values, color=colors[label(arm)], linewidth=2.0 if "chloff" in arm else 1.4,
                        linestyle="--" if "_t40_" in arm else "-", label=label(arm) if r == 0 else None)
            tg = targets[targets.target_reef == reef]
            ax.scatter(tg.report_year, tg.target_density_ha, color="#eb6834", s=26, zorder=5,
                       label="diver/cull target" if r == 0 else None)
            ax.set_yscale("symlog", linthresh=1.0)
            ax.grid(color="#e4e3df", linewidth=.8)
            if col == 0:
                ax.set_ylabel(f"{reef}\nadults / ha")
        axes[0, col].set_title(f"supply: {ds} outbreak layer", loc="left", fontsize=10)
        axes[0, col].legend(fontsize=7, frameon=False, ncol=2)
    fig.tight_layout()
    fig.savefig(paths["chl_survival.png"], dpi=160, facecolor="#fcfcfb")
    plt.close(fig)

    report = {"status": "bounded_diagnostic_not_promoted", "reference_replay_max_rel_delta": replay,
              "summary": json.loads(result.to_json(orient="records")),
              "inputs": {f"{ds}/{p.name}/{n}": sha256(p / n) for ds, (run, rescore, gate) in runs.items()
                         for p, n in ((run, "metadata.toml"), (run, "trajectories.csv"),
                                      (rescore, "detectability_scores.csv"), (gate, "gate_g2prime.csv"))}
                        | {"targets_reef_year.csv": sha256(TARGETS)},
              "script_sha256": sha256(Path(__file__)),
              "outputs": {n: sha256(p) for n, p in paths.items() if n != "chl_survival_audit.json"}}
    paths["chl_survival_audit.json"].write_text(json.dumps(report, indent=2, sort_keys=True, default=str) + "\n")
    pd.set_option("display.width", 320)
    pd.set_option("display.max_columns", 60)
    print(result[["dataset", "label", "loss_v3", "loss_v2_manta", "diver_error", "loss_v2_1992_2024", "g2_pass_ahalf30",
                  "liz_diver_error", "nd_diver_error", "eyr_diver_error", "macg_diver_error"]].round(2).to_string(index=False))
    for short in SHORT.values():
        print(result[["dataset", "label", f"{short}_adults_2019_2026", f"{short}_recruits_2019_2020",
                      f"{short}_second_wave_max_adults_2008_2018", f"{short}_fast_coral_2026", f"{short}_total_site_2026"]]
              .round(2).to_string(index=False))
    print("targets:", {s: result[f"{s}_target_2020_2026"].iloc[0] for s in SHORT.values()})
    print("Replay:", replay)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("run_original", "rescore_original", "gate_original", "run_corrected", "rescore_corrected", "gate_corrected"):
        parser.add_argument(name, type=Path)
    a = parser.parse_args()
    main({"original": (a.run_original.resolve(), a.rescore_original.resolve(), a.gate_original.resolve()),
          "corrected": (a.run_corrected.resolve(), a.rescore_corrected.resolve(), a.gate_corrected.resolve())})
