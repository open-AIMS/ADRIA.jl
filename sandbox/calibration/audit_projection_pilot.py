"""Audit a 1991-2100 projection run (crash_mechanism_protocol.md, "Preregistered: projection design
to 2100", 2026-10-19): pass checks 1-5 and the preregistered metrics.

Checks: (1) 1991-2026 equal the _lf04_ arms of 20261019T_51 exactly (every table); (2) no failures;
(3) closed supply has zero external immigration after 2026, and wave supply's 2027-2100 mean external
pelagic coefficient is reported against 1991-2026; (5) post-2026 adults above 2 x the arm's own 1991-2026
maximum are flagged with reef and year. Check 4 (DHW bias correction) is read from the input-set metadata.

Metrics per arm x reef x period (2027-50, 2051-75, 2076-2100): outbreak episodes (>= 2 consecutive years
with adults >= 14.7/ha = 0.22 COTS/tow at q 0.015), years in outbreak, mean episode peak and duration, mean
total and fast coral, years with coral < 0.10 and < 0.05, external share of settlement.

    python audit_projection_pilot.py <run>
Outputs go to <run>: projection_metrics.csv, projection_summary.csv, projection_audit.json, projection.png
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from audit_1991_switching_repair import sha256  # noqa: E402

HERE = Path(__file__).resolve().parent
REF = HERE / "runs" / "20261019T_51_lizfast_power_1991_v1"
PROJ = HERE.parents[1] / "sandbox" / "data" / "derived" / "projection_inputs_v0.5"
OUTBREAK_HA = 0.22 / 0.015
PERIODS = {"2027-2050": (2027, 2050), "2051-2075": (2051, 2075), "2076-2100": (2076, 2100)}
SUPPLY = {"W": "resampled waves", "M": "long-term mean", "C": "closed", "R": "repeat 2012-26"}
COLORS = {"W": "#2a78d6", "M": "#eb6834", "C": "#1baf7a", "R": "#eda100"}
STYLES = {"ssp126": ":", "ssp245": "-", "ssp370": "-.", "ssp585": "--"}
TABLES = ("trajectories", "flows", "mechanisms", "trajectories_extra", "flows_extra")
CAL_REEFS = ("Lizard Island Reef", "MacGillivray Reef", "North Direction Reef", "Eyrie Reef")


def parts(arm: str) -> dict:
    m = re.search(r"_pj(\d{3})g(\d{2})([WMCR])_", arm)
    return {"ssp": "ssp" + m.group(1), "gcm": int(m.group(2)), "supply": m.group(3),
            "seed": int(re.search(r"_s(\d+)_", arm).group(1))}


def episodes(years: np.ndarray, adults: np.ndarray) -> list[tuple[int, int, float]]:
    out, start = [], None
    above = adults >= OUTBREAK_HA
    for i, a in enumerate(np.append(above, False)):
        if a and start is None:
            start = i
        elif not a and start is not None:
            if i - start >= 2:
                out.append((int(years[start]), i - start, float(adults[start:i].max())))
            start = None
    return out


def main(run: Path) -> None:
    paths = {n: run / n for n in ("projection_metrics.csv", "projection_summary.csv", "projection_audit.json", "projection.png")}
    if any(p.exists() for p in paths.values()):
        raise FileExistsError("Refusing to overwrite audit")
    failures = pd.read_csv(run / "failures.csv")
    data = {n: pd.read_csv(run / f"{n}.csv") for n in TABLES}
    arms = sorted(data["trajectories"].treatment.unique())
    # (1) replay of 1991-2026
    replay = {}
    for name in TABLES:
        ref = pd.read_csv(REF / f"{name}.csv")
        cols = [c for c in ref if c not in ("treatment", "reef_name", "year")]
        for arm in arms:
            base = re.sub(r"_pj\d{3}g\d{2}[WMCR]_", "_", arm)
            a = data[name][(data[name].treatment == arm) & (data[name].year <= 2026)].sort_values(["reef_name", "year"])
            b = ref[ref.treatment == base].sort_values(["reef_name", "year"])
            if len(a) != len(b) or len(b) == 0:
                raise ValueError(f"Replay rows differ for {arm}/{name}")
            replay[f"{arm}/{name}"] = float(np.nanmax(np.abs(a[cols].to_numpy(float) - b[cols].to_numpy(float))))
    traj = pd.concat([data["trajectories"], data["trajectories_extra"]], ignore_index=True)
    flows = pd.concat([data["flows"], data["flows_extra"]], ignore_index=True)
    mech = data["mechanisms"]
    rows, flags, supply_check = [], [], {}
    for arm in arms:
        p = parts(arm)
        f_arm = flows[flows.treatment == arm]
        post = f_arm[f_arm.year > 2026]
        if p["supply"] == "C":
            supply_check[arm] = {"post2026_external_immigration_max": float(post.external_immigration.abs().max())}
        if p["supply"] == "W":
            pre = f_arm[f_arm.year <= 2026]
            supply_check[arm] = {"external_pelagic_mean_ratio_post_pre": float(post.external_pelagic.mean() / pre.external_pelagic.mean())}
        for reef in traj[traj.treatment == arm].reef_name.unique():
            x = traj[(traj.treatment == arm) & (traj.reef_name == reef)].set_index("year").sort_index()
            hist_max = float(x.adults_ha.loc[1991:2026].max())
            post_x = x.loc[2027:]
            if post_x.adults_ha.max() > 2 * hist_max:
                flags.append({"arm": arm, "reef": reef, "year": int(post_x.adults_ha.idxmax()),
                              "adults_ha": float(post_x.adults_ha.max()), "hist_max_1991_2026": hist_max})
            fl = f_arm[f_arm.reef_name == reef].set_index("year")
            m = mech[(mech.treatment == arm) & (mech.reef_name == reef)].set_index("year")
            for period, (y0, y1) in PERIODS.items():
                w = x.loc[y0:y1]
                eps = episodes(w.index.to_numpy(), w.adults_ha.to_numpy())
                ff = fl.loc[y0:y1]
                total = (ff.external_immigration + ff.internal_immigration + ff.local_retention + ff.background_immigration).sum()
                rows.append(p | {"treatment": arm, "reef_name": reef, "period": period,
                                 "outbreak_episodes": len(eps), "years_in_outbreak": int((w.adults_ha >= OUTBREAK_HA).sum()),
                                 "mean_episode_peak_ha": float(np.mean([e[2] for e in eps])) if eps else 0.0,
                                 "mean_episode_years": float(np.mean([e[1] for e in eps])) if eps else 0.0,
                                 "max_adults_ha": float(w.adults_ha.max()),
                                 "coral_mean": float(w.coral_cover.mean()),
                                 "years_coral_below_010": int((w.coral_cover < 0.10).sum()),
                                 "years_coral_below_005": int((w.coral_cover < 0.05).sum()),
                                 "fast_coral_mean": float(m.fast_coral.loc[y0:y1].mean()) if len(m) else np.nan,
                                 "external_share": float(ff.external_immigration.sum() / total) if total > 0 else np.nan})
    metrics = pd.DataFrame(rows)
    metrics.to_csv(paths["projection_metrics.csv"], index=False)
    num = [c for c in metrics if c not in ("ssp", "gcm", "supply", "seed", "treatment", "reef_name", "period")]
    summary = metrics.groupby(["ssp", "gcm", "supply", "reef_name", "period"])[num].mean().reset_index()
    summary.to_csv(paths["projection_summary.csv"], index=False)

    meta = json.loads((PROJ / "metadata.json").read_text())
    gcms = sorted(metrics.gcm.unique())
    names = pd.read_csv(PROJ / "dhw" / "dhw_scenarios.csv")
    dhw = [r for r in meta["dhw"]["summary"] if r["ssp"] in metrics.ssp.unique() and
           r["gcm"] in set(names[names.scenario.isin(gcms)].gcm)]

    years = np.arange(1991, 2101)
    fig, axes = plt.subplots(2, 4, figsize=(17, 7), sharex=True)
    for j, reef in enumerate(CAL_REEFS):
        for (ssp, sup), g in traj[traj.treatment.isin(arms) & (traj.reef_name == reef)].assign(
                ssp=lambda d: d.treatment.map(lambda a: parts(a)["ssp"]),
                sup=lambda d: d.treatment.map(lambda a: parts(a)["supply"])).groupby(["ssp", "sup"]):
            for i, var in enumerate(("adults_ha", "coral_cover")):
                piv = g.pivot(index="year", columns="treatment", values=var).reindex(years)
                ax = axes[i, j]
                ax.plot(piv.index, piv.mean(axis=1), color=COLORS[sup], linestyle=STYLES.get(ssp, "-"), linewidth=1.6,
                        label=f"{SUPPLY[sup]}, {ssp.upper().replace('SSP', 'SSP')}")
                ax.fill_between(piv.index, piv.min(axis=1), piv.max(axis=1), color=COLORS[sup], alpha=.10, linewidth=0)
        axes[0, j].set_title(reef, loc="left", fontsize=10)
        axes[0, j].axhline(OUTBREAK_HA, color="#b6b5b0", linewidth=.8)
    axes[0, 0].set_ylabel("Adult COTS ha$^{-1}$")
    axes[1, 0].set_ylabel("Coral cover (fraction of habitable area)")
    for ax in axes.flat:
        ax.axvline(2026.5, color="#52514e", linewidth=.8, linestyle=":")
        ax.grid(color="#e4e3df", linewidth=.8)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_xlim(1991, 2100)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=6, frameon=False, fontsize=8)
    fig.suptitle("Projection pilot: seed mean, ribbon = seed range; dotted line = start of projection; grey = outbreak level",
                 x=0.01, ha="left", fontsize=10)
    fig.tight_layout(rect=(0, 0.06, 1, 0.97))
    fig.savefig(paths["projection.png"], dpi=150, facecolor="#fcfcfb")
    plt.close(fig)

    report = {"status": "projection_pilot_not_promoted",
              "protocol": "crash_mechanism_protocol.md, Preregistered: projection design to 2100 (2026-10-19)",
              "check1_replay_max_abs": max(replay.values()), "check2_failures": int(len(failures)),
              "check3_supply": supply_check, "check4_dhw": dhw, "check5_flags": flags,
              "outbreak_threshold_ha": OUTBREAK_HA,
              "inputs": {n: sha256(run / f"{n}.csv") for n in TABLES} | {"metadata.toml": sha256(run / "metadata.toml"),
                                                                       "projection_inputs_metadata": sha256(PROJ / "metadata.json")},
              "script_sha256": sha256(Path(__file__)),
              "outputs": {n: sha256(p) for n, p in paths.items() if n != "projection_audit.json"}}
    paths["projection_audit.json"].write_text(json.dumps(report, indent=2, default=str) + "\n")
    pd.set_option("display.width", 250)
    cal = summary[summary.reef_name.isin(CAL_REEFS)]
    print(cal[["ssp", "supply", "reef_name", "period", "outbreak_episodes", "years_in_outbreak", "mean_episode_peak_ha",
               "coral_mean", "years_coral_below_005", "external_share"]].round(2).to_string(index=False))
    print("Replay max:", report["check1_replay_max_abs"], "failures:", report["check2_failures"])
    print("Supply:", json.dumps(supply_check)[:800])
    print("Flags:", len(flags), flags[:6])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    main(parser.parse_args().run.resolve())
