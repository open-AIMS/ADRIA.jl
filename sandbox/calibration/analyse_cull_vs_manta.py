"""Cull catch vs manta tows on the Lizard v0.3 domain (crash_mechanism_protocol.md,
"Preregistered cull-catch vs manta comparison", 2026-10-14).

Source: CSIRO workbook `260529-COTS-Manta-Cull-RHIS-Lawrence-CSIRO.xlsx` (data licence;
outputs stay in the git-ignored sandbox/data/derived/cull_vs_manta_lizard/; only summary
statistics go into the protocol).

Pairing (cull and surveillance vessels never share a voyage, so pairing is spatial-temporal):
- a cull visit = all cull dives at one cull site in one voyage (vessel x voyage), dated by the
  median dive date; cull CPUE = COTS culled per diver-hour of bottom time;
- paired manta = Manta-sheet tows whose path overlaps the same cull site (CullSiteName) within
  +/- W days of the visit (W = 30 primary, 14 sensitivity).
Fits (Poisson likelihood on paired manta counts given tows, visit level): manta COTS/tow =
  linear a*C, linear + background a*C + b, power a*C^k + b (k > 1: manta hyperdepleted
relative to cull catch; k < 1: hyperstable). Secondary: LTMP manta reef-years (reef_cots.csv,
reporting year) vs reef-year cull CPUE.

    python analyse_cull_vs_manta.py [--output <dir>]
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import gammaln

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "sandbox" / "data"
WORKBOOK = DATA / "260529-COTS-Manta-Cull-RHIS-Lawrence-CSIRO.xlsx"
SITES = DATA / "Lizard_Historical_v0.3" / "site_to_reef.csv"
LTMP = DATA / "reef_cots.csv"
DEFAULT_OUT = DATA / "derived" / "cull_vs_manta_lizard"
COHORTS = ["Cohort1", "Cohort2", "Cohort3", "Cohort4"]
WINDOWS = (30, 14)
LTMP_REEFS = {"14-114": "Macgillivray Reef", "14-118": "Eyrie Reef", "14-143": "North Direction Island",
              "14-116b": "Lizard Isles", "14-116c": "Lizard Isles", "14-116d": "Lizard Isles"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def poisson_fit(counts, tows, cull, model: str) -> dict:
    counts, tows, cull = (np.asarray(v, dtype=float) for v in (counts, tows, cull))

    def rate(theta):
        if model == "linear":
            return np.exp(theta[0]) * cull
        if model == "linear_background":
            return np.exp(theta[0]) * cull + np.exp(theta[1])
        return np.exp(theta[0]) * cull ** np.exp(theta[2]) + np.exp(theta[1])

    def nll(theta):
        mu = np.maximum(rate(theta) * tows, 1e-12)
        return float(-(counts * np.log(mu) - mu - gammaln(counts + 1)).sum())

    starts = {"linear": [np.log(0.02)], "linear_background": [np.log(0.02), np.log(0.02)],
              "power": [np.log(0.02), np.log(0.02), 0.0]}[model]
    best = min((minimize(nll, np.array(starts) + shift, method="Nelder-Mead",
                         options={"maxiter": 20000, "xatol": 1e-8, "fatol": 1e-8})
                for shift in (0.0, 0.5, -0.5)), key=lambda r: r.fun)
    k = len(starts)
    params = {"a": float(np.exp(best.x[0]))}
    if model != "linear":
        params["b"] = float(np.exp(best.x[1]))
    if model == "power":
        params["k"] = float(np.exp(best.x[2]))
    return {"model": model, "params": params, "nll": float(best.fun), "aic": float(2 * k + 2 * best.fun),
            "converged": bool(best.success)}


def pair_visits(visits: pd.DataFrame, manta: pd.DataFrame, window: int) -> pd.DataFrame:
    groups = manta[manta.CullSiteName.notna()].groupby("CullSiteName")
    rows = []
    for v in visits.itertuples(index=False):
        if v.CullSiteName not in groups.groups:
            continue
        t = groups.get_group(v.CullSiteName)
        t = t[(t.date - v.date).abs() <= pd.Timedelta(days=window)]
        if len(t):
            rows.append({"CullSiteName": v.CullSiteName, "reef_id": v.reef_id, "visit_date": v.date.date().isoformat(),
                         "year": v.date.year, "cull_cots": v.culled, "cull_hours": v.hours,
                         "cull_cpue_hr": v.culled / v.hours, "manta_cots": int(t.CrownOfThornsStarfishCount.sum()),
                         "manta_tows": len(t)})
    return pd.DataFrame(rows).assign(manta_cpue=lambda d: d.manta_cots / d.manta_tows)


def main(out: Path) -> None:
    if out.exists():
        raise FileExistsError(f"Refusing to overwrite {out}")
    reefs = set(pd.read_csv(SITES).reef_id)
    manta = pd.read_excel(WORKBOOK, sheet_name="Manta")
    manta = manta[manta.ReefLabel.isin(reefs)].copy()
    manta["date"] = pd.to_datetime(manta.SurveyTime.astype(str).str[:10], format="%Y-%m-%d")
    cull = pd.read_excel(WORKBOOK, sheet_name="Cull")
    cull = cull[cull.ReefLabel.isin(reefs)].copy()
    cull["date"] = pd.to_datetime(cull.SurveyDate.astype(str).str[:10], format="%Y-%m-%d")
    cull["culled"] = cull[COHORTS].sum(axis=1)
    visits = (cull.groupby(["CullSiteName", "VesselName_ID", "VoyageTitle"])
              .agg(reef_id=("ReefLabel", "first"), date=("date", "median"), culled=("culled", "sum"),
                   minutes=("Bottomtime", "sum"), dives=("culled", "size")).reset_index())
    visits = visits[visits.minutes > 0].assign(hours=lambda d: d.minutes / 60)

    out.mkdir(parents=True)
    summary = {"visits": int(len(visits)), "windows": {}}
    for window in WINDOWS:
        p = pair_visits(visits, manta, window)
        p.to_csv(out / f"pairs_w{window:02d}.csv", index=False)
        p["bin"] = pd.qcut(p.cull_cpue_hr.rank(method="first"), 6, labels=False)
        bins = p.groupby("bin").agg(cull_cpue_hr_median=("cull_cpue_hr", "median"),
                                    manta_cpue_mean=("manta_cots", lambda x: x.sum() / p.loc[x.index, "manta_tows"].sum()),
                                    manta_zero_share=("manta_cots", lambda x: float((x == 0).mean())),
                                    visits=("manta_tows", "size"), tows=("manta_tows", "sum"))
        bins.to_csv(out / f"bins_w{window:02d}.csv")
        fits = [poisson_fit(p.manta_cots, p.manta_tows, p.cull_cpue_hr, m)
                for m in ("linear", "linear_background", "power")]
        positive = p[p.cull_cots > 0]
        summary["windows"][f"pm{window}d"] = {
            "pairs": int(len(p)), "manta_tows": int(p.manta_tows.sum()),
            "share_cull_positive_manta_zero": float((positive.manta_cots == 0).mean()),
            "bins": json.loads(bins.reset_index().to_json(orient="records")),
            "fits": fits, "best_by_aic": min(fits, key=lambda f: f["aic"])["model"]}

    # Secondary: LTMP manta reef-years vs reef x reporting-year cull CPUE.
    cull["report_year"] = cull.date.dt.year + (cull.date.dt.month >= 9).astype(int)
    c = cull[cull.ReefLabel.isin(LTMP_REEFS)].assign(ltmp=lambda d: d.ReefLabel.map(LTMP_REEFS))
    reef_year = c.groupby(["ltmp", "report_year"]).agg(culled=("culled", "sum"), minutes=("Bottomtime", "sum")).reset_index()
    reef_year["cull_cpue_hr"] = reef_year.culled / (reef_year.minutes / 60)
    obs = pd.read_csv(LTMP)
    ltmp = reef_year.merge(obs[["reef_name", "year", "cotsptow", "tows"]], left_on=["ltmp", "report_year"],
                           right_on=["reef_name", "year"], how="left")
    ltmp.to_csv(out / "ltmp_reef_year.csv", index=False)
    summary["ltmp_reef_years"] = json.loads(ltmp[["ltmp", "report_year", "culled", "cull_cpue_hr", "cotsptow"]]
                                            .round(3).to_json(orient="records"))

    p = pd.read_csv(out / "pairs_w30.csv")
    fig, ax = plt.subplots(figsize=(7.5, 5))
    ax.scatter(p.cull_cpue_hr, p.manta_cpue, s=10, alpha=.35, color="#3987e5", linewidths=0, label="cull visit (±30 d)")
    b = pd.read_csv(out / "bins_w30.csv")
    ax.plot(b.cull_cpue_hr_median, b.manta_cpue_mean, "o-", color="#0b0b0b", label="sextile means")
    grid = np.linspace(0, p.cull_cpue_hr.quantile(.98), 200)
    for f in summary["windows"]["pm30d"]["fits"]:
        a, bb, k = f["params"]["a"], f["params"].get("b", 0.0), f["params"].get("k", 1.0)
        ax.plot(grid, a * grid ** k + bb, linewidth=1.4, label=f"{f['model']} (AIC {f['aic']:.0f})")
    ax.set_xlim(0, p.cull_cpue_hr.quantile(.98))
    ax.set_ylim(0, max(.6, b.manta_cpue_mean.max() * 1.5))
    ax.set_xlabel("Cull catch (COTS per diver-hour)")
    ax.set_ylabel("Manta COTS per tow over the same cull site")
    ax.grid(color="#e4e3df", linewidth=.8)
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(out / "cull_vs_manta.png", dpi=160, facecolor="#fcfcfb")
    plt.close(fig)

    meta = {"status": "derived_analysis_not_promoted",
            "licence": "CSIRO workbook under data licence: derived outputs are internal (git-ignored), never published",
            "source_sha256": sha256(WORKBOOK), "site_to_reef_sha256": sha256(SITES), "ltmp_sha256": sha256(LTMP),
            "script_sha256": sha256(Path(__file__)), "summary": summary,
            "outputs": {f.name: sha256(f) for f in sorted(out.iterdir()) if f.name != "metadata.json"}}
    (out / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk != "bins"} for k, v in summary["windows"].items()}, indent=1))
    print(ltmp[["ltmp", "report_year", "culled", "cull_cpue_hr", "cotsptow"]].round(3).to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    main(parser.parse_args().output.resolve())
