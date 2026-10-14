"""Diver-anchored COTS density targets for the 2020s build-up (lead's decisions, 2026-10-14,
seventh round: the model must match the fifth-outbreak build-up 2020-2026 seen in SALAD and
cull data, which are more accurate than manta at low density).

Sources (licensed / internal; outputs go to the git-ignored
sandbox/data/derived/density_targets_2020s/; only summaries go into the protocol):
- SALAD 2021-2026 (Pratchett lab master workbook, sheet Tracks, Lizard region): diver swims
  with recorded area (5 m wide); density = sum(COTS) / sum(area) per reef x survey month.
  "Eagle" is Eyrie Reef. Densities are as recorded (SALAD itself misses some COTS).
- Cull dives (CSIRO workbook, sheet Cull): catch per diver-hour per reef x month.

Anchor: SALAD reef-months paired with cull on the same reef within +/- WINDOW_DAYS
(all cull dives on that reef in the window, pooled). Fit SALAD counts ~ Poisson(area_ha x c x
CPUE^k) (k = 1 and free k), giving density (COTS ha^-1) from cull catch per diver-hour.

Targets: per reef x LTMP reporting year (Sep-Dec -> year + 1), 2019-2026: SALAD density
(area-weighted) where surveyed, else the cull-derived density, with the source flagged.

    python build_2020s_density_targets.py [--output <dir>]
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import gammaln

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "sandbox" / "data"
SALAD = DATA / "COTS SALAD MASTER 2021-2026.xlsx"
CSIRO = DATA / "260529-COTS-Manta-Cull-RHIS-Lawrence-CSIRO.xlsx"
DEFAULT_OUT = DATA / "derived" / "density_targets_2020s"
WINDOW_DAYS = 60
# SALAD reef name -> (target reef, GBRMPA reef IDs in the domain)
REEFS = {
    "Lizard": ("Lizard Island Reef", ["14-116b", "14-116c", "14-116d"]),
    "Eagle": ("Eyrie Reef", ["14-118"]),
    "MacGillivray": ("MacGillivray Reef", ["14-114"]),
    "North Direction": ("North Direction Reef", ["14-143"]),
    "South Direction": ("South Direction Reef", ["14-147a", "14-147b"]),
    "Martin": ("Martin Reef", ["14-123"]),
}
COHORTS = ["Cohort1", "Cohort2", "Cohort3", "Cohort4"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def report_year(dates: pd.Series) -> pd.Series:
    return dates.dt.year + (dates.dt.month >= 9).astype(int)


def fit(counts, area, cpue, free_k: bool) -> dict:
    counts, area, cpue = (np.asarray(v, float) for v in (counts, area, cpue))

    def nll(theta):
        k = np.exp(theta[1]) if free_k else 1.0
        mu = np.maximum(area * np.exp(theta[0]) * cpue ** k, 1e-12)
        return float(-(counts * np.log(mu) - mu - gammaln(counts + 1)).sum())

    start = [np.log(3.0), 0.0] if free_k else [np.log(3.0)]
    r = min((minimize(nll, np.array(start) + s, method="Nelder-Mead",
                      options={"maxiter": 20000, "xatol": 1e-9, "fatol": 1e-9}) for s in (0.0, 0.7, -0.7)),
            key=lambda r: r.fun)
    out = {"c_density_per_cpue": float(np.exp(r.x[0])), "k": float(np.exp(r.x[1])) if free_k else 1.0,
           "nll": float(r.fun), "aic": float(2 * len(start) + 2 * r.fun), "converged": bool(r.success)}
    # Profile 95% interval for c (k fixed at its estimate).
    k = out["k"]
    grid = np.exp(np.linspace(np.log(out["c_density_per_cpue"]) - 2, np.log(out["c_density_per_cpue"]) + 2, 801))
    prof = [-(counts * np.log(np.maximum(area * g * cpue ** k, 1e-12)) - area * g * cpue ** k - gammaln(counts + 1)).sum()
            for g in grid]
    ok = grid[np.array(prof) <= r.fun + 1.92]
    out["c_profile_95"] = [float(ok.min()), float(ok.max())]
    return out


def main(out: Path) -> None:
    if out.exists():
        raise FileExistsError(f"Refusing to overwrite {out}")
    tracks = pd.read_excel(SALAD, "Tracks")
    tracks = tracks[tracks.Region.astype(str).str.contains("Lizard", case=False) & tracks.Reef.isin(REEFS)].copy()
    tracks["date"] = pd.to_datetime(tracks.Date, errors="coerce")
    tracks = tracks.dropna(subset=["date"])
    tracks = tracks[tracks["AREA(5MWIDE)"] > 0]
    tracks["month"] = tracks.date.dt.to_period("M")
    salad = (tracks.groupby(["Reef", "month"])
             .agg(date=("date", "median"), tracks=("Track", "size"), area_ha=("AREA(5MWIDE)", lambda a: a.sum() / 1e4),
                  cots=("No. COTS", "sum"), scars=("No. scars", "sum")).reset_index())
    salad["density_ha"] = salad.cots / salad.area_ha
    salad["target_reef"] = salad.Reef.map(lambda r: REEFS[r][0])

    cull = pd.read_excel(CSIRO, "Cull")
    cull["date"] = pd.to_datetime(cull.SurveyDate.astype(str).str[:10], format="%Y-%m-%d")
    cull["culled"] = cull[COHORTS].sum(axis=1)
    label_to_reef = {rid: REEFS[name][0] for name in REEFS for rid in REEFS[name][1]}
    cull = cull[cull.ReefLabel.isin(label_to_reef) & (cull.Bottomtime > 0)].copy()
    cull["target_reef"] = cull.ReefLabel.map(label_to_reef)

    pairs = []
    for s in salad.itertuples(index=False):
        c = cull[(cull.target_reef == s.target_reef) & ((cull.date - s.date).abs() <= pd.Timedelta(days=WINDOW_DAYS))]
        if len(c):
            pairs.append({"target_reef": s.target_reef, "salad_month": str(s.month), "salad_area_ha": s.area_ha,
                          "salad_cots": int(s.cots), "salad_density_ha": s.density_ha, "cull_dives": len(c),
                          "cull_hours": c.Bottomtime.sum() / 60, "cull_cots": int(c.culled.sum()),
                          "cull_cpue_hr": c.culled.sum() / (c.Bottomtime.sum() / 60)})
    pairs = pd.DataFrame(pairs)
    usable = pairs[pairs.cull_cpue_hr > 0]
    fits = {"linear": fit(usable.salad_cots, usable.salad_area_ha, usable.cull_cpue_hr, False),
            "power": fit(usable.salad_cots, usable.salad_area_ha, usable.cull_cpue_hr, True)}
    best = min(fits.values(), key=lambda f: f["aic"])

    # Targets by reef x reporting year.
    salad["report_year"] = report_year(salad.date)
    s_year = (salad.groupby(["target_reef", "report_year"])
              .agg(salad_cots=("cots", "sum"), salad_area_ha=("area_ha", "sum"), salad_tracks=("tracks", "sum"))
              .reset_index())
    s_year["salad_density_ha"] = s_year.salad_cots / s_year.salad_area_ha
    cull["report_year"] = report_year(cull.date)
    c_year = (cull.groupby(["target_reef", "report_year"])
              .agg(cull_cots=("culled", "sum"), cull_hours=("Bottomtime", lambda m: m.sum() / 60), cull_dives=("culled", "size"))
              .reset_index())
    c_year["cull_cpue_hr"] = c_year.cull_cots / c_year.cull_hours
    c_year["cull_density_ha"] = best["c_density_per_cpue"] * c_year.cull_cpue_hr ** best["k"]
    targets = s_year.merge(c_year, on=["target_reef", "report_year"], how="outer").sort_values(["target_reef", "report_year"])
    targets = targets[targets.report_year.between(2019, 2026)]
    targets["target_density_ha"] = targets.salad_density_ha.fillna(targets.cull_density_ha)
    targets["target_source"] = np.where(targets.salad_density_ha.notna(), "SALAD", "cull_converted")

    out.mkdir(parents=True)
    salad.assign(month=salad.month.astype(str)).to_csv(out / "salad_reef_month.csv", index=False)
    pairs.to_csv(out / "salad_cull_pairs.csv", index=False)
    targets.to_csv(out / "targets_reef_year.csv", index=False)
    meta = {"status": "derived_target_not_promoted",
            "licence": "SALAD and CSIRO workbooks are internal/licensed: derived outputs git-ignored, never published",
            "sources": {"salad_sha256": sha256(SALAD), "csiro_sha256": sha256(CSIRO)},
            "pairing": f"same reef, cull within +/-{WINDOW_DAYS} days of the SALAD survey month median date",
            "pairs": int(len(pairs)), "pairs_with_cull_catch": int(len(usable)), "fits": fits,
            "selected_fit": "power" if best is fits["power"] else "linear",
            "units": "density COTS ha^-1 (SALAD as recorded); cull CPUE COTS per diver-hour; reporting year Sep-Dec -> +1",
            "script_sha256": sha256(Path(__file__)),
            "outputs": {f.name: sha256(f) for f in sorted(out.iterdir()) if f.suffix == ".csv"}}
    (out / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
    pd.set_option("display.width", 250)
    print(pairs.round(2).to_string(index=False))
    print(json.dumps(fits, indent=1))
    print(targets[["target_reef", "report_year", "salad_density_ha", "cull_cpue_hr", "cull_density_ha",
                   "target_density_ha", "target_source"]].round(2).to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    main(parser.parse_args().output.resolve())
