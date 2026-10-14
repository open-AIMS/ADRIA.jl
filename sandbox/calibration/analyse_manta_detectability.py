"""Paired evidence on manta-tow detectability against diver-based density
(crash_mechanism_protocol.md, 2026-10-14). Lizard Island region only.

Diver densities (independent of manta):
  * SALAD 2021-2026 (Pratchett lab master workbook, sheet Tracks): reef-year density =
    sum(COTS) / sum(area) (recorded) and sum(COTS + scars) / sum(area) (inferred).
    "Eagle" in SALAD is Eyrie Reef.
  * Pratchett (2005, Coral Reefs 24:453-462) Lizard belt transects (200 m^2, searched
    twice, <0.8% missed), digitised from Fig. 2 and the text (per 200 m^2 x 50 = per ha).
    The two open-circle points (Sweatman et al. 1998) are excluded as not independent.
Manta CPUE (COTS per tow):
  * LTMP (sandbox/data/reef_cots.csv), dated by survey_date (NOT the reporting year).
  * Control Program / EOTR manta (Lawrence/CSIRO workbook, sheet Manta), per reef-year.
Cull catch rate (COTS per hour of bottom time; sheet Cull), reported, not fitted.

Pairing: same reef; a diver estimate is paired with a manta survey on that reef within
PAIR_WINDOW_MONTHS. Models for manta counts (negative binomial, log-likelihood on counts
with tows as exposure), diver density A as the true-density proxy:
  linear      mu = q A
  saturating  mu = q A^2 / (A + A_half)          (the preregistered detection curve)
  power       mu = a A^b                          (b > 1: hyperdepletion)

Derived tables go to the git-ignored sandbox/data/derived (licensed/internal sources);
the fit summary and plot go to a new run directory.

Run with base miniconda Python:
    python analyse_manta_detectability.py <run_dir>
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
SALAD = DATA / "COTS SALAD MASTER 2021-2026.xlsx"
CSIRO = DATA / "260529-COTS-Manta-Cull-RHIS-Lawrence-CSIRO.xlsx"
LTMP = DATA / "reef_cots.csv"
DERIVED_ROOT = DATA / "derived" / "manta_detectability_pairs_lizard"
PAIR_WINDOW_MONTHS = 6
REEFS = {  # canonical -> (SALAD name, LTMP name, CSIRO ReefName regex)
    "Lizard": ("Lizard", "Lizard Isles", r"^Lizard Island Reef"),
    "Eyrie": ("Eagle", "Eyrie Reef", r"^Eyrie Reef"),
    "MacGillivray": ("MacGillivray", "Macgillivray Reef", r"^MacGillivray Reef"),
    "North Direction": ("North Direction", "North Direction Island", r"^North Direction Reef"),
    "South Direction": ("South Direction", "South Direction Island", r"^South Direction Reef"),
    "Martin": ("Martin", "Martin Reef", r"^Martin Reef"),
}
# Pratchett 2005 Fig. 2 / text: (decimal date, starfish per 200 m^2). Digitised (+-0.05).
PRATCHETT_2005 = [(1996.79, 0.80), (1996.96, 1.05), (1997.12, 0.96), (1997.46, 1.01),
                  (1997.87, 0.86), (1998.12, 1.10), (1998.46, 0.92), (1998.87, 0.27),
                  (1999.04, 0.06), (1999.71, 0.01)]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def decimal_date(series: pd.Series) -> pd.Series:
    d = pd.to_datetime(series, errors="coerce")
    return d.dt.year + (d.dt.dayofyear - 1) / 365.25


def diver_estimates() -> pd.DataFrame:
    tracks = pd.read_excel(SALAD, "Tracks")
    tracks = tracks[tracks.Region.astype(str).str.contains("Lizard", case=False)].copy()
    tracks["date"] = decimal_date(tracks.Date)
    rows = []
    for reef, (salad, _, _) in REEFS.items():
        x = tracks[tracks.Reef == salad].dropna(subset=["date"])
        for _, g in x.groupby(x.date.round(0) if False else x.Date.astype(str).str[:7]):
            area_ha = g["AREA(5MWIDE)"].sum() / 1e4
            if area_ha <= 0:
                continue
            rows.append({"reef": reef, "method": "SALAD", "date": float(g.date.mean()),
                         "area_ha": area_ha, "cots": int(g["No. COTS"].sum()),
                         "scars": int(g["No. scars"].sum()),
                         "density_recorded_ha": g["No. COTS"].sum() / area_ha,
                         "density_inferred_ha": (g["No. COTS"].sum() + g["No. scars"].sum()) / area_ha})
    for date, per200 in PRATCHETT_2005:
        rows.append({"reef": "Lizard", "method": "Pratchett2005_belt", "date": date, "area_ha": np.nan,
                     "cots": np.nan, "scars": np.nan, "density_recorded_ha": per200 * 50.0,
                     "density_inferred_ha": per200 * 50.0})
    return pd.DataFrame(rows)


def manta_surveys() -> pd.DataFrame:
    ltmp = pd.read_csv(LTMP)
    rows = []
    for reef, (_, name, _) in REEFS.items():
        for r in ltmp[ltmp.reef_name == name].itertuples(index=False):
            rows.append({"reef": reef, "source": "LTMP", "date": float(r.survey_date),
                         "tows": int(r.tows), "cots": int(r.cots), "report_year": int(r.year)})
    manta = pd.read_excel(CSIRO, "Manta", usecols=["SurveyTime", "ReefName", "CrownOfThornsStarfishCount"])
    manta["date"] = decimal_date(manta.SurveyTime)
    manta["month"] = pd.to_datetime(manta.SurveyTime, errors="coerce").dt.strftime("%Y-%m")
    for reef, (_, _, pattern) in REEFS.items():
        x = manta[manta.ReefName.astype(str).str.contains(pattern, regex=True)]
        for month, g in x.groupby("month"):
            rows.append({"reef": reef, "source": "ControlProgram", "date": float(g.date.mean()),
                         "tows": int(len(g)), "cots": int(g.CrownOfThornsStarfishCount.sum()), "report_year": np.nan})
    return pd.DataFrame(rows)


def cull_rates() -> pd.DataFrame:
    cull = pd.read_excel(CSIRO, "Cull", usecols=["SurveyDate", "ReefName", "Bottomtime",
                                                  "Cohort1", "Cohort2", "Cohort3", "Cohort4"])
    cull["date"] = decimal_date(cull.SurveyDate)
    cull["n"] = cull[["Cohort1", "Cohort2", "Cohort3", "Cohort4"]].sum(axis=1)
    rows = []
    for reef, (_, _, pattern) in REEFS.items():
        x = cull[cull.ReefName.astype(str).str.contains(pattern, regex=True)]
        for year, g in x.groupby(np.floor(x.date)):
            minutes = g.Bottomtime.sum()
            if minutes > 0:
                rows.append({"reef": reef, "year": int(year), "dives": len(g), "cots": int(g.n.sum()),
                             "cots_per_hour": g.n.sum() / (minutes / 60.0)})
    return pd.DataFrame(rows)


def pair(diver: pd.DataFrame, manta: pd.DataFrame) -> pd.DataFrame:
    """One pair per manta survey: the diver estimate on the same reef nearest in time
    (within the window), so no manta count is used twice."""
    rows = []
    window = PAIR_WINDOW_MONTHS / 12.0
    for m in manta.itertuples(index=False):
        d = diver[(diver.reef == m.reef) & ((diver.date - m.date).abs() <= window)]
        if d.empty:
            continue
        nearest = d.iloc[(d.date - m.date).abs().argsort()[:1]].iloc[0]
        rows.append({"reef": m.reef, "diver_method": nearest.method, "manta_source": m.source,
                     "diver_date": nearest.date, "manta_date": m.date,
                     "gap_months": 12 * (m.date - nearest.date),
                     "density_recorded_ha": nearest.density_recorded_ha,
                     "density_inferred_ha": nearest.density_inferred_ha,
                     "tows": int(m.tows), "manta_cots": int(m.cots), "manta_cpue": m.cots / m.tows})
    return pd.DataFrame(rows).sort_values(["diver_method", "reef", "manta_date"]).reset_index(drop=True)


def nb_nll(mu: np.ndarray, counts: np.ndarray, k: float) -> float:
    mu = np.maximum(mu, 1e-12)
    return -float(np.sum(gammaln(counts + k) - gammaln(k) - gammaln(counts + 1) +
                         k * np.log(k / (k + mu)) + counts * np.log(mu / (k + mu))))


def fit(pairs: pd.DataFrame, density: str) -> dict:
    A = pairs[density].to_numpy(float)
    n = pairs.tows.to_numpy(float)
    y = pairs.manta_cots.to_numpy(float)
    models = {
        "linear": (lambda p: n * np.exp(p[0]) * A, [np.log(0.015)]),
        "saturating": (lambda p: n * np.exp(p[0]) * A * A / (A + np.exp(p[1])), [np.log(0.02), np.log(5.0)]),
        "power": (lambda p: n * np.exp(p[0]) * A ** p[1], [np.log(0.01), 1.0]),
    }
    out = {}
    for name, (mean, start) in models.items():
        objective = lambda p: nb_nll(mean(p[:-1]), y, np.exp(p[-1]))
        best = None
        for k0 in (np.log(0.3), np.log(1.0), np.log(5.0)):
            r = minimize(objective, np.array(start + [k0]), method="Nelder-Mead",
                         options={"maxiter": 20000, "xatol": 1e-9, "fatol": 1e-10})
            best = r if best is None or r.fun < best.fun else best
        p = best.x
        params = {"linear": {"q": float(np.exp(p[0]))},
                  "saturating": {"q_max": float(np.exp(p[0])), "A_half": float(np.exp(p[1]))},
                  "power": {"a": float(np.exp(p[0])), "b": float(p[1])}}[name]
        params["nb_k"] = float(np.exp(p[-1]))
        out[name] = {"params": params, "nll": float(best.fun), "n_params": len(p),
                     "aic": float(2 * best.fun + 2 * len(p))}
    # Profile likelihood 95% interval for A_half (saturating).
    grid = np.exp(np.linspace(np.log(0.05), np.log(500), 120))
    prof = []
    for h in grid:
        objective = lambda p: nb_nll(n * np.exp(p[0]) * A * A / (A + h), y, np.exp(p[1]))
        r = min((minimize(objective, np.array([np.log(0.02), k0]), method="Nelder-Mead") for k0 in (np.log(0.3), np.log(2.0))),
                key=lambda r: r.fun)
        prof.append(r.fun)
    prof = np.array(prof)
    inside = grid[prof - prof.min() <= 1.92]
    out["saturating"]["A_half_profile_95"] = [float(inside.min()), float(inside.max())]
    out["saturating"]["A_half_profile_hits_grid_bound"] = bool(inside.min() <= grid[0] or inside.max() >= grid[-1])
    return out


def main(run: Path) -> None:
    DERIVED = DERIVED_ROOT / run.name
    if run.exists() or DERIVED.exists():
        raise FileExistsError("Refusing to overwrite outputs")
    diver, manta, cull = diver_estimates(), manta_surveys(), cull_rates()
    pairs = pair(diver, manta)
    DERIVED.mkdir(parents=True)
    run.mkdir(parents=True)
    for name, table in (("diver_estimates.csv", diver), ("manta_surveys.csv", manta),
                        ("cull_rates.csv", cull), ("pairs.csv", pairs)):
        table.to_csv(DERIVED / name, index=False)
    fits = {}
    for subset, rows in (("all", pairs), ("LTMP", pairs[pairs.manta_source == "LTMP"]),
                         ("ControlProgram", pairs[pairs.manta_source == "ControlProgram"]),
                         ("SALAD_only", pairs[pairs.diver_method == "SALAD"])):
        if len(rows) >= 4:
            fits[subset] = {density: fit(rows, density) for density in
                            ("density_recorded_ha", "density_inferred_ha")} | {"n_pairs": int(len(rows))}

    fig, axis = plt.subplots(figsize=(8, 6))
    styles = {("SALAD", "LTMP"): ("#2a78d6", "o"), ("SALAD", "ControlProgram"): ("#86b6ef", "s"),
              ("Pratchett2005_belt", "LTMP"): ("#eb6834", "D"),
              ("Pratchett2005_belt", "ControlProgram"): ("#eb6834", "^")}
    for (method, source), group in pairs.groupby(["diver_method", "manta_source"]):
        color, marker = styles.get((method, source), ("#52514e", "o"))
        axis.scatter(group.density_inferred_ha, group.manta_cpue.clip(lower=2e-3), color=color, marker=marker,
                     s=40, edgecolor="white", linewidth=1, label=f"{method} vs {source} manta (n={len(group)})", zorder=3)
    a = np.logspace(-1, 2.3, 200)
    axis.plot(a, 0.015 * a, color="#52514e", linestyle="--", linewidth=1.5, label="frozen q = 0.015")
    if "all" in fits:
        s = fits["all"]["density_inferred_ha"]["saturating"]["params"]
        axis.plot(a, s["q_max"] * a * a / (a + s["A_half"]), color="#0d366b", linewidth=2,
                  label=f"fit: q_max {s['q_max']:.3f}, A½ {s['A_half']:.1f} /ha")
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlabel("Diver density (COTS/ha; SALAD inferred incl. scars, or belt transects)")
    axis.set_ylabel("Manta COTS/tow (zeros at 0.002)")
    axis.grid(color="#e4e3df", linewidth=.8)
    axis.legend(fontsize=7, frameon=False, loc="upper left")
    axis.set_title("Lizard region: manta CPUE against paired diver density", loc="left", fontsize=10)
    fig.tight_layout()
    fig.savefig(run / "manta_vs_diver_density.png", dpi=170, facecolor="#fcfcfb")
    plt.close(fig)

    report = {
        "status": "evidence_analysis_not_promoted",
        "pair_window_months": PAIR_WINDOW_MONTHS,
        "n_pairs_by_source": pairs.groupby(["diver_method", "manta_source"]).size().rename("n").reset_index().to_dict(orient="records"),
        "fits": fits,
        "cull_rates_summary": cull.groupby("reef").agg(years=("year", "nunique"), cots=("cots", "sum"),
                                                       mean_cots_per_hour=("cots_per_hour", "mean")).reset_index().to_dict(orient="records"),
        "derived_dir": str(DERIVED),
        "inputs": {p.name: sha256(p) for p in (SALAD, CSIRO, LTMP, Path(__file__))},
        "pairing": "one pair per manta survey: nearest diver estimate on the same reef within the window",
        "outputs": {name: sha256(DERIVED / name) for name in ("diver_estimates.csv", "manta_surveys.csv", "cull_rates.csv", "pairs.csv")}
                   | {"manta_vs_diver_density.png": sha256(run / "manta_vs_diver_density.png")},
        "caveats": ["SALAD searches follow feeding scars (targeted, not random); densities are reef-section, not reef-wide",
                    "culling between surveys lowers density; pairs within 6 months, not same day",
                    "Pratchett 2005 values digitised from Fig. 2 (+-0.05 per 200 m^2)",
                    "diver detectability <1 (SALAD ~0.5 recorded; belt ~0.7-0.9); inferred density partly corrects SALAD"],
    }
    (run / "manta_detectability_fit.json").write_text(json.dumps(report, indent=2, default=float) + "\n")
    pd.set_option("display.width", 220)
    print(pairs.round(3).to_string(index=False))
    print(json.dumps(fits, indent=1, default=float))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    main(parser.parse_args().run.resolve())
