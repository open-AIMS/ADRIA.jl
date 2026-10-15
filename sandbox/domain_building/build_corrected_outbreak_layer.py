"""Correct the GBR outbreak-probability prediction layer with diver-anchored evidence, 2018-2026
(lead's decisions, 2026-10-14/15: correct the prediction layer by spatial interpolation of
cull and SALAD residuals; 15 COTS/ha is the diver-equivalent of 0.22 COTS per manta tow).

Inputs (licensed / internal; outputs go to the git-ignored sandbox/data/derived):
- predictions: sandbox/data/gbrPredsAdj_20262408.csv (1991-2025) and
  gbrPreds2026Adj_20262408.csv (2026); columns year, reefName "<name> (<GBRMPA id>)", outbrProb, x, y;
- cull dives (CSIRO workbook, sheet Cull; GBR-wide), reporting year (Sep-Dec -> year + 1);
- SALAD diver densities (build_2020s_density_targets.py output, Lizard region).

Outbreak indicator I at an observed reef-year:
- SALAD density >= 15 COTS/ha where surveyed;
- otherwise cull catch >= C15 COTS per diver-hour with >= MIN_HOURS of bottom time, where
  C15 = (15 / 7.874)^(1 / 0.549) from the SALAD anchor (density = 7.874 x CPUE^0.549).
Residual r = I - p_pred. Per year, residuals are spread to every prediction reef by
inverse-distance weighting within R_KM, shrunk toward zero:
    r_hat = sum(w_i r_i) / (sum(w_i) + w0),  w_i = 1 / max(d_i, 1 km)^2,  w0 = 1 / R_KM^2
so a single observation at distance R contributes half its residual and reefs beyond R are
unchanged. Corrected p = clamp(p + r_hat, 0, 1). Years before 2018 are copied unchanged.

    python build_corrected_outbreak_layer.py [--output <dir>]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "sandbox" / "data"
PREDS = DATA / "gbrPredsAdj_20262408.csv"
PREDS_2026 = DATA / "gbrPreds2026Adj_20262408.csv"
CSIRO = DATA / "260529-COTS-Manta-Cull-RHIS-Lawrence-CSIRO.xlsx"
TARGETS = DATA / "derived" / "density_targets_2020s" / "targets_reef_year.csv"
SALAD_IDS = {"Lizard Island Reef": ["14-116b", "14-116c", "14-116d"], "Eyrie Reef": ["14-118"],
             "MacGillivray Reef": ["14-114"], "North Direction Reef": ["14-143"],
             "South Direction Reef": ["14-147a", "14-147b"], "Martin Reef": ["14-123"]}
DEFAULT_OUT = DATA / "derived" / "outbreak_layer_corrected_2018_2026"
YEARS = range(2018, 2027)
DENSITY_THRESHOLD = 15.0
ANCHOR_C, ANCHOR_K = 7.874, 0.549
C15 = (DENSITY_THRESHOLD / ANCHOR_C) ** (1 / ANCHOR_K)
MIN_HOURS = 2.0
R_KM = 50.0


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def gbrmpa_id(name: str) -> str | None:
    m = re.search(r"\(([0-9]{2}-[0-9]{3}[a-z]?)\)\s*$", str(name))
    return m.group(1) if m else None


def haversine_km(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = map(np.radians, (lat1, lon1, lat2, lon2))
    a = np.sin((lat2 - lat1) / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2
    return 2 * 6371.0 * np.arcsin(np.sqrt(a))


def main(out: Path) -> None:
    if out.exists():
        raise FileExistsError(f"Refusing to overwrite {out}")
    preds = pd.concat([pd.read_csv(PREDS), pd.read_csv(PREDS_2026)[["year", "reefName", "outbrProb", "x", "y"]]],
                      ignore_index=True)
    if preds.duplicated(["year", "reefName"]).any():
        raise ValueError("Duplicate prediction reef-years")
    preds["reef_id"] = preds.reefName.map(gbrmpa_id)

    cull = pd.read_excel(CSIRO, sheet_name="Cull")
    cull["date"] = pd.to_datetime(cull.SurveyDate.astype(str).str[:10], format="%Y-%m-%d")
    cull["report_year"] = cull.date.dt.year + (cull.date.dt.month >= 9).astype(int)
    cull["culled"] = cull[["Cohort1", "Cohort2", "Cohort3", "Cohort4"]].sum(axis=1)
    c = (cull[cull.ReefLabel.notna() & cull.report_year.isin(YEARS)]
         .groupby(["ReefLabel", "report_year"]).agg(culled=("culled", "sum"), hours=("Bottomtime", lambda m: m.sum() / 60))
         .reset_index().rename(columns={"ReefLabel": "reef_id", "report_year": "year"}))
    c = c[c.hours >= MIN_HOURS]
    c["cpue_hr"] = c.culled / c.hours
    c["indicator"] = (c.cpue_hr >= C15).astype(float)
    c["source"] = "cull"
    salad = pd.read_csv(TARGETS)
    salad = salad[salad.target_source.eq("SALAD") & salad.report_year.isin(YEARS)]
    s_rows = [{"reef_id": rid, "year": int(r.report_year), "indicator": float(r.salad_density_ha >= DENSITY_THRESHOLD),
               "source": "SALAD", "density_ha": float(r.salad_density_ha)}
              for r in salad.itertuples(index=False) for rid in SALAD_IDS[r.target_reef]]
    obs = pd.concat([pd.DataFrame(s_rows), c[["reef_id", "year", "indicator", "source", "cpue_hr"]]], ignore_index=True)
    obs = obs.sort_values("source", key=lambda s: s.map({"SALAD": 0, "cull": 1})).drop_duplicates(["reef_id", "year"])
    obs = obs.merge(preds[["year", "reef_id", "outbrProb", "x", "y"]], on=["reef_id", "year"], how="inner")
    obs["residual"] = obs.indicator - obs.outbrProb

    corrected = preds.copy()
    corrected["residual_hat"] = 0.0
    for year in YEARS:
        o = obs[obs.year == year]
        rows = corrected.year == year
        if o.empty or not rows.any():
            continue
        targets = corrected.loc[rows]
        d = haversine_km(targets.y.to_numpy()[:, None], targets.x.to_numpy()[:, None],
                         o.y.to_numpy()[None, :], o.x.to_numpy()[None, :])
        w = np.where(d <= R_KM, 1.0 / np.maximum(d, 1.0) ** 2, 0.0)
        r_hat = (w * o.residual.to_numpy()[None, :]).sum(axis=1) / (w.sum(axis=1) + 1.0 / R_KM ** 2)
        corrected.loc[rows, "residual_hat"] = r_hat
    corrected["outbrProb_original"] = corrected.outbrProb
    corrected["outbrProb"] = np.clip(corrected.outbrProb + corrected.residual_hat, 0.0, 1.0)

    out.mkdir(parents=True)
    corrected[["year", "reefName", "outbrProb", "x", "y"]].to_csv(out / "gbrPredsCorrected_1991_2026.csv", index=False)
    corrected.to_csv(out / "corrected_layer_diagnostics.csv", index=False)
    obs.to_csv(out / "observed_indicators.csv", index=False)
    north = corrected[(corrected.y > -17.5) & (corrected.y < -12)]
    summary = (north.groupby("year")[["outbrProb_original", "outbrProb"]].mean().round(3)
               .join(obs.groupby("year").agg(n_obs=("indicator", "size"), n_outbreak=("indicator", "sum"))))
    meta = {"status": "derived_input_not_promoted",
            "licence": "inputs internal/licensed; outputs git-ignored, never published",
            "indicator": f"SALAD density >= {DENSITY_THRESHOLD}/ha, else cull CPUE >= {C15:.3f}/diver-hour (anchor {ANCHOR_C} x CPUE^{ANCHOR_K}) with >= {MIN_HOURS} h",
            "interpolation": f"IDW 1/d^2 within {R_KM} km, shrinkage w0 = 1/{R_KM}^2; clamp [0, 1]; years {YEARS.start}-{YEARS.stop - 1}",
            "inputs": {p.name: sha256(p) for p in (PREDS, PREDS_2026, CSIRO, TARGETS)},
            "observed_reef_years": int(len(obs)), "observed_outbreak_reef_years": int(obs.indicator.sum()),
            "northern_mean_probability": json.loads(summary.reset_index().to_json(orient="records")),
            "script_sha256": sha256(Path(__file__)),
            "outputs": {f.name: sha256(f) for f in sorted(out.iterdir()) if f.suffix == ".csv"}}
    (out / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(f"C15 = {C15:.3f} COTS per diver-hour")
    print(summary.to_string())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    main(parser.parse_args().output.resolve())
