"""Survey-timing tracker for the calibration reefs (crash_mechanism_protocol.md, 2026-10-14).

The scoring convention is unchanged: observations are compared with model year Y by their
LTMP reporting year Y (the AIMS financial year July(Y-1)-June(Y)). This script records, for
every calibration-reef survey (manta COTS, manta hard coral, photo transects), when it was
actually made, and whether it predates the bleaching summer of the model year it is
compared with. Model year Y applies the DHW of summer Y (peaking Feb-Apr Y), so a survey
made before about May of year Y cannot show year-Y bleaching that the model already includes.

Effective DHW per reef-year uses the domain DHW (scenario 1) and ADRIA's
effective_dhw_at_depth (kappa_base 0.04, mixing scale 12, depth_med), averaged over the
reef's sites.

Run with base miniconda Python (xarray, geopandas):
    python build_survey_timing.py <run_dir>
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "sandbox" / "data"
DOMAIN = DATA / "Lizard_Historical_v0.3"
OUT_DIR = DATA / "derived" / "ltmp_survey_timing"
REEFS = {"Lizard Island Reef": ("Lizard Isles", "14-116"), "MacGillivray Reef": ("Macgillivray Reef", "14-114"),
         "North Direction Reef": ("North Direction Island", "14-143"), "Eyrie Reef": ("Eyrie Reef", "14-118")}
BLEACHING_SEASON_END = 0.33  # decimal year (~end of April): later surveys can show that summer's bleaching
DHW_FLAGS = (3.0, 4.0)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def effective_dhw(surface: np.ndarray, depth: np.ndarray) -> np.ndarray:
    kappa = 0.04 / (1.0 + surface / 12.0)
    return surface * np.exp(-kappa * np.maximum(0.0, depth - 2.0))


def reef_dhw() -> pd.DataFrame:
    ds = xr.open_dataset(DOMAIN / "DHWs" / "dhw_RCPhistorical.nc")
    sites = gpd.read_file(DOMAIN / "spatial" / "lizard_cluster.gpkg")
    dhw = ds.dhw.sel(scenario=1).values
    years = ds.time.values.astype(int)
    rows = []
    for reef, (_, gid) in REEFS.items():
        idx = np.where(sites.GBRMPA_ID.astype(str).str.startswith(gid))[0]
        depth = sites.depth_med.values[idx].astype(float)
        for yi, year in enumerate(years):
            rows.append({"reef_name": reef, "year": int(year),
                         "effective_dhw": float(effective_dhw(dhw[idx, yi].astype(float), depth).mean())})
    return pd.DataFrame(rows)


def surveys() -> pd.DataFrame:
    cots = pd.read_csv(DATA / "reef_cots.csv")
    manta = pd.read_csv(DATA / "reef_manta.csv")
    photo = pd.read_csv(DATA / "reef_photo_transect.csv", low_memory=False)
    rows = []
    for reef, (obs, _) in REEFS.items():
        for r in cots[cots.reef_name == obs].itertuples(index=False):
            rows.append(("cots_manta", reef, int(r.year), float(r.survey_date)))
        m = manta[(manta.reef_name == obs) & (manta.data_type == "manta") & (manta.domain_category == "reef") &
                  (manta.variable == "HC") & (manta.purpose == "MANTA") & (manta.project_code == "LTMP") &
                  (manta.depth == 9)]
        for (year, date), _ in m.groupby(["report_year", "date"]):
            rows.append(("coral_manta", reef, int(year), float(date)))
        p = photo[(photo.reef_name == obs) & (photo.purpose == "COMPOSITION") & (photo.variable == "HARD CORAL") &
                  (photo.depth == 9)]
        for (year, date), _ in p.groupby(["report_year", "date"]):
            rows.append(("photo_transect", reef, int(year), float(date)))
    return pd.DataFrame(rows, columns=["dataset", "reef_name", "report_year", "survey_date"])


def main(run: Path) -> None:
    if run.exists() or OUT_DIR.exists():
        raise FileExistsError("Refusing to overwrite outputs")
    dhw = reef_dhw().set_index(["reef_name", "year"]).effective_dhw
    table = surveys()
    table["survey_month"] = (np.floor((table.survey_date % 1) * 12) + 1).astype(int)
    table["lag_years_vs_model_year"] = table.report_year - table.survey_date  # >0: survey made before model year Y ends
    table["survey_before_model_year_bleaching"] = table.survey_date < table.report_year + BLEACHING_SEASON_END
    table["model_year_effective_dhw"] = [float(dhw.get((r, y), np.nan)) for r, y in zip(table.reef_name, table.report_year)]
    table["prev_year_effective_dhw"] = [float(dhw.get((r, y - 1), np.nan)) for r, y in zip(table.reef_name, table.report_year)]
    for level in DHW_FLAGS:
        tag = f"{level:g}"
        # The model has applied bleaching (DHW >= level) that the survey had not yet seen.
        table[f"model_ahead_bleaching_ge{tag}"] = table.survey_before_model_year_bleaching & (table.model_year_effective_dhw >= level)
    table = table.sort_values(["dataset", "reef_name", "report_year"]).reset_index(drop=True)
    OUT_DIR.mkdir(parents=True)
    run.mkdir(parents=True)
    table.to_csv(OUT_DIR / "survey_timing.csv", index=False)
    dhw.reset_index().to_csv(OUT_DIR / "reef_effective_dhw.csv", index=False)
    flagged = table[table.model_ahead_bleaching_ge3]
    summary = {
        "status": "tracking_only_scoring_unchanged",
        "convention": "observations scored at LTMP reporting year Y (July Y-1 to June Y); model year Y includes summer-Y bleaching",
        "bleaching_season_end_decimal": BLEACHING_SEASON_END,
        "median_lag_years_by_dataset": table.groupby("dataset").lag_years_vs_model_year.median().round(3).to_dict(),
        "share_surveys_before_model_year_bleaching": float(table.survey_before_model_year_bleaching.mean()),
        "flagged_ge3": flagged[["dataset", "reef_name", "report_year", "survey_date", "model_year_effective_dhw"]]
            .round(3).to_dict(orient="records"),
        "outputs": {name: sha256(OUT_DIR / name) for name in ("survey_timing.csv", "reef_effective_dhw.csv")},
        "inputs": {p.name: sha256(p) for p in (DATA / "reef_cots.csv", DATA / "reef_manta.csv",
                                               DATA / "reef_photo_transect.csv",
                                               DOMAIN / "DHWs" / "dhw_RCPhistorical.nc", Path(__file__))},
        "derived_dir": str(OUT_DIR),
    }
    (run / "survey_timing_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    pd.set_option("display.width", 200)
    print(json.dumps({k: summary[k] for k in ("median_lag_years_by_dataset", "share_surveys_before_model_year_bleaching")}, indent=1))
    print(flagged[["dataset", "reef_name", "report_year", "survey_date", "survey_month", "model_year_effective_dhw",
                   "prev_year_effective_dhw"]].round(2).to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    main(parser.parse_args().run.resolve())
