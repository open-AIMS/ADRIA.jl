"""Full 1993-2026 reef wet-season chlorophyll-a series on the in-water scale (lead's decisions,
2026-10-15, thirteenth round; crash_mechanism_protocol.md, "Preregistered full chlorophyll
series (in-water scale) and Ext re-test").

- 2003-2024: MODIS-Aqua OC3 reef season mean and p90 (extract_chla_wet_season.py) converted to
  in-water with the same-day matchup fit log10(iw) = a + b log10(OC3)
  (build_chla_inwater_matchups.py). Primary: OLS (preregistered). Sensitivity: reduced major
  axis (b_rma = b / r), because daily pixel noise attenuates the OLS slope.
- 2025-2026: VIIRS season values mapped to the MODIS scale per reef as
  MODIS = MODIS reef climatology x exp(s x VIIRS log anomaly), s = the pooled VIIRS-to-MODIS
  anomaly slope (extract_chla_viirs_wet_season.py), then converted.
- 1993-2002: Lizard-region (within 0.6 deg of Lizard Island) in-water season mean and p90 over
  all samples (all depths; seasons with >= 16 samples), times each reef's mean ratio to the
  domain median over 2003-2026 (converted satellite). 1991-1992: missing (neutral downstream).
Anomalies: log(value / reef geometric mean over 1993-2026), for mean and p90.

    python build_chla_index_full.py [--output <dir>]
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DERIVED = ROOT / "sandbox" / "data" / "derived"
MODIS = DERIVED / "chla_wet_season_lizard" / "chla_wet_season_reef_year.csv"
VIIRS_DIR = DERIVED / "chla_index_lizard_2003_2026"
MATCH = DERIVED / "chla_inwater_matchups"
INWATER = ROOT / "sandbox" / "data" / "GBR_chlorophyll_monitoring_1992-2009.csv"
DEFAULT_OUT = DERIVED / "chla_index_lizard_1993_2026_inwater"
LIZARD = (-14.68, 145.45)
MIN_SAMPLES = 16


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(out: Path) -> None:
    if out.exists():
        raise FileExistsError(f"Refusing to overwrite {out}")
    meta_match = json.loads((MATCH / "metadata.json").read_text())
    ols = meta_match["fits"]["all_depths"]
    pairs = pd.read_csv(MATCH / "matchups.csv")
    lx, ly = np.log10(pairs.oc3), np.log10(pairs.inwater_all)
    b_rma = ols["b"] / ols["r"]
    fits = {"ols": (ols["a"], ols["b"]), "rma": (float(ly.mean() - b_rma * lx.mean()), float(b_rma))}

    modis = pd.read_csv(MODIS, dtype={"RME_UNIQUE_ID": str})
    viirs = pd.read_csv(VIIRS_DIR / "viirs_wet_season_reef_year.csv", dtype={"RME_UNIQUE_ID": str})
    slope = json.loads((VIIRS_DIR / "metadata.json").read_text())["calibration"]
    clim = modis.groupby("RME_UNIQUE_ID")[["mean", "p90"]].agg(lambda v: float(np.exp(np.nanmean(np.log(v)))))
    v = viirs[viirs.year >= 2025].join(clim.rename(columns={"mean": "cm", "p90": "cp"}), on="RME_UNIQUE_ID")
    v["mean"] = v.cm * np.exp(slope["mean"]["slope_no_intercept"] * v.log_anom_mean)
    v["p90"] = v.cp * np.exp(slope["p90"]["slope_no_intercept"] * v.log_anom_p90)
    # MODIS ends 2024-12-25: its 2025-2026 rows are empty and must not shadow the VIIRS rows.
    sat = pd.concat([modis[modis.year <= 2024][["RME_UNIQUE_ID", "reef_id", "year", "mean", "p90"]].assign(source="MODIS"),
                     v[["RME_UNIQUE_ID", "year", "mean", "p90"]].assign(source="VIIRS")], ignore_index=True)
    sat["reef_id"] = sat.RME_UNIQUE_ID.map(modis.drop_duplicates("RME_UNIQUE_ID").set_index("RME_UNIQUE_ID").reef_id)

    d = pd.read_csv(INWATER, skiprows=26, low_memory=False)
    d["date"] = pd.to_datetime(d.SAMPLE_TIME, dayfirst=True, errors="coerce")
    d["CHL_A"] = pd.to_numeric(d.CHL_A, errors="coerce")
    dist = np.hypot(d.LATITUDE - LIZARD[0], (d.LONGITUDE - LIZARD[1]) * np.cos(np.radians(-LIZARD[0])))
    w = d[(dist < 0.6) & (d.CHL_A > 0) & d.date.dt.month.isin([12, 1, 2, 3, 4])].copy()
    w["year"] = w.date.dt.year + (w.date.dt.month == 12).astype(int)
    region = w.groupby("year").CHL_A.agg(n="size", mean="mean", p90=lambda x: x.quantile(0.9))
    region = region[(region.n >= MIN_SAMPLES) & region.index.isin(range(1993, 2003))]

    out.mkdir(parents=True)
    summary = {}
    frames = []
    for name, (a, b) in fits.items():
        s = sat.copy()
        s["mean_iw"] = 10 ** a * s["mean"] ** b
        s["p90_iw"] = 10 ** a * s["p90"] ** b
        dom = s.groupby("year")[["mean_iw", "p90_iw"]].median()
        ratio = s.join(dom.rename(columns={"mean_iw": "dm", "p90_iw": "dp"}), on="year")
        offs = ratio.assign(om=np.log(ratio.mean_iw / ratio.dm), op=np.log(ratio.p90_iw / ratio.dp)) \
                    .groupby("RME_UNIQUE_ID")[["om", "op"]].mean()
        early = pd.MultiIndex.from_product([offs.index, region.index], names=["RME_UNIQUE_ID", "year"]).to_frame(index=False)
        early = early.join(offs, on="RME_UNIQUE_ID").join(region, on="year")
        early["mean_iw"] = early["mean"] * np.exp(early.om)
        early["p90_iw"] = early["p90"] * np.exp(early.op)
        early["source"] = "in-water regional x reef offset"
        early["reef_id"] = early.RME_UNIQUE_ID.map(sat.drop_duplicates("RME_UNIQUE_ID").set_index("RME_UNIQUE_ID").reef_id)
        full = pd.concat([early[["RME_UNIQUE_ID", "reef_id", "year", "mean_iw", "p90_iw", "source"]],
                          s[["RME_UNIQUE_ID", "reef_id", "year", "mean_iw", "p90_iw", "source"]]], ignore_index=True)
        gm = full.groupby("RME_UNIQUE_ID")[["mean_iw", "p90_iw"]].agg(lambda x: float(np.exp(np.nanmean(np.log(x)))))
        full = full.join(gm.rename(columns={"mean_iw": "gm_mean", "p90_iw": "gm_p90"}), on="RME_UNIQUE_ID")
        full["log_anom_mean"] = np.log(full.mean_iw / full.gm_mean)
        full["log_anom_p90"] = np.log(full.p90_iw / full.gm_p90)
        full["conversion"] = name
        frames.append(full)
        med = full.groupby("year")[["mean_iw", "p90_iw", "log_anom_p90"]].median()
        summary[name] = {"a": a, "b": b, "domain_median_by_year": json.loads(med.round(4).to_json(orient="index"))}
    table = pd.concat(frames, ignore_index=True).sort_values(["conversion", "RME_UNIQUE_ID", "year"])
    if table.duplicated(["conversion", "RME_UNIQUE_ID", "year"]).any():
        raise ValueError("Duplicate reef-years in the chl series")
    table.to_csv(out / "chla_index_reef_year.csv", index=False)
    region.to_csv(out / "inwater_region_wet_season.csv")
    meta = {"status": "derived_input_not_promoted", "fits": summary, "region_samples": region.n.to_dict(),
            "inputs": {"modis_sha256": sha256(MODIS), "viirs_sha256": sha256(VIIRS_DIR / "viirs_wet_season_reef_year.csv"),
                       "matchups_sha256": sha256(MATCH / "matchups.csv"), "inwater_sha256": sha256(INWATER)},
            "script_sha256": sha256(Path(__file__)),
            "outputs": {"chla_index_reef_year.csv": sha256(out / "chla_index_reef_year.csv")}}
    (out / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
    show = table.groupby(["year", "conversion"])[["p90_iw", "log_anom_p90"]].median().unstack("conversion")
    print(show.round(3).to_string())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    main(parser.parse_args().output.resolve())
