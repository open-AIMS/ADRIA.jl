"""Build the opt-in wave-exposure settlement modifier (COTS_SETTLEMENT_MODIFIER_CSV;
crash_mechanism_protocol.md, "Preregistered: wave-exposure settlement mechanism", 2026-10-19).

Input: sandbox/data/derived/never_cots_diagnostic/site_waves.csv (analyse_never_cots_reefs.py; site
means of the GBR ReefGuidance processed wave rasters, Hs in m, peak period Tp in s, power index
Hs^2 Tp). Sites without raster values (29, including all of Linnet) take the values of the nearest
site with values (haversine on site centroids); the fill is recorded per site.

Exposure E (per site), by --metric:
- tp_reef: area-weighted mean Tp of the site's reef (s), shared by all its sites;
- power_reef_max: the most exposed site's power index on the reef, shared by all its sites;
- power_site: the site's own power index.
Multiplier (dimensionless, in [floor, 1]): linear ramp
    m = 1 - (1 - floor) * clamp((E - e0) / (e1 - e0), 0, 1)
so sites at or below e0 are unchanged and sites at or above e1 settle at `floor`.

    python build_wave_settlement_modifier.py --metric tp_reef --e0 6.0 --e1 7.0 --floor 0.25 --name tp60_70_f025

Output: sandbox/data/derived/wave_settlement_modifier/<name>/{settlement_modifier.csv, reefs.csv,
metadata.json} (git-ignored; never overwritten).
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "sandbox" / "data"
WAVES = DATA / "derived" / "never_cots_diagnostic" / "site_waves.csv"
SITE_MAP = DATA / "Lizard_Historical_v0.4" / "site_to_reef.csv"
OUT_ROOT = DATA / "derived" / "wave_settlement_modifier"
METRICS = ("tp_reef", "power_reef_max", "power_site")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def haversine_km(lat0, lon0, lat, lon):
    p0, p = np.radians(lat0), np.radians(lat)
    a = np.sin((p - p0) / 2) ** 2 + np.cos(p0) * np.cos(p) * np.sin(np.radians(lon - lon0) / 2) ** 2
    return 6371.0 * 2 * np.arcsin(np.sqrt(a))


def filled_waves() -> pd.DataFrame:
    w = pd.read_csv(WAVES)
    if w.site_id.duplicated().any():
        raise ValueError("Duplicate sites in site_waves.csv")
    valid = w.dropna(subset=["hs_mean", "tp_mean", "power_mean"]).reset_index(drop=True)
    w["filled_from"] = ""
    w["fill_distance_km"] = 0.0
    for i in w.index[w.tp_mean.isna() | w.power_mean.isna() | w.hs_mean.isna()]:
        d = haversine_km(w.at[i, "y_coord"], w.at[i, "x_coord"], valid.y_coord.to_numpy(), valid.x_coord.to_numpy())
        k = int(np.argmin(d))
        for c in ("hs_mean", "hs_p90", "tp_mean", "power_mean"):
            w.at[i, c] = valid.at[k, c]
        w.at[i, "filled_from"] = valid.at[k, "site_id"]
        w.at[i, "fill_distance_km"] = float(d[k])
    return w


def main(metric: str, e0: float, e1: float, floor: float, name: str) -> None:
    if metric not in METRICS:
        raise ValueError(f"metric must be one of {METRICS}")
    if not (e1 > e0 and 0.0 <= floor <= 1.0):
        raise ValueError("Need e1 > e0 and floor in [0, 1]")
    out = OUT_ROOT / name
    if out.exists():
        raise FileExistsError(f"Refusing to overwrite {out}")
    w = filled_waves()
    sites = pd.read_csv(SITE_MAP)
    if set(sites.site_id) != set(w.site_id) or len(sites) != len(w):
        raise ValueError("site_waves.csv and the v0.4 site map have different sites")
    reef_tp = w.groupby("GBRMPA_ID").apply(lambda g: np.average(g.tp_mean, weights=g.area), include_groups=False)
    reef_pmax = w.groupby("GBRMPA_ID").power_mean.max()
    w["exposure"] = (w.GBRMPA_ID.map(reef_tp) if metric == "tp_reef" else
                     w.GBRMPA_ID.map(reef_pmax) if metric == "power_reef_max" else w.power_mean)
    # Clipped to [floor, 1] so floating-point rounding (e.g. 1 - 0.9 < 0.1) cannot leave the range.
    w["multiplier"] = np.clip(1.0 - (1.0 - floor) * np.clip((w.exposure - e0) / (e1 - e0), 0.0, 1.0), floor, 1.0)
    if not (np.isfinite(w.multiplier).all() and w.multiplier.between(floor, 1.0).all()):
        raise ValueError("Invalid multipliers")
    out.mkdir(parents=True)
    w[["site_id", "multiplier"]].to_csv(out / "settlement_modifier.csv", index=False)
    reefs = w.groupby(["GBRMPA_ID", "reef_name"]).apply(lambda g: pd.Series({
        "sites": len(g), "filled_sites": int((g.filled_from != "").sum()),
        "tp_reef_mean_s": np.average(g.tp_mean, weights=g.area), "power_reef_max": g.power_mean.max(),
        "exposure_area_mean": np.average(g.exposure, weights=g.area),
        "multiplier_area_mean": np.average(g.multiplier, weights=g.area),
        "multiplier_min": g.multiplier.min()}), include_groups=False).reset_index()
    reefs.to_csv(out / "reefs.csv", index=False)
    fills = w[w.filled_from != ""][["site_id", "reef_name", "filled_from", "fill_distance_km"]]
    meta = {"mechanism": "wave-exposure settlement modifier (COTS_SETTLEMENT_MODIFIER_CSV)",
            "metric": metric, "e0": e0, "e1": e1, "floor": floor,
            "units": {"tp_reef": "s", "power_reef_max": "m^2 s (Hs^2 Tp index)", "power_site": "m^2 s (Hs^2 Tp index)"}[metric],
            "form": "m = 1 - (1 - floor) * clamp((E - e0) / (e1 - e0), 0, 1); time-invariant; per site",
            "inputs": {"site_waves.csv": sha256(WAVES), "site_to_reef.csv": sha256(SITE_MAP)},
            "raster_source": "GBR-ReefGuidance-processed_2025-11-19 {Cairns-Cooktown,FarNorthern}_waves_{Hs,Tp}.tif (via analyse_never_cots_reefs.py)",
            "gap_fill": {"rule": "nearest site with raster values (haversine)", "sites": len(fills),
                         "reefs": sorted(set(fills.reef_name)), "max_distance_km": float(fills.fill_distance_km.max()) if len(fills) else 0.0,
                         "per_site": fills.to_dict(orient="records")},
            "dimensions": {"sites": len(w), "reefs": len(reefs)},
            "summary": {"sites_below_1": int((w.multiplier < 1).sum()), "sites_at_floor": int((w.multiplier == floor).sum()),
                        "reefs_mean_below_0.99": int((reefs.multiplier_area_mean < 0.99).sum())},
            "outputs": {n: sha256(out / n) for n in ("settlement_modifier.csv", "reefs.csv")},
            "script_sha256": sha256(Path(__file__))}
    (out / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps({k: meta[k] for k in ("metric", "e0", "e1", "floor", "summary")}, indent=1))
    print(json.dumps({k: v for k, v in meta["gap_fill"].items() if k != "per_site"}, indent=1))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metric", required=True, choices=METRICS)
    parser.add_argument("--e0", type=float, required=True)
    parser.add_argument("--e1", type=float, required=True)
    parser.add_argument("--floor", type=float, required=True)
    parser.add_argument("--name", required=True)
    a = parser.parse_args()
    main(a.metric, a.e0, a.e1, a.floor, a.name)
