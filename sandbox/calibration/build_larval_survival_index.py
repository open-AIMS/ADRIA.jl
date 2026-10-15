"""Site x year larval-survival multiplier tables for the opt-in ADRIA switch
COTS_LARVAL_SURVIVAL_CSV (crash_mechanism_protocol.md, "Preregistered chlorophyll
larval-survival switch and 2020s experiment", 2026-10-15).

    multiplier = min(CAP, exp(b * log_anom)),  log_anom = log(wet-season chl-a / reef normal)

from the reef wet-season chl-a index (chla_index_lizard_2003_2026: MODIS 2003-2024, VIIRS
calibrated 2025-2026). Seasons are labelled by spawning year (Dec(year-1)-Apr(year)) and apply
in that model year. Before 2003 and where the index is missing the multiplier is 1 (neutral;
lead's decision until a discharge predictor is tested). Reef values are expanded to the reef's
sites (site_to_reef.csv UNIQUE_ID). Variants: {mean, p90} x b {1.5, 3}.

    python build_larval_survival_index.py [--output <dir>]
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
INDEX = ROOT / "sandbox" / "data" / "derived" / "chla_index_lizard_2003_2026" / "chla_index_reef_year.csv"
SITES = ROOT / "sandbox" / "data" / "Lizard_Historical_v0.4" / "site_to_reef.csv"
DEFAULT_OUT = ROOT / "sandbox" / "data" / "derived" / "larval_survival_index_lizard"
CAP = 5.0
YEARS = range(1991, 2027)
VARIANTS = {"mb15": ("log_anom_mean", 1.5), "mb30": ("log_anom_mean", 3.0),
            "pb15": ("log_anom_p90", 1.5), "pb30": ("log_anom_p90", 3.0)}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(out: Path) -> None:
    if out.exists():
        raise FileExistsError(f"Refusing to overwrite {out}")
    index = pd.read_csv(INDEX, dtype={"RME_UNIQUE_ID": str})
    sites = pd.read_csv(SITES, dtype={"UNIQUE_ID": str, "site_id": str})
    out.mkdir(parents=True)
    summary = {}
    for name, (column, b) in VARIANTS.items():
        reef = index.set_index(["RME_UNIQUE_ID", "year"])[column]
        reef = reef[np.isfinite(reef)]
        grid = pd.MultiIndex.from_product([sites.site_id, list(YEARS)], names=["site_id", "year"]).to_frame(index=False)
        grid = grid.merge(sites[["site_id", "UNIQUE_ID"]], on="site_id")
        anomaly = pd.Series([reef.get((u, y), np.nan) for u, y in zip(grid.UNIQUE_ID, grid.year)])
        grid["multiplier"] = np.where(np.isfinite(anomaly), np.minimum(CAP, np.exp(b * anomaly)), 1.0)
        path = out / f"larval_survival_{name}.csv"
        grid[["site_id", "year", "multiplier"]].to_csv(path, index=False)
        summary[name] = {"column": column, "b": b, "cap": CAP,
                         "domain_median_by_year": grid.groupby("year").multiplier.median().round(3).to_dict(),
                         "share_capped": float((grid.multiplier >= CAP).mean()), "sha256": sha256(path)}
    meta = {"status": "derived_input_not_promoted", "index": str(INDEX.relative_to(ROOT)), "index_sha256": sha256(INDEX),
            "sites": str(SITES.relative_to(ROOT)), "sites_sha256": sha256(SITES),
            "formula": f"min({CAP}, exp(b * log_anom)); 1 before 2003 and where missing", "variants": summary,
            "script_sha256": sha256(Path(__file__))}
    (out / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(pd.DataFrame({k: v["domain_median_by_year"] for k, v in summary.items()}).loc[2009:2026].to_string())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    main(parser.parse_args().output.resolve())
