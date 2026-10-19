"""Year-specific chl survival multiplier for external pelagic supply (COTS_EXTERNAL_SURVIVAL_CSV;
crash_mechanism_protocol.md, "Preregistered: chl survival on external supply", 2026-10-19,
twenty-fourth round).

External source reefs lie outside the domain and have no chl index of their own, so the regional
proxy is the area-weighted mean over domain sites of the `aprmak` larval-survival table
(build_larval_survival_index_corrected.py; multiplier = min(5, exp(3 x p90 log anomaly)), 1990s
corrected by k): m(y) = sum_s area_s x multiplier(s, y) / sum_s area_s. Dimensionless; years
1991-1992 are 1 (no index).

    python build_external_chl_survival.py [--output <dir>]
Output: sandbox/data/derived/external_chl_survival_aprmak/{external_survival.csv, metadata.json}
(git-ignored; never overwritten).
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "sandbox" / "data"
TABLE = DATA / "derived" / "larval_survival_index_lizard_full_1990s_corrected" / "larval_survival_ap_rma_k.csv"
SITES = DATA / "Lizard_Historical_v0.4" / "site_to_reef.csv"
DEFAULT_OUT = DATA / "derived" / "external_chl_survival_aprmak"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(out: Path) -> None:
    if out.exists():
        raise FileExistsError(f"Refusing to overwrite {out}")
    t = pd.read_csv(TABLE)
    sites = pd.read_csv(SITES)[["site_id", "area_m2"]]
    if set(t.site_id) != set(sites.site_id):
        raise ValueError("Survival table and site map have different sites")
    t = t.merge(sites, on="site_id")
    if t.groupby("year").site_id.nunique().min() != len(sites):
        raise ValueError("Survival table is missing site-years")
    m = (t.assign(w=t.multiplier * t.area_m2).groupby("year").w.sum() / t.groupby("year").area_m2.sum()).rename("multiplier")
    if not (m.notna().all() and (m >= 0).all()):
        raise ValueError("Invalid multipliers")
    out.mkdir(parents=True)
    m.reset_index().to_csv(out / "external_survival.csv", index=False)
    meta = {"mechanism": "chl survival on external pelagic supply (COTS_EXTERNAL_SURVIVAL_CSV)",
            "form": "area-weighted mean over domain sites of the aprmak multiplier, per calendar year",
            "units": "dimensionless", "years": [int(m.index.min()), int(m.index.max())],
            "mean_1993_2026": float(m.loc[1993:].mean()),
            "values": {str(int(y)): round(float(v), 4) for y, v in m.items()},
            "inputs": {"larval_survival_ap_rma_k.csv": sha256(TABLE), "site_to_reef.csv": sha256(SITES)},
            "output": sha256(out / "external_survival.csv"), "script_sha256": sha256(Path(__file__))}
    (out / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps({k: meta[k] for k in ("years", "mean_1993_2026")}), m.round(3).to_dict())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    main(parser.parse_args().output.resolve())
