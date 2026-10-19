"""Why is COTS self-recruitment so low in the model? (lead's question, 2026-10-15, sixteenth round).

In ADRIA/COTSMod (separated recruitment) a reef's own-reef settlement is
    own production per ha x C[i, i] x settlement,
where C is the Owen hydrodynamic connectivity (rows sources, columns sinks; the fraction of a
reef's larvae that settle on each reef). No parameter scales the diagonal (immigration scalar,
settlement probability = 1). This script reports, per calibration reef and averaged over the
five hydrodynamic years:
- self_retention: C[i, i], the fraction of the reef's larvae that return to it;
- settle_anywhere: sum_j C[i, j], the fraction settling on any GBR reef;
- to_domain: the fraction settling on the 113 domain reefs (including itself);
- sink shares under equal per-ha production everywhere: of larvae arriving at i, the share from
  i itself, from other domain reefs, and from non-domain reefs (sum_j C[j, i] A_j / A_i);
- the same with the model's actual production asymmetry: non-domain outbreak reefs at
  P = 1 produce `boundary x k x Ext` per ha versus local `fecundity x adults` (both from run
  metadata), shown as the adult density a domain reef would need to match one outbreak reef.
Lizard Island Reef is three domain reefs (14-116b, c, d): settlement among them is "internal",
so its own-reef retention is also reported for the three combined.

    python analyse_self_recruitment.py [--output <dir>]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import tomllib
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "sandbox" / "data"
OWEN = DATA / "OwenHydro" / "connectivityMatsYearlyGlobal"
DOMAIN = DATA / "Lizard_Historical_v0.4"
RME = DATA / "rme_ml_2025_06_05" / "data_files"
RUN_META = ROOT / "sandbox/calibration/runs/20261015T_36_init_cull_original_1991_v1/metadata.toml"
CAL = {"Lizard Island Reef": ["14-116b", "14-116c", "14-116d"], "MacGillivray Reef": ["14-114"],
       "North Direction Reef": ["14-143"], "Eyrie Reef": ["14-118"]}
ALIASES = {"11-325": "10-441", "11-244e": "11-288", "11-244f": "11-303", "11-244g": "11-310",
           "11-244h": "11-311", "20198": "20-198"}
DEFAULT_OUT = DATA / "derived" / "self_recruitment_diagnostic"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main(out: Path) -> None:
    if out.exists():
        raise FileExistsError(f"Refusing to overwrite {out}")
    import geopandas as gpd
    shape = gpd.read_file(OWEN / "reefShapefiles" / "gbrShapeLL.shp")
    ids = [re.search(r"\(([^()]*)\)\s*$", str(n)).group(1).strip() for n in shape.reefName]
    index = {r: k for k, r in enumerate(ids)}
    local = sorted({index[r] for r in pd.read_csv(DOMAIN / "site_to_reef.csv", dtype=str).reef_id})
    in_domain = np.zeros(len(ids), bool)
    in_domain[local] = True
    idl = pd.read_csv(RME / "id" / "id_list_2024_12_01.csv", comment="#", header=None, dtype={0: str})
    area_old = {str(r.iloc[0]): float(r.iloc[1]) * 100.0 * (1.0 - float(r.iloc[2])) for _, r in idl.iterrows()}
    area = np.array([area_old.get(r, area_old.get(ALIASES.get(r, ""), 0.0)) for r in ids])
    mats = sorted(OWEN.glob("conMatCots*.csv"))
    rows = []
    for path in mats:
        c = np.loadtxt(path, delimiter=",", dtype=np.float64)
        c[~np.isfinite(c).all(axis=1)] = 0.0
        for cal, members in CAL.items():
            k = [index[r] for r in members]
            a = area[k]
            # Larvae from the group (area-weighted production) and arrivals at the group.
            out_frac = (a[:, None] * c[k, :]).sum(axis=0) / a.sum()          # per-larva fate of group larvae
            arrivals = (area[:, None] * c[:, k]).sum(axis=1)                 # equal per-ha production everywhere
            single = [float(c[j, j]) for j in k]
            rows.append({"calibration_reef": cal, "hydro": path.name,
                         "self_retention_single_reef_mean": float(np.average(single, weights=a)),
                         "self_retention_group": float(out_frac[k].sum()),
                         "settle_anywhere": float(out_frac.sum()), "to_domain": float(out_frac[in_domain].sum()),
                         "sink_share_self_group": float(arrivals[k].sum() / arrivals.sum()),
                         "sink_share_self_single_reef": float(np.average([area[j] * c[j, j] / (area * c[:, j]).sum() for j in k], weights=a)),
                         "sink_share_other_domain": float(arrivals[in_domain].sum() - arrivals[k].sum()) / arrivals.sum(),
                         "sink_share_external": float(arrivals[~in_domain].sum() / arrivals.sum())})
    d = pd.DataFrame(rows)
    summary = d.drop(columns="hydro").groupby("calibration_reef").mean()
    meta_run = tomllib.loads(RUN_META.read_text())
    fec = float(meta_run["fecundity_pre_pelagic_recruits_adult_year"])
    ext_prod = float(meta_run["unit_migration"]["migrated"]["external_production_pre_pelagic_recruits_ha_year"])
    production = {"local_fecundity_per_adult_year": fec, "external_outbreak_production_ha_year_ext1": ext_prod,
                  "adults_ha_equivalent_of_one_outbreak_reef_ext1": ext_prod / fec,
                  "adults_ha_equivalent_of_one_outbreak_reef_ext5": 5 * ext_prod / fec,
                  "note": "local production is further reduced by body condition^2, the gradual Allee term "
                          "exp(-(3/A)^0.5) and Ricker exp(-b A); external is multiplied by P (0-1)"}
    out.mkdir(parents=True)
    d.to_csv(out / "self_recruitment_by_hydro_year.csv", index=False)
    summary.to_csv(out / "self_recruitment_summary.csv")
    meta = {"status": "diagnostic", "production": production,
            "inputs": {str(p.relative_to(ROOT)): sha256(p) for p in [*mats, DOMAIN / "site_to_reef.csv", RUN_META]},
            "script_sha256": sha256(Path(__file__)),
            "outputs": {n: sha256(out / n) for n in ("self_recruitment_by_hydro_year.csv", "self_recruitment_summary.csv")}}
    (out / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
    pd.set_option("display.width", 250)
    print(summary.T.to_string(float_format=lambda v: f"{v:.4g}"))
    print(json.dumps(production, indent=1))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    main(parser.parse_args().output.resolve())
