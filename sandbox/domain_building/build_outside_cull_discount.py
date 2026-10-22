"""Outside-reef cull discount on external COTS supply (crash_mechanism_protocol.md, 2026-10-22,
thirty-third round; lead's agreed "source-reef cull discount").

External supply at local reef i (calendar year y, hydrodynamic year h) is linear in the outbreak
probability of each outside source reef j:
    B_i(y, h) = sum_j C_h[j, i] * A_j / A_i * P_j(y)          (pelagic form: also x survival s_j)
Culling at source j removes a fraction f_j(y) of its adults, so its larval output falls by f_j(y):
    B'_i(y, h) = B_i(y, h) - sum_j C_h[j, i] * A_j / A_i * P_j(y) * f_j(y)

Fraction removed, observed reporting years (CSIRO cull sheet, sheet Cull; outside-domain Owen source reefs):
    CPUE = COTS culled / diver-hours;  D = 7.874 * CPUE^0.549 (SALAD-cull anchor, COTS/ha);
    f = min(F_MAX, COTS culled / (D * A_j)),  A_j = RME habitable area (ha), F_MAX = 0.9.
Reef-years with < 2 diver-hours are not converted (f = 0). Using the whole habitable area makes f
conservative where culls concentrate on outbreak sites.

Scenarios (projection years use each reef's mean f over reporting years 2023-2025, f_recent, which
assumes the recent catch per hour continues; doubled effort doubles f, capped at F_MAX):
    B  observed culls to 2025, none after                 (baseline with discount)
    K  B + carry recent effort on 2026-2030               (f = f_recent)
    D  double effort 2020-2030                            (2020-25: min(F_MAX, 2 f_obs); 2026-30: min(F_MAX, 2 f_recent))
    L  observed to 2025, double recent effort 2026-2030   (2026-30: min(F_MAX, 2 f_recent))
Projected years take P_j from the supply year map (resampled source year) of each wave replicate.

Validation: the undiscounted boundary is reconstructed from the Owen matrices and each dataset's
probability layer and must equal the stored boundary (1991-2100 incl. resampled years) to 1e-9.

    python sandbox/domain_building/build_outside_cull_discount.py [--replicates 5]
Output (licensed inputs, git-ignored): sandbox/data/derived/projection_inputs_v0.5_outside_cull/
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
from datetime import datetime
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_lizard_cots_connectivity_owen as owen  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "sandbox" / "data"
V04 = DATA / "Lizard_Historical_v0.4"
CSIRO = DATA / "260529-COTS-Manta-Cull-RHIS-Lawrence-CSIRO.xlsx"
OUT = DATA / "derived" / "projection_inputs_v0.5_outside_cull"
DATASETS = {
    "owen_global_2018_2023_corr_w3": (DATA / "derived" / "outbreak_layer_correction_w3_x5_2018_2026" / "gbrPredsFaceValue_1991_2026.csv",
                                      DATA / "derived" / "projection_inputs_v0.5" / "supply"),
    "owen_global_2018_2023_corr_w05_graded": (DATA / "derived" / "outbreak_layer_correction_graded_w05_x5_2018_2026" / "gbrPredsFaceValue_1991_2026.csv",
                                              DATA / "derived" / "projection_inputs_v0.5_supply_owen_global_2018_2023_corr_w05_graded" / "supply"),
    "owen_global_2018_2023": ([DATA / "gbrPredsAdj_20262408.csv", DATA / "gbrPreds2026Adj_20262408.csv"],
                              DATA / "derived" / "projection_inputs_v0.5_supply_owen_global_2018_2023" / "supply"),
}
SCENARIOS = {"B": "observed outside culls to 2025, none after",
             "K": "observed to 2025, recent outside effort carried on 2026-2030",
             "D": "outside effort doubled 2020-2030",
             "L": "observed to 2025, recent outside effort doubled 2026-2030"}
F_MAX, MIN_HOURS, RECENT = 0.9, 2.0, (2023, 2025)
D_C, D_K = 7.874, 0.549


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_layer(src) -> pd.DataFrame:
    paths = src if isinstance(src, list) else [src]
    d = pd.concat([pd.read_csv(p)[["year", "reefName", "outbrProb"]] for p in paths], ignore_index=True)
    d["reef_id"] = d.reefName.map(owen.reef_id_from_name)
    return d.pivot(index="reef_id", columns="year", values="outbrProb")


def geometry():
    """Owen global order, local and external masks, areas and survival exactly as the connectivity builder."""
    shape = gpd.read_file(owen.REEF_SHAPE)
    global_ids = [owen.reef_id_from_name(n) for n in shape["reefName"]]
    lookup = {r: i for i, r in enumerate(global_ids)}
    sites = pd.read_csv(V04 / "site_to_reef.csv", dtype=str)
    local_idx = sorted({lookup[r] for r in sites["reef_id"]})
    local_ids = [global_ids[i] for i in local_idx]
    id_list = pd.read_csv(owen.RME_ID_LIST, comment="#", header=None, dtype={0: str})
    area = {str(r.iloc[0]): float(r.iloc[1]) * 100.0 * (1.0 - float(r.iloc[2])) for _, r in id_list.iterrows()}
    survival = pd.read_csv(owen.DEFAULT_SURVIVAL, dtype={"id": str})
    s_lookup = {r: i for i, r in enumerate(survival["id"].astype(str))}
    old = [r if r in area else owen.NEW_TO_OLD_ALIASES.get(r) for r in global_ids]
    mapped = np.array([o is not None and o in s_lookup and o in area for o in old])
    local_mask = np.zeros(len(global_ids), bool)
    local_mask[local_idx] = True
    ext = ~local_mask & mapped
    full_area = np.array([area[o] if m else 0.0 for o, m in zip(old, mapped)])
    surv_years = sorted(int(c) for c in survival.columns if c != "id")
    matrices = {}
    for path in sorted(owen.SOURCE_ROOT.glob("conMatCots*.csv"), key=owen.matrix_year):
        h = owen.matrix_year(path)
        print(f"Reading {path.name}", flush=True)
        m = np.loadtxt(path, delimiter=",", dtype=np.float64)
        m[np.sum(~np.isfinite(m), axis=1) == len(global_ids), :] = 0.0
        sy = h if str(h) in survival.columns else max(surv_years)
        s = np.array([float(survival.loc[s_lookup[o], str(sy)]) if mk else 0.0 for o, mk in zip(old, mapped)])
        local_areas = full_area[local_idx]
        pressure = m[np.ix_(ext, local_idx)] * (full_area[ext, None] / local_areas[None, :])
        matrices[h] = (pressure, pressure * s[ext, None])
    return global_ids, local_ids, np.array(global_ids)[ext].tolist(), full_area[ext], matrices


def cull_fractions(ext_ids: list[str], ext_area: np.ndarray) -> pd.DataFrame:
    c = pd.read_excel(CSIRO, sheet_name="Cull")
    t = pd.to_datetime(c.SurveyDate)
    c["year"] = t.dt.year + (t.dt.month >= 9).astype(int)
    c["reef_id"] = c.ReefName.str.extract(r"\((\d\d-\d{3}[a-z]?)\)\s*$")[0]
    c["hours"] = c.Bottomtime.fillna(0.0).astype(float) / 60.0
    c["cots"] = c[["Cohort1", "Cohort2", "Cohort3", "Cohort4"]].fillna(0).sum(axis=1)
    g = c.groupby(["reef_id", "year"]).agg(hours=("hours", "sum"), cots=("cots", "sum")).reset_index()
    area = dict(zip(ext_ids, ext_area))
    g = g[g.reef_id.isin(area)].copy()
    g["area_ha"] = g.reef_id.map(area)
    g["cpue"] = np.where(g.hours >= MIN_HOURS, g.cots / g.hours.clip(lower=1e-9), np.nan)
    g["density_ha"] = D_C * g.cpue ** D_K
    g["f_obs"] = np.where((g.hours >= MIN_HOURS) & (g.cots > 0),
                          np.minimum(F_MAX, g.cots / (g.density_ha * g.area_ha)), 0.0)
    return g


def scenario_f(g: pd.DataFrame, ext_ids: list[str], years: list[int], code: str) -> np.ndarray:
    idx = {r: i for i, r in enumerate(ext_ids)}
    f = np.zeros((len(ext_ids), len(years)))
    col = {y: k for k, y in enumerate(years)}
    for r in g.itertuples():
        if r.year <= 2025 and r.year in col:
            f[idx[r.reef_id], col[r.year]] = r.f_obs
    recent = g[g.year.between(*RECENT)].groupby("reef_id").f_obs.sum() / (RECENT[1] - RECENT[0] + 1)
    fr = np.zeros(len(ext_ids))
    for rid, v in recent.items():
        fr[idx[rid]] = v
    for y in years:
        k = col[y]
        if code == "D" and 2020 <= y <= 2025:
            f[:, k] = np.minimum(F_MAX, 2 * f[:, k])
        if 2026 <= y <= 2030:
            f[:, k] = {"B": 0.0 * fr, "K": fr, "D": np.minimum(F_MAX, 2 * fr), "L": np.minimum(F_MAX, 2 * fr)}[code]
    return f


def main(replicates: int) -> None:
    if OUT.exists():
        raise FileExistsError(f"Refusing to overwrite {OUT}")
    global_ids, local_ids, ext_ids, ext_area, matrices = geometry()
    cull = cull_fractions(ext_ids, ext_area)
    OUT.mkdir(parents=True)
    cull.to_csv(OUT / "outside_cull_fraction.csv", index=False)
    report = {"status": "derived_input_not_promoted",
              "protocol": "crash_mechanism_protocol.md, 2026-10-22 (thirty-third round): outside-reef cull discount",
              "formula": "B' = B - sum_j C[j,i] A_j/A_i P_j f_j (pelagic: x s_j); f = min(0.9, culled / (7.874 CPUE^0.549 x A_j))",
              "scenarios": SCENARIOS, "f_max": F_MAX, "min_hours": MIN_HOURS, "recent_years": RECENT,
              "first_discount_year": int(cull[cull.f_obs > 0].year.min()),
              "outside_reefs_culled": int(cull[cull.f_obs > 0].reef_id.nunique()),
              "f_obs_summary_by_year": cull[cull.f_obs > 0].groupby("year").f_obs.agg(["count", "mean", "max"]).round(4)
              .to_dict(orient="index"), "datasets": {}}
    for ds, (layer_src, supply_root) in DATASETS.items():
        layer = load_layer(layer_src)
        P = layer.reindex(index=ext_ids).fillna(0.0)
        year_map = pd.read_csv(supply_root / "supply_year_map.csv")
        recon_err, summary = 0.0, {}
        for r in range(1, replicates + 1):
            name = f"wave_r{r:02d}"
            src = supply_root / name / "annual_forcing"
            b = pd.read_csv(src / "external_outbreak_boundary.csv", dtype={"reef_id": str})
            years = sorted(b.calendar_year.unique())
            ym = year_map[year_map.scenario == name].set_index("year").source_year.astype(int).to_dict()
            source_year = [y if y <= 2026 else ym[y] for y in years]
            Pm = P[source_year].to_numpy()                               # ext x years
            local_pos = {rid: k for k, rid in enumerate(local_ids)}
            yi = {y: k for k, y in enumerate(years)}
            # Reconstruct the undiscounted boundary and check it.
            for h, (pot_w, pel_w) in matrices.items():
                rows = b[b.hydrodynamic_year == h]
                li = rows.reef_id.map(local_pos).to_numpy()
                yk = rows.calendar_year.map(yi).to_numpy()
                pot = (pot_w.T @ Pm)[li, yk]
                pel = (pel_w.T @ Pm)[li, yk]
                scale = max(1e-12, float(rows.probability_area_connectivity.abs().max()))
                recon_err = max(recon_err, float(np.max(np.abs(pot - rows.probability_area_connectivity.to_numpy()))) / scale,
                                float(np.max(np.abs(pel - rows.survival_weighted_probability_area_connectivity.to_numpy()))) / scale)
            if recon_err > 1e-9:
                raise ValueError(f"{ds}: reconstructed boundary differs from the stored one ({recon_err:.2e})")
            for code in SCENARIOS:
                f = scenario_f(cull, ext_ids, years, code)
                dst = OUT / ds / code / "supply" / name / "annual_forcing"
                dst.mkdir(parents=True)
                for fpath in src.iterdir():
                    if fpath.name != "external_outbreak_boundary.csv":
                        shutil.copy2(fpath, dst / fpath.name)
                lines = (src / "external_outbreak_boundary.csv").read_text().splitlines()
                out_lines, new = [lines[0]], b.copy()
                for h, (pot_w, pel_w) in matrices.items():
                    sel = (b.hydrodynamic_year == h).to_numpy()
                    rows = b[sel]
                    li = rows.reef_id.map(local_pos).to_numpy()
                    yk = rows.calendar_year.map(yi).to_numpy()
                    dpot = (pot_w.T @ (Pm * f))[li, yk]
                    dpel = (pel_w.T @ (Pm * f))[li, yk]
                    new.loc[sel, "probability_area_connectivity"] = np.maximum(0.0, rows.probability_area_connectivity.to_numpy() - dpot)
                    new.loc[sel, "survival_weighted_probability_area_connectivity"] = np.maximum(
                        0.0, rows.survival_weighted_probability_area_connectivity.to_numpy() - dpel)
                changed = ((new.probability_area_connectivity != b.probability_area_connectivity) |
                           (new.survival_weighted_probability_area_connectivity != b.survival_weighted_probability_area_connectivity)).to_numpy()
                for k, line in enumerate(lines[1:]):
                    if changed[k]:
                        r_ = new.iloc[k]
                        out_lines.append(f"{r_.reef_id},{int(r_.calendar_year)},{int(r_.hydrodynamic_year)},"
                                         f"{repr(float(r_.probability_area_connectivity))},"
                                         f"{repr(float(r_.survival_weighted_probability_area_connectivity))}")
                    else:
                        out_lines.append(line)
                (dst / "external_outbreak_boundary.csv").write_text("\n".join(out_lines) + "\n")
                if r == 1:
                    ratio = (new.groupby("calendar_year").survival_weighted_probability_area_connectivity.sum() /
                             b.groupby("calendar_year").survival_weighted_probability_area_connectivity.sum())
                    summary[code] = {"supply_ratio_by_year_r01": ratio.loc[2012:2030].round(4).to_dict(),
                                     "first_changed_year": int(new.calendar_year[changed].min()) if changed.any() else None}
        report["datasets"][ds] = {"layer": [str(p.relative_to(ROOT)) for p in (layer_src if isinstance(layer_src, list) else [layer_src])],
                                  "reconstruction_max_rel_error": recon_err, "scenario_summary": summary,
                                  "supply_root": str(supply_root.relative_to(ROOT))}
        print(ds, "reconstruction error", recon_err, flush=True)
    report["inputs"] = {"csiro": sha256(CSIRO)} | {
        Path(p).name: sha256(Path(p)) for v in DATASETS.values() for p in (v[0] if isinstance(v[0], list) else [v[0]])}
    report["script_sha256"] = sha256(Path(__file__))
    report["created"] = datetime.now().isoformat(timespec="seconds")
    (OUT / "metadata.json").write_text(json.dumps(report, indent=2, default=str) + "\n")
    print(json.dumps({k: report[k] for k in ("first_discount_year", "outside_reefs_culled")}, indent=1))
    for ds, v in report["datasets"].items():
        for code, s in v["scenario_summary"].items():
            print(ds[-22:], code, {y: s["supply_ratio_by_year_r01"].get(y) for y in (2014, 2020, 2023, 2024, 2025, 2026, 2027, 2029)})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--replicates", type=int, default=5)
    main(parser.parse_args().replicates)
