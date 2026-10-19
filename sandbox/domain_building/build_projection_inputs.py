"""Build the 1991-2100 projection input set v0.5 for the Lizard W-wave-lf04 model
(crash_mechanism_protocol.md, "Preregistered: projection design to 2100", 2026-10-19).

v0.4 (2 GB, mostly connectivity) is not copied. This builder writes only the inputs that
change, and the projection runner swaps them into the loaded v0.4 domain:

  dhw/dhw_<ssp>.nc            scenario (GCM) x location x time 1985-2100; 1985-2026 = v0.4 exactly;
                              2027-2100 = RME CMIP6 annual max DHW, quantile delta mapping per reef
  dhw/dhw_scenarios.csv       ssp, scenario index, gcm, source file, sha256
  dhw/qdm_diagnostics.csv     reference-period fit and warming-signal preservation per GCM x SSP x reef
  env/env_year_map.csv        replicate, year, source_year (2027-2100; iid uniform on 1993-2024)
  env/cyclone_mortality_r<rr>.csv   site x 1985-2100 (v0.4 text kept; 2027+ = source-year mortality)
  env/larval_survival_r<rr>.csv     site_id, year, multiplier (aprmak; 1991-2026 rows kept verbatim)
  mellin/background_mortality.csv   site x 1991-2100 (2027+ = expected loss at 1% simulated incidence)
  supply/<scenario>/annual_forcing/ corr_w3 annual forcing, files copied byte for byte except the
                              external outbreak boundary, extended to 2100 (rows <= 2026 verbatim)
  supply/supply_year_map.csv  scenario, year, source_year
  metadata.json               parameters, source and output hashes, checks

Supply scenarios: wave_r<rr> (circular 15-year block bootstrap of 1991-2026), mean (1991-2026 mean),
closed (zero), repeat (2012-2026 cycled). Units: DHW degC-weeks; cyclone and background mortality are
annual fractions; boundary as in corr_w3 (dimensionless multiplier of outbreak production, COTS ha^-1 yr^-1).

Run with base miniconda Python (pandas, scipy; NetCDF3 written as in build_lizard_domain_v4.py):
    python sandbox/domain_building/build_projection_inputs.py [--replicates 10]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import netcdf_file

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "sandbox" / "data"
V04 = DATA / "Lizard_Historical_v0.4"
RME_DHW = DATA / "rme_ml_2025_06_05" / "data_files" / "dhw_csv"
CHL = DATA / "derived" / "larval_survival_index_lizard_full_1990s_corrected" / "larval_survival_ap_rma_k.csv"
MELLIN = DATA / "derived" / "mellin_background_mortality_1991_2026" / "background_mortality.csv"
MELLIN_COEF = DATA / "derived" / "mellin_disturbance_export" / "mellin_coefficients.csv"
DATASET = "owen_global_2018_2023_corr_w3"
FORCING = V04 / "cots_connectivity" / "datasets" / DATASET / "annual_forcing"
OUT = DATA / "derived" / "projection_inputs_v0.5"

HIST_END, PROJ_END = 2026, 2100
PROJ_YEARS = np.arange(HIST_END + 1, PROJ_END + 1)
REF = (2000, 2026)
WINDOW = 21
ENV_POOL = (1993, 2024)
BLOCK = 15
REPEAT = (2012, 2026)
SEED = 20261019
SIMULATED_INCIDENCE = 0.01
SSP_CODES = {"119": "ssp119", "126": "ssp126", "245": "ssp245", "370": "ssp370", "585": "ssp585"}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_rme_name(name: str) -> tuple[str, str]:
    stem = name.replace("_Dec_2022.csv", "").replace("_annual_DHW_max.csv", "")
    gcm, ssp = stem.rsplit("_", 1)
    return gcm.replace("-", "_"), SSP_CODES[ssp.replace("SSP", "")]


def qdm(obs_ref: np.ndarray, mod_ref: np.ndarray, mod_all: pd.Series, years: np.ndarray) -> np.ndarray:
    """Additive quantile delta mapping (Cannon et al. 2015), centred 21-year window, floored at 0."""
    first, last = int(mod_all.index.min()), int(mod_all.index.max())
    out = np.empty(len(years))
    for i, t in enumerate(years):
        lo = min(max(t - WINDOW // 2, first), last - WINDOW + 1)
        window = mod_all.loc[lo:lo + WINDOW - 1].to_numpy()
        x = float(mod_all.loc[t])
        tau = (np.sum(window < x) + 0.5 * np.sum(window == x)) / len(window)
        tau = min(max(tau, 0.0), 1.0)
        out[i] = max(0.0, np.quantile(obs_ref, tau) + x - np.quantile(mod_ref, tau))
    return out


def build_dhw(site_reef: pd.Series, out: Path, report: dict) -> None:
    with netcdf_file(V04 / "DHWs" / "dhw_RCPhistorical.nc", "r", mmap=False) as hist:
        h = np.array(hist.variables["dhw"][:], dtype=np.float32)  # scenario x location x time
        hist_years = hist.variables["time"][:].astype(int)
        locations = np.array(hist.variables["location"][:], dtype=np.int32)
    if not np.array_equal(h[0], h[1]):
        raise ValueError("v0.4 DHW scenarios differ; expected duplicates")
    if list(hist_years) != list(range(1985, HIST_END + 1)):
        raise ValueError("Unexpected v0.4 DHW years")
    obs = pd.DataFrame(h[0].T.astype(float), index=hist_years, columns=site_reef.values)
    obs_reef = obs.T.groupby(level=0).first().T  # identical within reef (checked below)
    if not np.allclose(obs.T.groupby(level=0).std(ddof=0).fillna(0).to_numpy(), 0.0):
        raise ValueError("Observed DHW varies within a reef; per-reef QDM is invalid")
    files = sorted(RME_DHW.glob("*.csv"))
    by_ssp: dict[str, list[tuple[str, Path]]] = {}
    for f in files:
        gcm, ssp = parse_rme_name(f.name)
        by_ssp.setdefault(ssp, []).append((gcm, f))
    reefs = sorted(set(site_reef))
    rows, diag = [], []
    (out / "dhw").mkdir(parents=True)
    for ssp, items in sorted(by_ssp.items()):
        items.sort()
        cube = np.empty((len(items), len(site_reef), len(hist_years) + len(PROJ_YEARS)), dtype=np.float32)
        for s, (gcm, f) in enumerate(items):
            g = pd.read_csv(f).set_index("id")
            g.columns = g.columns.astype(int)
            missing = sorted(set(reefs) - set(g.index))
            if missing:
                raise ValueError(f"{f.name} lacks reefs {missing[:5]}")
            corrected = {}
            for reef in reefs:
                series = g.loc[reef].astype(float)
                obs_ref = obs_reef.loc[REF[0]:REF[1], reef].to_numpy()
                mod_ref = series.loc[REF[0]:REF[1]].to_numpy()
                corrected[reef] = qdm(obs_ref, mod_ref, series, PROJ_YEARS)
                ref_fit = qdm(obs_ref, mod_ref, series, np.arange(REF[0], REF[1] + 1))
                late = slice(2081, 2100)
                diag.append({"ssp": ssp, "gcm": gcm, "reef_id": reef,
                             "obs_ref_mean": obs_ref.mean(), "raw_ref_mean": mod_ref.mean(),
                             "corr_ref_mean": ref_fit.mean(),
                             "ref_quantile_max_abs_err": float(np.max(np.abs(
                                 np.quantile(ref_fit, np.linspace(0.05, 0.95, 19)) -
                                 np.quantile(obs_ref, np.linspace(0.05, 0.95, 19))))),
                             "raw_change_2081_2100": series.loc[late].mean() - mod_ref.mean(),
                             "corr_change_2081_2100": corrected[reef][-20:].mean() - obs_ref.mean(),
                             "corr_mean_2081_2100": corrected[reef][-20:].mean()})
            cube[s, :, :len(hist_years)] = h[0]
            cube[s, :, len(hist_years):] = np.stack([corrected[r] for r in site_reef.values]).astype(np.float32)
            rows.append({"ssp": ssp, "scenario": s + 1, "gcm": gcm, "source_file": f.name, "sha256": sha256(f)})
        if not (np.all(np.isfinite(cube)) and np.all(cube >= 0.0)):
            raise ValueError(f"Invalid DHW cube for {ssp}")
        with netcdf_file(out / "dhw" / f"dhw_{ssp}.nc", "w") as dataset:
            dataset.createDimension("time", cube.shape[2])
            dataset.createDimension("location", cube.shape[1])
            dataset.createDimension("scenario", cube.shape[0])
            dataset.createVariable("time", "i4", ("time",))[:] = np.concatenate([hist_years, PROJ_YEARS])
            dataset.createVariable("location", "i4", ("location",))[:] = locations
            dataset.createVariable("scenario", "i4", ("scenario",))[:] = np.arange(1, cube.shape[0] + 1)
            variable = dataset.createVariable("dhw", "f4", ("scenario", "location", "time"))
            variable[:] = cube
            variable.units = "degree heating weeks (annual maximum)"
            dataset.history = (f"1985-2026 = Lizard v0.4 dhw_RCPhistorical.nc; 2027-2100 = RME CMIP6 {ssp} annual max DHW, "
                               "additive QDM per reef (ref 2000-2026, 21-yr window), floored at 0")
        with netcdf_file(out / "dhw" / f"dhw_{ssp}.nc", "r", mmap=False) as check:
            written = np.array(check.variables["dhw"][:])
        if not np.array_equal(written[:, :, :len(hist_years)], np.broadcast_to(h[0], written[:, :, :len(hist_years)].shape)):
            raise ValueError("Written DHW does not equal v0.4 before 2027")
    pd.DataFrame(rows).to_csv(out / "dhw" / "dhw_scenarios.csv", index=False)
    d = pd.DataFrame(diag)
    d.to_csv(out / "dhw" / "qdm_diagnostics.csv", index=False)
    summary = d.groupby(["ssp", "gcm"])[["obs_ref_mean", "raw_ref_mean", "corr_ref_mean", "ref_quantile_max_abs_err",
                                         "raw_change_2081_2100", "corr_change_2081_2100", "corr_mean_2081_2100"]].mean()
    report["dhw"] = {"reference": REF, "window_years": WINDOW, "method": "additive QDM, floor 0",
                     "max_ref_quantile_abs_err": float(d.ref_quantile_max_abs_err.max()),
                     "summary": summary.round(3).reset_index().to_dict(orient="records")}


def wide_extend(src: Path, dst: Path, key: str, new_cols: list[str], values: dict[str, list[str]]) -> None:
    """Copy a wide site x year CSV verbatim and append year columns of text tokens (row order kept)."""
    lines = src.read_text().splitlines()
    header = lines[0].split(",")
    if header[0] != key:
        raise ValueError(f"{src} first column is not {key}")
    out = [lines[0] + "," + ",".join(new_cols)]
    for line in lines[1:]:
        site = line.split(",", 1)[0]
        out.append(line + "," + ",".join(values[site]))
    dst.write_text("\n".join(out) + "\n")


def raw_wide(path: Path) -> dict[str, dict[str, str]]:
    """Site -> year column -> the source text token (pandas' float parser is not exactly round-trip)."""
    lines = path.read_text().splitlines()
    header = lines[0].split(",")
    return {parts[0]: dict(zip(header[1:], parts[1:])) for parts in (line.split(",") for line in lines[1:])}


def build_env(sites: list[str], replicates: int, out: Path, report: dict) -> None:
    (out / "env").mkdir(parents=True)
    cyc_src = V04 / "cyclones" / "cyclone_mortality.csv"
    cyc = raw_wide(cyc_src)
    chl_lines = CHL.read_text().splitlines()
    if chl_lines[0] != "site_id,year,multiplier":
        raise ValueError("Unexpected chl table header")
    chl_tok = {}
    for line in chl_lines[1:]:
        site, year, value = line.split(",")
        chl_tok[(site, int(year))] = value
    if max(y for _, y in chl_tok) != HIST_END or str(HIST_END) not in next(iter(cyc.values())):
        raise ValueError("Unexpected historical end year in cyclone or chl tables")
    chl_sites = list(dict.fromkeys(site for site, _ in chl_tok))
    pool = np.arange(ENV_POOL[0], ENV_POOL[1] + 1)
    maps = []
    for r in range(1, replicates + 1):
        rng = np.random.default_rng(SEED + 1000 * r)
        source = rng.choice(pool, size=len(PROJ_YEARS), replace=True)
        maps += [{"replicate": r, "year": int(y), "source_year": int(s)} for y, s in zip(PROJ_YEARS, source)]
        values = {site: [cyc[site][str(s)] for s in source] for site in cyc}
        wide_extend(cyc_src, out / "env" / f"cyclone_mortality_r{r:02d}.csv", "site_id",
                    [str(y) for y in PROJ_YEARS], values)
        extra = [f"{site},{y},{chl_tok[(site, int(s))]}" for site in chl_sites for y, s in zip(PROJ_YEARS, source)]
        (out / "env" / f"larval_survival_r{r:02d}.csv").write_text("\n".join(chl_lines + extra) + "\n")
    pd.DataFrame(maps).to_csv(out / "env" / "env_year_map.csv", index=False)
    if set(cyc) != set(sites):
        raise ValueError("Cyclone table sites differ from the domain")
    report["env"] = {"pool": ENV_POOL, "draw": "iid uniform per projected year, shared by cyclones and chl",
                     "seed_rule": f"{SEED} + 1000 * replicate", "replicates": replicates,
                     "cyclone_source": str(cyc_src.relative_to(ROOT)), "chl_source": str(CHL.relative_to(ROOT))}


def build_mellin(out: Path, report: dict) -> None:
    (out / "mellin").mkdir(parents=True)
    coef = pd.read_csv(MELLIN_COEF).set_index("term").log_cover_effect_mean
    expected = 1.0 - (1.0 - SIMULATED_INCIDENCE * (1.0 - np.exp(coef["disease"]))) * \
        (1.0 - SIMULATED_INCIDENCE * (1.0 - np.exp(coef["unknown"])))
    table = pd.read_csv(MELLIN, dtype={"site_id": str}).set_index("site_id")
    modal = float(table.stack().round(12).mode().iloc[0])
    if abs(modal - expected) > 1e-12:
        raise ValueError(f"Expected-incidence value {expected} is not the table's unobserved value {modal}")
    values = {site: [repr(float(expected))] * len(PROJ_YEARS) for site in table.index}
    wide_extend(MELLIN, out / "mellin" / "background_mortality.csv", "site_id", [str(y) for y in PROJ_YEARS], values)
    report["mellin"] = {"rule": "expected loss at 1% simulated incidence (build_mellin_background_mortality.jl)",
                        "value_2027_2100": expected}


def build_supply(replicates: int, out: Path, report: dict) -> None:
    src = FORCING / "external_outbreak_boundary.csv"
    lines = src.read_text().splitlines()
    b = pd.read_csv(src, dtype={"reef_id": str})
    cols = ["probability_area_connectivity", "survival_weighted_probability_area_connectivity"]
    years = np.arange(int(b.calendar_year.min()), HIST_END + 1)
    if int(b.calendar_year.max()) != HIST_END or len(years) != 36:
        raise ValueError("Unexpected boundary years")
    if lines[0] != "reef_id,calendar_year,hydrodynamic_year," + ",".join(cols):
        raise ValueError("Unexpected boundary header")
    key = {}
    for line in lines[1:]:
        reef, year, hydro, v0, v1 = line.split(",")
        key[(reef, int(year), int(hydro))] = (v0, v1)
    reefs, hydros = sorted(b.reef_id.unique()), sorted(b.hydrodynamic_year.unique())
    scenarios: dict[str, np.ndarray | str] = {}
    for r in range(1, replicates + 1):
        rng = np.random.default_rng(SEED + 2000 * r)
        seq: list[int] = []
        while len(seq) < len(PROJ_YEARS):
            start = int(rng.integers(0, len(years)))
            seq += [int(years[(start + j) % len(years)]) for j in range(BLOCK)]
        scenarios[f"wave_r{r:02d}"] = np.array(seq[:len(PROJ_YEARS)])
    # Mean from exactly parsed tokens (Python float() is correctly rounded; pandas' default parser is not).
    mean = {(reef, h): tuple(float(np.mean([float(key[(reef, int(y), int(h))][c]) for y in years])) for c in (0, 1))
            for reef in reefs for h in hydros}
    rep = np.arange(REPEAT[0], REPEAT[1] + 1)
    scenarios["repeat"] = np.array([rep[i % len(rep)] for i in range(len(PROJ_YEARS))])
    scenarios["mean"], scenarios["closed"] = "mean", "closed"
    year_map = []
    for name, source in scenarios.items():
        d = out / "supply" / name / "annual_forcing"
        d.mkdir(parents=True)
        for f in FORCING.iterdir():
            if f.name not in ("external_outbreak_boundary.csv",):
                shutil.copy2(f, d / f.name)
        extra = []
        for i, y in enumerate(PROJ_YEARS):
            for reef in reefs:
                for h in hydros:
                    if isinstance(source, str):
                        v = mean[(reef, h)] if source == "mean" else (0.0, 0.0)
                        v = (repr(float(v[0])), repr(float(v[1])))
                    else:
                        v = key[(reef, int(source[i]), int(h))]
                    extra.append(f"{reef},{y},{h},{v[0]},{v[1]}")
            year_map.append({"scenario": name, "year": int(y),
                             "source_year": int(source[i]) if not isinstance(source, str) else source})
        (d / "external_outbreak_boundary.csv").write_text("\n".join(lines + extra) + "\n")
    pd.DataFrame(year_map).to_csv(out / "supply" / "supply_year_map.csv", index=False)
    report["supply"] = {"dataset": DATASET, "scenarios": list(scenarios), "block_years": BLOCK,
                        "pool": [int(years[0]), int(years[-1])], "repeat": REPEAT,
                        "seed_rule": f"{SEED} + 2000 * replicate",
                        "copied_files_sha256": {f.name: sha256(f) for f in FORCING.iterdir()}}


def main(replicates: int) -> None:
    if OUT.exists():
        raise FileExistsError(f"Refusing to overwrite {OUT}")
    site_map = pd.read_csv(V04 / "site_to_reef.csv", dtype={"site_id": str, "reef_id": str})
    if list(site_map.indexV2) != list(range(1, len(site_map) + 1)):
        raise ValueError("site_to_reef.csv is not in domain location order")
    site_reef = site_map.set_index("site_id").reef_id
    OUT.mkdir(parents=True)
    report: dict = {"status": "projection_inputs_not_promoted",
                    "protocol": "crash_mechanism_protocol.md, Preregistered: projection design to 2100 (2026-10-19)",
                    "replay": "every table equals the v0.4 / W-wave-lf04 input in 1991-2026 (text kept verbatim, DHW array equal)",
                    "years": [1985, PROJ_END], "projected_years": [int(PROJ_YEARS[0]), int(PROJ_YEARS[-1])]}
    build_dhw(site_reef, OUT, report)
    build_env(list(site_reef.index), replicates, OUT, report)
    build_mellin(OUT, report)
    build_supply(replicates, OUT, report)
    report["inputs"] = {str(p.relative_to(ROOT)): sha256(p) for p in (
        V04 / "DHWs" / "dhw_RCPhistorical.nc", V04 / "cyclones" / "cyclone_mortality.csv", V04 / "site_to_reef.csv",
        CHL, MELLIN, MELLIN_COEF, FORCING / "external_outbreak_boundary.csv")}
    report["inputs"]["rme_dhw_csv"] = "per file in dhw/dhw_scenarios.csv"
    report["script_sha256"] = sha256(Path(__file__))
    report["created"] = datetime.now().isoformat(timespec="seconds")
    report["outputs"] = {str(p.relative_to(OUT)).replace("\\", "/"): sha256(p)
                         for p in sorted(OUT.rglob("*")) if p.is_file()}
    (OUT / "metadata.json").write_text(json.dumps(report, indent=2, default=str) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k not in ("outputs", "dhw")}, indent=1, default=str)[:3000])
    print(pd.DataFrame(report["dhw"]["summary"]).to_string())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--replicates", type=int, default=10)
    main(parser.parse_args().replicates)
