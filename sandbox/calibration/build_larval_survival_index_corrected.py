"""1990s chlorophyll peak bias correction for the p90-anomaly larval-survival table (lead's
decisions, 2026-10-15, fourteenth round; crash_mechanism_protocol.md, "Lead's decisions
(2026-10-15, fourteenth round) and next steps", step 1).

The 1993-2002 segment of the full series (build_chla_index_full.py) is the Lizard-region
in-water p90 spliced onto the satellite era. No high-chl season falls in the 2003-2008
overlap, so the splice cannot calibrate the peaks, and the 1990s domain-median p90 anomaly
varies ~4.5x more than the satellite era. Correction (RMA conversion, p90 anomaly only):

    k = mean(domain-median anomaly 2011, 2024) / mean(domain-median anomaly 1996, 1999)
    log_anom_p90 <- k * log_anom_p90   for 1993-2002 (both signs; reef long-term mean kept)

which brings the 1996 and 1999 peaks in line with the 2011 and 2024 satellite peaks.
Satellite-era rows (2003-2026) are unchanged and must reproduce the uncorrected ap_rma
table exactly (checked). Multiplier = min(5, exp(3 log_anom_p90)); 1 where missing (1991-92).

    python build_larval_survival_index_corrected.py [--output <dir>] [--figure <png>]
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

ROOT = Path(__file__).resolve().parents[2]
DERIVED = ROOT / "sandbox" / "data" / "derived"
INDEX = DERIVED / "chla_index_lizard_1993_2026_inwater" / "chla_index_reef_year.csv"
UNCORRECTED = DERIVED / "larval_survival_index_lizard_full" / "larval_survival_ap_rma.csv"
SITES = ROOT / "sandbox" / "data" / "Lizard_Historical_v0.4" / "site_to_reef.csv"
DEFAULT_OUT = DERIVED / "larval_survival_index_lizard_full_1990s_corrected"
DEFAULT_FIG = ROOT / "sandbox" / "calibration" / "figures" / "larval_survival_1990s_correction.png"
CAP, B = 5.0, 3.0
INWATER_YEARS = (1993, 2002)
PEAKS_1990S, PEAKS_SATELLITE = (1996, 1999), (2011, 2024)
YEARS = range(1991, 2027)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(out: Path, figure: Path) -> None:
    if out.exists():
        raise FileExistsError(f"Refusing to overwrite {out}")
    idx = pd.read_csv(INDEX, dtype={"RME_UNIQUE_ID": str})
    sites = pd.read_csv(SITES, dtype={"UNIQUE_ID": str, "site_id": str})
    rma = idx[idx.conversion == "rma"].dropna(subset=["log_anom_p90"]).copy()
    if rma.duplicated(["RME_UNIQUE_ID", "year"]).any():
        raise ValueError("Duplicate reef-years in the chl series")
    med = rma.groupby("year").log_anom_p90.median()
    k = float(med.loc[list(PEAKS_SATELLITE)].mean() / med.loc[list(PEAKS_1990S)].mean())
    early = rma.year.between(*INWATER_YEARS)
    if set(rma[early].source) != {"in-water regional x reef offset"} or (rma[~early].source.str.contains("in-water")).any():
        raise ValueError("The in-water segment is not exactly 1993-2002")
    rma["log_anom_p90_uncorrected"] = rma.log_anom_p90
    rma.loc[early, "log_anom_p90"] = k * rma.loc[early, "log_anom_p90"]
    rma["p90_iw_corrected"] = rma.gm_p90 * np.exp(rma.log_anom_p90)
    rma["k"] = np.where(early, k, 1.0)

    values = rma.set_index(["RME_UNIQUE_ID", "year"]).log_anom_p90.map(lambda a: float(np.minimum(CAP, np.exp(B * a))))
    grid = pd.MultiIndex.from_product([sites.site_id, list(YEARS)], names=["site_id", "year"]).to_frame(index=False)
    grid = grid.merge(sites[["site_id", "UNIQUE_ID"]], on="site_id")
    m = pd.Series([values.get((u, y), np.nan) for u, y in zip(grid.UNIQUE_ID, grid.year)], dtype=float)
    table = grid.assign(multiplier=np.where(np.isfinite(m), m, 1.0))[["site_id", "year", "multiplier"]]

    # Control: satellite-era and missing years reproduce the uncorrected ap_rma table exactly.
    old = pd.read_csv(UNCORRECTED, dtype={"site_id": str}, float_precision="round_trip")
    both = table.merge(old, on=["site_id", "year"], suffixes=("", "_old"), validate="one_to_one")
    if len(both) != len(table) or len(old) != len(table):
        raise ValueError("Site x year grid differs from the uncorrected table")
    keep = ~both.year.between(*INWATER_YEARS)
    control_max_abs = float((both.multiplier[keep] - both.multiplier_old[keep]).abs().max())
    if control_max_abs != 0.0:
        raise ValueError(f"Satellite-era multipliers changed (max |delta| {control_max_abs})")

    out.mkdir(parents=True)
    index_path, table_path = out / "chla_index_p90_rma_corrected.csv", out / "larval_survival_ap_rma_k.csv"
    rma[["RME_UNIQUE_ID", "reef_id", "year", "source", "conversion", "p90_iw", "p90_iw_corrected", "gm_p90",
         "log_anom_p90_uncorrected", "k", "log_anom_p90"]].to_csv(index_path, index=False)
    table.to_csv(table_path, index=False)

    reef_new = values.groupby(level="year").median()
    reef_old = rma.set_index(["RME_UNIQUE_ID", "year"]).log_anom_p90_uncorrected \
                  .map(lambda a: float(np.minimum(CAP, np.exp(B * a)))).groupby(level="year").median()
    med_new = rma.groupby("year").log_anom_p90.median()
    sd = {"1993_2002_uncorrected": float(med.loc[INWATER_YEARS[0]:INWATER_YEARS[1]].std()),
          "1993_2002_corrected": float(med_new.loc[INWATER_YEARS[0]:INWATER_YEARS[1]].std()),
          "2003_2026": float(med.loc[2003:].std())}

    fig, axes = plt.subplots(1, 2, figsize=(14, 4.6))
    for ax, (old_s, new_s, ylab, ref) in zip(axes, ((med, med_new, "Domain-median p90 log anomaly", 0.0),
                                                    (reef_old, reef_new, "Reef-median survival multiplier", 1.0))):
        ax.axvspan(INWATER_YEARS[0] - 0.5, INWATER_YEARS[1] + 0.5, color="#e4e3df", alpha=0.6, linewidth=0,
                   label="in-water era (regional, spliced)")
        ax.axhline(ref, color="#52514e", linewidth=0.8, linestyle=":")
        ax.plot(old_s.index, old_s.values, color="#86b6ef", linewidth=1.5, marker="o", markersize=3.5,
                label="before correction")
        ax.plot(new_s.index, new_s.values, color="#0d366b", linewidth=2, marker="o", markersize=3.5,
                label=f"1990s scaled by k = {k:.2f}")
        top = new_s.loc[list(PEAKS_SATELLITE)].mean()
        ax.axhline(top, color="#eb6834", linewidth=0.9, linestyle="--", label="mean of 2011 and 2024 peaks")
        for y in PEAKS_1990S + PEAKS_SATELLITE:
            ax.annotate(f"{new_s.loc[y]:.2f}", (y, new_s.loc[y]), textcoords="offset points", xytext=(0, 7),
                        ha="center", fontsize=7.5, color="#0b0b0b")
        ax.set_xlabel("Wet season (spawning year)")
        ax.set_ylabel(ylab)
        ax.grid(color="#e4e3df", linewidth=0.8)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_title("A  Chlorophyll p90 anomaly (RMA, in-water scale)", loc="left", fontsize=10)
    axes[1].set_title("B  Larval-survival multiplier min(5, exp(3 x anomaly))", loc="left", fontsize=10)
    axes[1].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    figure.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(figure, dpi=160, facecolor="#fcfcfb")
    plt.close(fig)

    meta = {"status": "derived_input_not_promoted", "k": k,
            "k_definition": f"mean domain-median p90 log anomaly {PEAKS_SATELLITE} / mean {PEAKS_1990S} (RMA)",
            "applied_to": f"log_anom_p90, RMA conversion, seasons {INWATER_YEARS[0]}-{INWATER_YEARS[1]}, both signs",
            "domain_median_anomaly_sd": sd, "cap": CAP, "b": B,
            "control_satellite_era_max_abs_delta": control_max_abs,
            "reef_median_multiplier": {"uncorrected": reef_old.round(3).to_dict(), "corrected": reef_new.round(3).to_dict()},
            "inputs": {str(p.relative_to(ROOT)): sha256(p) for p in (INDEX, UNCORRECTED, SITES)},
            "script_sha256": sha256(Path(__file__)), "figure": str(figure.relative_to(ROOT)),
            "outputs": {p.name: sha256(p) for p in (index_path, table_path)} | {figure.name: sha256(figure)}}
    (out / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(f"k = {k:.4f}; domain-median anomaly sd {sd}")
    print(pd.DataFrame({"anom_old": med, "anom_new": med_new, "mult_old": reef_old, "mult_new": reef_new})
          .loc[1993:2002].round(3).to_string())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--figure", type=Path, default=DEFAULT_FIG)
    a = parser.parse_args()
    main(a.output.resolve(), a.figure.resolve())
