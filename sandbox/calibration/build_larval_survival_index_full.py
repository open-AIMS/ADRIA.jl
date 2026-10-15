"""Site x year larval-survival multiplier tables from the full 1993-2026 in-water-scale
chlorophyll series (build_chla_index_full.py), for COTS_LARVAL_SURVIVAL_CSV, plus a figure of
the survival forms (crash_mechanism_protocol.md, thirteenth round).

Variants (cap 5; 1 where the index is missing, i.e. 1991-1992):
- am_ols : anomaly form on the season mean, OLS conversion:  min(5, exp(3 log_anom_mean))
- ap_rma : anomaly form on the season p90, RMA conversion:   min(5, exp(3 log_anom_p90))
- ab_rma : absolute relaxed-Fabricius form on the in-water-equivalent p90 (RMA):
           S(c) = f + (1 - f) logistic(-0.48 + 2.02 log2 c), f = 0.1;
           multiplier = min(5, S(c) / S(c_ref)), c_ref = domain median p90 1993-2026.
RMA is used for p90 because over the 2003-2007 overlap it reproduces the in-water p90 level
(0.33-0.36 vs 0.32) whereas OLS compresses it (0.26-0.28); OLS and RMA both match the mean.

    python build_larval_survival_index_full.py [--output <dir>] [--figure <png>]
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
INDEX = ROOT / "sandbox" / "data" / "derived" / "chla_index_lizard_1993_2026_inwater" / "chla_index_reef_year.csv"
SITES = ROOT / "sandbox" / "data" / "Lizard_Historical_v0.4" / "site_to_reef.csv"
FAB = ROOT / "sandbox/data/external/COTS_Model/Data/ChlorophyllModel/Fabricius.csv"
DEFAULT_OUT = ROOT / "sandbox" / "data" / "derived" / "larval_survival_index_lizard_full"
DEFAULT_FIG = ROOT / "sandbox" / "calibration" / "figures" / "larval_survival_full_series.png"
CAP, B, FLOOR = 5.0, 3.0, 0.1
INTERCEPT, SLOPE = -0.48, 2.02  # logit fit to Fabricius et al. 2010 (plot_larval_survival_curve.py)
YEARS = range(1991, 2027)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relaxed(c, f=FLOOR):
    return f + (1 - f) / (1 + np.exp(-(INTERCEPT + SLOPE * np.log2(c))))


def main(out: Path, figure: Path) -> None:
    if out.exists():
        raise FileExistsError(f"Refusing to overwrite {out}")
    idx = pd.read_csv(INDEX, dtype={"RME_UNIQUE_ID": str})
    sites = pd.read_csv(SITES, dtype={"UNIQUE_ID": str, "site_id": str})
    rma = idx[idx.conversion == "rma"]
    c_ref = float(rma.groupby("year").p90_iw.median().median())
    reef_values = {
        "am_ols": idx[idx.conversion == "ols"].set_index(["RME_UNIQUE_ID", "year"]).log_anom_mean.dropna().map(lambda a: float(np.minimum(CAP, np.exp(B * a)))),
        "ap_rma": rma.set_index(["RME_UNIQUE_ID", "year"]).log_anom_p90.dropna().map(lambda a: float(np.minimum(CAP, np.exp(B * a)))),
        "ab_rma": rma.set_index(["RME_UNIQUE_ID", "year"]).p90_iw.dropna().map(lambda c: float(np.minimum(CAP, relaxed(c) / relaxed(c_ref)))),
    }
    for name, values in reef_values.items():
        if values.index.duplicated().any():
            raise ValueError(f"Duplicate reef-years in {name}")
    out.mkdir(parents=True)
    grid = pd.MultiIndex.from_product([sites.site_id, list(YEARS)], names=["site_id", "year"]).to_frame(index=False)
    grid = grid.merge(sites[["site_id", "UNIQUE_ID"]], on="site_id")
    summary = {}
    medians = {}
    for name, values in reef_values.items():
        values = values[np.isfinite(values)]
        m = pd.to_numeric(pd.Series([values.get((u, y), np.nan) for u, y in zip(grid.UNIQUE_ID, grid.year)]),
                          errors="coerce").astype(float)
        table = grid.assign(multiplier=np.where(np.isfinite(m), m, 1.0))
        path = out / f"larval_survival_{name}.csv"
        table[["site_id", "year", "multiplier"]].to_csv(path, index=False)
        reef_median = values.groupby(level="year").median()
        medians[name] = reef_median
        summary[name] = {"reef_median_by_year": reef_median.round(3).to_dict(), "sha256": sha256(path)}

    # Figure: the absolute forms and the implied multipliers by season.
    fab = pd.read_csv(FAB)
    fab = fab[fab.keep.astype(str).str.upper().str.startswith("T")]
    fig, axes = plt.subplots(1, 2, figsize=(14, 4.8))
    ax = axes[0]
    c = np.exp(np.linspace(np.log(0.08), np.log(3), 300))
    for f, color in ((0.0, "#eb6834"), (0.05, "#86b6ef"), (0.1, "#3987e5"), (0.2, "#0d366b")):
        ax.plot(c, relaxed(c, f), color=color, linewidth=2.2 if f == FLOOR else 1.5,
                label=f"floor f = {f:g}" + (" (original)" if f == 0 else " (used)" if f == FLOOR else ""))
    ax.scatter(fab.chl, fab.surv / 100, s=18, color="#0b0b0b", zorder=3, label="Fabricius et al. 2010 (lab)")
    q = np.quantile(rma.p90_iw.dropna(), [0.05, 0.95])
    ax.axvspan(q[0], q[1], color="#cde2fb", alpha=0.6, linewidth=0, label="90% of reef-seasons (in-water p90)")
    ax.axvline(c_ref, color="#52514e", linewidth=0.8, linestyle=":")
    ax.set_xscale("log", base=2)
    ax.set_xticks([0.125, 0.25, 0.5, 1, 2])
    ax.set_xticklabels(["0.125", "0.25", "0.5", "1", "2"])
    ax.set_xlabel("In-water-equivalent chlorophyll-a (ug/L, wet-season p90)")
    ax.set_ylabel("Larval survival (proportion)")
    ax.set_title("A  Absolute form: relaxed Fabricius curve", loc="left", fontsize=10)
    ax.legend(frameon=False, fontsize=8)
    ax = axes[1]
    styles = {"am_ols": ("#86b6ef", "anomaly, season mean (OLS)"), "ap_rma": ("#0d366b", "anomaly, season p90 (RMA)"),
              "ab_rma": ("#eb6834", "absolute relaxed, p90 (RMA)")}
    for name, (color, lab) in styles.items():
        s = medians[name]
        ax.plot(s.index, s.values, color=color, linewidth=1.8, marker="o", markersize=3.5, label=lab)
    ax.axhline(1, color="#52514e", linewidth=0.8, linestyle=":")
    ax.axvspan(1992.5, 2002.5, color="#e4e3df", alpha=0.6, linewidth=0, label="in-water era (regional)")
    ax.set_xlabel("Wet season (spawning year)")
    ax.set_ylabel("Reef-median survival multiplier")
    ax.set_title("B  Implied multipliers, 1993-2026", loc="left", fontsize=10)
    ax.legend(frameon=False, fontsize=8)
    for ax in axes:
        ax.grid(color="#e4e3df", linewidth=0.8)
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    figure.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(figure, dpi=160, facecolor="#fcfcfb")
    meta = {"status": "derived_input_not_promoted", "index": str(INDEX.relative_to(ROOT)), "index_sha256": sha256(INDEX),
            "sites_sha256": sha256(SITES), "cap": CAP, "b": B, "floor": FLOOR, "c_ref_p90_iw": c_ref,
            "variants": summary, "figure": str(figure.relative_to(ROOT)), "script_sha256": sha256(Path(__file__))}
    (out / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(f"c_ref = {c_ref:.3f}")
    print(pd.DataFrame(medians).round(2).to_string())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--figure", type=Path, default=DEFAULT_FIG)
    a = parser.parse_args()
    main(a.output.resolve(), a.figure.resolve())
