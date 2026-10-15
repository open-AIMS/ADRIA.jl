"""Plot the proposed chlorophyll-driven larval-survival multiplier (lead's decision, 2026-10-15,
eleventh round) against the lead's original Fabricius-based survival function.

Proposed (relative, per reef and wet season):
    multiplier = min(CAP, (chl / reef_normal)^b) = min(CAP, exp(b * log_anom)),  b in {0, 1.5, 3}
Original (absolute; COTSModel_ChlorophyllModel.R): logit(survival) ~ log2(chl), Fabricius et al.
(2010) laboratory data (sandbox/data/external/COTS_Model/Data/ChlorophyllModel/Fabricius.csv,
rows with keep = TRUE).

    python plot_larval_survival_curve.py <output.png>
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
FAB = ROOT / "sandbox/data/external/COTS_Model/Data/ChlorophyllModel/Fabricius.csv"
INDEX = ROOT / "sandbox/data/derived/chla_index_lizard_2003_2026/chla_index_reef_year.csv"
CAP = 5.0
B_VALUES = (0.0, 1.5, 3.0)
COLORS = {0.0: "#9a9893", 1.5: "#3987e5", 3.0: "#0d366b"}
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"


def main(out: Path) -> None:
    fab = pd.read_csv(FAB)
    fab = fab[fab.keep.astype(str).str.upper().str.startswith("T")]
    y = np.clip(fab.surv / 100, 0.001, 0.999)
    slope, intercept = np.polyfit(np.log2(fab.chl), np.log(y / (1 - y)), 1)
    index = pd.read_csv(INDEX)
    index = index[np.isfinite(index.log_anom_mean)]

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8))
    # A: proposed relative multiplier
    ax = axes[0]
    ratio = np.exp(np.linspace(np.log(0.4), np.log(2.5), 300))
    for b in B_VALUES:
        ax.plot(ratio, np.minimum(CAP, ratio ** b), color=COLORS[b], linewidth=2.2, label=f"b = {b:g}" + (" (off)" if b == 0 else ""))
    q = np.exp(np.quantile(index.log_anom_mean, [0.05, 0.95]))
    ax.axvspan(q[0], q[1], color="#cde2fb", alpha=0.5, linewidth=0, label="90% of reef-seasons (mean)")
    ax.axhline(1, color=MUTED, linewidth=0.8, linestyle=":")
    ax.set_xscale("log", base=2)
    ax.set_xticks([0.5, 0.75, 1, 1.5, 2])
    ax.set_xticklabels(["0.5x", "0.75x", "1x", "1.5x", "2x"])
    ax.set_xlabel("Wet-season chl-a relative to the reef's normal")
    ax.set_ylabel("Larval survival multiplier")
    ax.set_title("A  Proposed: gradual, capped at 5x", loc="left", fontsize=10, color=INK)
    ax.legend(frameon=False, fontsize=8)
    # B: original absolute curve
    ax = axes[1]
    chl = np.exp(np.linspace(np.log(0.02), np.log(6), 300))
    surv = 1 / (1 + np.exp(-(intercept + slope * np.log2(chl))))
    ax.plot(chl, surv, color="#eb6834", linewidth=2.2, label="Original logit fit")
    ax.scatter(fab.chl, fab.surv / 100, s=22, color=INK, zorder=3, label="Fabricius et al. 2010 (lab)")
    ax.axvspan(0.2, 0.5, color="#cde2fb", alpha=0.5, linewidth=0, label="Typical northern GBR in situ")
    ax.set_xscale("log", base=2)
    ax.set_xlabel("Chlorophyll-a (ug/L, absolute)")
    ax.set_ylabel("Larval survival (proportion)")
    ax.set_title("B  Original: near-zero below ~0.5 ug/L", loc="left", fontsize=10, color=INK)
    ax.legend(frameon=False, fontsize=8)
    # C: domain-median multiplier by season
    ax = axes[2]
    med = index.groupby("year")[["log_anom_mean", "log_anom_p90"]].median()
    for b in (1.5, 3.0):
        ax.plot(med.index, np.minimum(CAP, np.exp(b * med.log_anom_mean)), color=COLORS[b], linewidth=2.0,
                marker="o", markersize=4, label=f"b = {b:g}, wet-season mean")
        ax.plot(med.index, np.minimum(CAP, np.exp(b * med.log_anom_p90)), color=COLORS[b], linewidth=1.4,
                linestyle="--", label=f"b = {b:g}, wet-season p90")
    ax.axhline(1, color=MUTED, linewidth=0.8, linestyle=":")
    ax.axvspan(2002.5, 2024.5, color="#ffffff", alpha=0)
    for year, label in ((2011, "2011 floods"), (2019, "2019 floods"), (2024, "2024 (Jasper)")):
        ax.annotate(label, (year, float(np.exp(3 * med.loc[year, "log_anom_p90"]))), xytext=(0, 8),
                    textcoords="offset points", ha="center", fontsize=8, color=MUTED)
    ax.set_xlabel("Wet season (spawning year); 2025-26 from VIIRS")
    ax.set_ylabel("Domain-median survival multiplier")
    ax.set_title("C  Implied by the reef chl-a index, 2003-2026", loc="left", fontsize=10, color=INK)
    ax.legend(frameon=False, fontsize=7, ncol=2)
    for ax in axes:
        ax.grid(color=GRID, linewidth=0.8)
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=160, facecolor="#fcfcfb")
    print(f"logit fit: intercept {intercept:.2f}, slope {slope:.2f} per doubling; saved {out}")
    print(np.minimum(CAP, np.exp(3 * med)).round(2).loc[[2011, 2019, 2020, 2024, 2025]].to_string())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    main(parser.parse_args().output.resolve())
