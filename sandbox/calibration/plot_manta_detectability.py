"""Explain the preregistered density-dependent manta detection curve
(crash_mechanism_protocol.md, 2026-10-13):

    CPUE = q * A * A / (A + A_half),  q = 0.015 COTS/tow per adult/ha

Left: relative detection d(A) = A / (A + A_half) against true adult density.
Right: expected manta COTS/tow against true density, with the frozen linear
model (A_half = 0) for comparison. Reference marks: model trough densities,
the 0.22 COTS/tow outbreak level and the observed trough counts."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

Q = 0.015
A_HALF = (1, 3, 10, 30)
# Reference ordinal blue ramp, steps 250/400/550/700 (validated: dataviz --ordinal).
RAMP = ("#86b6ef", "#3987e5", "#1c5cab", "#0d366b")
CONTROL = "#52514e"
INK, MUTED, GRID, BAND = "#0b0b0b", "#52514e", "#e4e3df", "#f0efec"
TROUGH = (2.0, 3.0)        # model trough adults/ha (runs 11-12)
OUTBREAK_CPUE = 0.22       # COTS/tow outbreak level
OBS_TROUGH_CPUE = 0.03     # observed trough counts 0-0.03 COTS/tow (LTMP 1999-2010)


def cpue(a: np.ndarray, a_half: float) -> np.ndarray:
    return Q * a if a_half == 0 else Q * a * a / (a + a_half)


def main(output: Path) -> None:
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite {output}")
    a = np.logspace(-1, np.log10(300), 400)
    fig, (left, right) = plt.subplots(1, 2, figsize=(13, 5.2))
    for axis in (left, right):
        axis.axvspan(*TROUGH, color=BAND, zorder=0)
        axis.set_xscale("log")
        axis.set_xlim(0.1, 300)
        axis.grid(color=GRID, linewidth=.8)
        axis.tick_params(colors=MUTED)
        for spine in axis.spines.values():
            spine.set_color(GRID)
        axis.set_xlabel("True adult density A (COTS / ha)", color=INK)

    left.axhline(1.0, color=CONTROL, linestyle="--", linewidth=2, label="frozen linear model (AÂ½ = 0)")
    for a_half, color in zip(A_HALF, RAMP):
        d = a / (a + a_half)
        left.plot(a, d, color=color, linewidth=2, label=f"AÂ½ = {a_half} /ha")
        left.annotate(f"AÂ½ {a_half}", (a_half, .5), xytext=(4, -12), textcoords="offset points",
                      fontsize=8, color=INK)
        left.plot([a_half], [.5], "o", color=color, markersize=8, markeredgecolor="white", markeredgewidth=2)
    left.axhline(.5, color=GRID, linewidth=1)
    left.text(2.45, 1.04, "model troughs\n2â€“3 /ha", ha="center", va="bottom", fontsize=8, color=MUTED)
    left.set_ylim(0, 1.15)
    left.set_ylabel("Relative manta detection d(A) = A / (A + AÂ½)", color=INK)
    left.set_title("Detection falls at low density, saturates at outbreak density",
                   loc="left", fontsize=11, color=INK)

    right.plot(a, cpue(a, 0), color=CONTROL, linestyle="--", linewidth=2, label="frozen linear model (AÂ½ = 0)")
    for a_half, color in zip(A_HALF, RAMP):
        right.plot(a, cpue(a, a_half), color=color, linewidth=2, label=f"AÂ½ = {a_half} /ha")
    right.axhspan(1e-4, OBS_TROUGH_CPUE, color=BAND, zorder=0)
    right.text(0.12, 0.012, "observed trough counts\n(0â€“0.03 COTS/tow)", fontsize=8, color=MUTED, va="center")
    right.axhline(OUTBREAK_CPUE, color=MUTED, linestyle=":", linewidth=1)
    right.text(0.12, OUTBREAK_CPUE * 1.15, "outbreak level 0.22 COTS/tow", fontsize=8, color=MUTED)
    # Worked example at a model trough of 3 adults/ha.
    for a_half, color, name in ((0, CONTROL, "linear"), (10, RAMP[2], "A½ 10")):
        y = float(cpue(np.array([3.0]), a_half)[0])
        right.plot([3.0], [y], "o", color=color, markersize=8, markeredgecolor="white", markeredgewidth=2, zorder=5)
        right.annotate(f"{name}: {y:.3f} at 3 /ha", (3.0, y), xytext=(8, -3), textcoords="offset points", fontsize=8, color=INK)
    right.set_yscale("log")
    right.set_ylim(1e-4, 5)
    right.set_ylabel("Expected manta COTS / tow (q = 0.015)", color=INK)
    right.set_title("Same model troughs read as near-zero manta counts",
                    loc="left", fontsize=11, color=INK)
    handles, labels = left.get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=5, fontsize=9, frameon=False,
               bbox_to_anchor=(.5, -.02))
    fig.suptitle("Density-dependent manta detection: CPUE = q Â· A Â· A / (A + AÂ½)",
                 x=.01, ha="left", fontsize=13, color=INK)
    fig.tight_layout(rect=(0, .06, 1, .95))
    fig.savefig(output, dpi=170, bbox_inches="tight", facecolor="#fcfcfb")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    main(parser.parse_args().output)
