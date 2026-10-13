"""Empirical-Bayes 1991 COTS starting surface (crash_mechanism_protocol.md, 2026-10-13).

Each donor reef's pooled 1989-1991 manta rate c/t is replaced by its Poisson-gamma
posterior mean (c + alpha) / (t + beta). The gamma prior is fitted by negative-binomial
marginal maximum likelihood to the regional donors (within RADIUS_KM of any domain
site). The IDW of build_1991_initial_surfaces.py is then re-run with unchanged weights
(tows / max(distance, 1 km)^2), radius and minimum donors. Only the COTS column changes;
the original is kept as cots_cpue_1989_1991_idw_unshrunk.

Run with base miniconda Python (scipy):
    python build_1991_eb_surface.py --source runs/20261007T199104_1991_initial_surfaces --output <new dir>
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import gammaln

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_1991_initial_surfaces import MIN_DONORS, RADIUS_KM, ROOT, haversine_km, idw, sha256  # noqa: E402

COTS_COLUMN = "cots_cpue_1989_1991_idw"


def fit_gamma_prior(counts: np.ndarray, tows: np.ndarray) -> tuple[float, float, float]:
    """Negative-binomial marginal MLE of a gamma(alpha, rate beta) prior on COTS/tow."""
    def nll(log_params: np.ndarray) -> float:
        a, b = np.exp(log_params)
        return -float(np.sum(gammaln(counts + a) - gammaln(a) - gammaln(counts + 1) +
                             a * np.log(b / (b + tows)) + counts * np.log(tows / (b + tows))))
    result = minimize(nll, np.array([0.0, np.log(50.0)]), method="Nelder-Mead",
                      options={"xatol": 1e-8, "fatol": 1e-10, "maxiter": 5000})
    if not result.success:
        raise RuntimeError(f"Prior fit failed: {result.message}")
    alpha, beta = np.exp(result.x)
    return float(alpha), float(beta), float(result.fun)


def build(source: Path, output: Path) -> None:
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite {output}")
    surface = pd.read_csv(source / "surface.csv")
    donors = pd.read_csv(source / "cots_donors.csv")
    if (donors.tows <= 0).any() or (donors.cots < 0).any():
        raise ValueError("Invalid donor counts or effort")
    lon = surface.longitude.to_numpy(dtype=float)
    lat = surface.latitude.to_numpy(dtype=float)
    distances = haversine_km(lon, lat, donors.longitude.to_numpy(dtype=float),
                             donors.latitude.to_numpy(dtype=float))
    # Reproduce the stored surface first, so the IDW call is known to be identical.
    check, _, _ = idw(distances, donors.cots_per_tow.to_numpy(dtype=float), donors.tows.to_numpy(dtype=float))
    stored = surface[COTS_COLUMN].to_numpy(dtype=float)
    if not np.allclose(check, stored, rtol=0, atol=1e-12, equal_nan=True):
        raise ValueError("Could not reproduce the source IDW surface")
    regional = (distances <= RADIUS_KM).any(axis=0)
    counts = donors.cots.to_numpy(dtype=float)
    tows = donors.tows.to_numpy(dtype=float)
    alpha, beta, nll = fit_gamma_prior(counts[regional], tows[regional])
    donors = donors.assign(regional_prior_donor=regional,
                           cots_per_tow_unshrunk=donors.cots_per_tow,
                           cots_per_tow_eb=(counts + alpha) / (tows + beta))
    eb, count, furthest = idw(distances, donors.cots_per_tow_eb.to_numpy(dtype=float), tows)
    if np.isnan(eb).any() or not np.isfinite(eb).all():
        raise ValueError("EB surface does not cover every site")
    out = surface.copy()
    out[f"{COTS_COLUMN}_unshrunk"] = stored
    out[COTS_COLUMN] = eb
    output.mkdir(parents=True)
    out.to_csv(output / "surface.csv", index=False)
    donors.to_csv(output / "cots_donors_eb.csv", index=False)
    metadata = {
        "status": "experimental_not_promoted",
        "method": "Poisson-gamma empirical Bayes per donor reef, then IDW unchanged",
        "posterior_mean": "(cots + alpha) / (tows + beta)",
        "prior": {"alpha": alpha, "beta_tows": beta, "mean_cots_per_tow": alpha / beta,
                  "fit": "negative-binomial marginal maximum likelihood (Nelder-Mead on log parameters)",
                  "negative_log_likelihood": nll,
                  "donors": int(regional.sum()),
                  "donor_rule": f"donor reefs within {RADIUS_KM} km of any domain site"},
        "interpolation": {"radius_km": RADIUS_KM, "min_distinct_donor_reefs": MIN_DONORS,
                          "cots_weight": "tows / max(distance_km, 1)^2 (unchanged)"},
        "source_surface_reproduced_max_abs": float(np.nanmax(np.abs(check - stored))),
        "years": [1989, 1990, 1991],
        "units": "COTS per manta tow",
        "columns": {COTS_COLUMN: "EB-shrunk IDW (used by the runner)",
                    f"{COTS_COLUMN}_unshrunk": "original IDW surface"},
        "inputs": {str(path.relative_to(ROOT)).replace("\\", "/"): sha256(path) for path in (
            source / "surface.csv", source / "cots_donors.csv", source / "metadata.json",
            Path(__file__), Path(__file__).with_name("build_1991_initial_surfaces.py"))},
        "outputs": {name: sha256(output / name) for name in ("surface.csv", "cots_donors_eb.csv")},
        "dimensions": {"sites": len(out), "donors": len(donors)},
    }
    (output / "metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    print(json.dumps(metadata["prior"], indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    build(args.source.resolve(), args.output.resolve())
