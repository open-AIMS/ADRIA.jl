"""Evaluate the amended qualitative gate G2' (crash_mechanism_protocol.md, 2026-10-13)
alongside the old gate, on detectability re-score outputs (rescore_detectability.jl:
raw and smoothed_3y scores plus detected-CPUE series at every A_half).

G2' on Lizard, MacGillivray and North Direction:
  1. two detected model peaks on each;
  2. all observed peaks matched under raw (7/7) and under smoothed_3y (6/6);
  3. smoothed model minimum between the first two peaks <= 0.5 x the smaller peak;
  4. mean model CPUE over the observed inter-wave years (strictly between the two
     observed raw peaks, observed CPUE < 0.22) <= upper 95% Poisson (Garwood) bound of
     the observed mean (sum counts / sum tows).

When the re-scores carry frozen score v2 columns (rescore_detectability.jl, 2026-10-14), the
summed v2 loss and v2 matched/observed peak counts are reported alongside; G2' is unchanged.

Run with base miniconda Python (scipy):
    python evaluate_gate_g2prime.py <output_dir> <rescore_dir> [<rescore_dir> ...]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2

sys.path.insert(0, str(Path(__file__).resolve().parent))
from audit_1991_switching_repair import REEFS, TWO_WAVE, sha256  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
OBS = ROOT / "sandbox/data/reef_cots.csv"
OUTBREAK_CPUE = 0.22
SEPARATION = 0.5


def peaks(text) -> list[float]:
    return [float(v) for v in str(text).split(";") if v and v != "nan"]


def observed_troughs(raw_scores: pd.DataFrame, observations: pd.DataFrame) -> dict:
    """Observed inter-wave years per two-wave reef, from the observed raw peaks."""
    out = {}
    for reef in TWO_WAVE:
        obs_peaks = sorted({int(y) for y in raw_scores[raw_scores.reef_name == reef].obs_peak_years.iloc[0].split(";")})
        if len(obs_peaks) < 2:
            raise ValueError(f"{reef} lacks two observed raw peaks")
        x = observations[observations.reef_name == REEFS[reef]]
        window = x[(x.year > obs_peaks[0]) & (x.year < obs_peaks[1]) & (x.cotsptow < OUTBREAK_CPUE)]
        n, tows = int(window.cots.sum()), float(window.tows.sum())
        out[reef] = {"observed_peaks": obs_peaks, "years": sorted(window.year.astype(int).tolist()),
                     "counts": n, "tows": tows, "observed_mean_cpue": n / tows,
                     "poisson_upper_95": float(chi2.ppf(.975, 2 * (n + 1)) / 2 / tows),
                     "poisson_lower_95": float(chi2.ppf(.025, 2 * n) / 2 / tows) if n > 0 else 0.0}
    return out


def main(output: Path, rescores: list[Path], trough_reefs: set[str] | None = None) -> None:
    # Lead's decision (2026-10-14): the trough test may be restricted to a subset of the
    # two-wave reefs (Lizard, MacGillivray); North Direction is then reported, not gated.
    trough_reefs = set(TWO_WAVE) if trough_reefs is None else trough_reefs
    if not trough_reefs <= set(TWO_WAVE):
        raise ValueError(f"Unknown trough reefs: {trough_reefs - set(TWO_WAVE)}")
    files = {name: output / name for name in ("gate_g2prime.csv", "gate_g2prime.json")}
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite {output}")
    scores = pd.concat([pd.read_csv(r / "detectability_scores.csv") for r in rescores], ignore_index=True)
    series = pd.concat([pd.read_csv(r / "detectability_series.csv") for r in rescores], ignore_index=True)
    observations = pd.read_csv(OBS)
    raw = scores[scores.observation_treatment == "raw"]
    smooth = scores[scores.observation_treatment == "smoothed_3y"]
    troughs = observed_troughs(raw, observations)
    key = ["run", "treatment", "a_half"]
    raw_idx = raw.set_index(key + ["reef_name"]).sort_index()
    has_v2 = "loss_v2" in scores.columns
    if has_v2:
        primary = scores[scores.score_v2_primary.astype(str).str.lower() == "true"]
        if primary.groupby(key).size().ne(len(REEFS)).any():
            raise ValueError("Expected one score v2 primary row per reef")
        v3_cols = {"summed_loss_v3": ("loss_v3", "sum"), "summed_diver_error": ("diver_error", "sum"),
                   "summed_loss_v2_manta": ("loss_v2_manta", "sum")} if "loss_v3" in scores.columns else {}
        v2 = primary.groupby(key).agg(
            **v3_cols,
            summed_loss_v2=("loss_v2", "sum"), matched_v2=("matched_peaks", "sum"),
            observed_v2=("obs_peak_years", lambda v: sum(len(str(x).split(";")) for x in v)))
    rows = []
    for (run, arm, a_half), group in smooth.groupby(key):
        row = {"run": run, "treatment": arm, "a_half": a_half,
               "summed_loss_smoothed": float(group.loss.sum()),
               "matched_smoothed": int(group.matched_peaks.sum()),
               "observed_smoothed": int(sum(len(str(v).split(";")) for v in group.obs_peak_years)),
               "matched_raw": int(raw_idx.loc[(run, arm, a_half)].matched_peaks.sum()),
               "observed_raw": int(sum(len(str(v).split(";")) for v in raw_idx.loc[(run, arm, a_half)].obs_peak_years))}
        if has_v2:
            row |= {k: (int(v) if k in ("matched_v2", "observed_v2") else float(v)) for k, v in v2.loc[(run, arm, a_half)].items()}
        two_waves, distinct, trough_ok, old_ratio = True, True, True, []
        for s in group[group.reef_name.isin(TWO_WAVE)].itertuples(index=False):
            short = {"Lizard Island Reef": "liz", "MacGillivray Reef": "macg", "North Direction Reef": "nd"}[s.reef_name]
            years = [int(y) for y in peaks(s.sim_peak_years)]
            heights = peaks(s.sim_peak_heights_cpue)
            x = series[(series.run == run) & (series.treatment == arm) & (series.a_half == a_half) &
                       (series.reef_name == s.reef_name) & (series.year >= 1992)].sort_values("year").reset_index(drop=True)
            cpue = x.set_index("year").detected_cpue
            ratio = np.nan
            if len(years) >= 2:
                sm = x.detected_cpue.rolling(3, center=True, min_periods=1).mean()
                ratio = float(sm[x.year.between(years[0], years[1])].min() / min(heights[:2]))
            else:
                two_waves = False
            window = troughs[s.reef_name]
            model_mean = float(cpue.loc[window["years"]].mean())
            reef_trough_ok = model_mean <= window["poisson_upper_95"]
            distinct &= bool(np.isfinite(ratio) and ratio <= SEPARATION)
            if s.reef_name in trough_reefs:
                trough_ok &= reef_trough_ok
            old_ratio.append(ratio)
            row |= {f"{short}_peaks": ";".join(f"{y}:{h:.2f}" for y, h in zip(years, heights)),
                    f"{short}_separation_ratio": ratio,
                    f"{short}_trough_model_mean_cpue": model_mean,
                    f"{short}_trough_observed_mean_cpue": window["observed_mean_cpue"],
                    f"{short}_trough_observed_upper_95": window["poisson_upper_95"],
                    f"{short}_trough_ok": reef_trough_ok}
        finite = [r for r in old_ratio if np.isfinite(r)]
        row["worst_separation_ratio"] = max(finite) if len(finite) == len(TWO_WAVE) else np.nan
        row["c1_two_waves"] = two_waves
        row["c2_peaks_matched_raw_and_smoothed"] = (row["matched_raw"] == row["observed_raw"] and
                                                    row["matched_smoothed"] == row["observed_smoothed"])
        row["c3_distinct_waves"] = two_waves and distinct
        row["c4_trough_consistent"] = trough_ok
        row["passes_g2prime"] = (row["c1_two_waves"] and row["c2_peaks_matched_raw_and_smoothed"] and
                                 row["c3_distinct_waves"] and row["c4_trough_consistent"])
        row["passes_old_gate_raw_count"] = (two_waves and row["matched_raw"] == row["observed_raw"] and
                                            bool(np.isfinite(row["worst_separation_ratio"]) and
                                                 row["worst_separation_ratio"] <= .10))
        rows.append(row)
    result = pd.DataFrame(rows)
    output.mkdir(parents=True)
    result.to_csv(files["gate_g2prime.csv"], index=False)
    report = {
        "status": "gate_evaluation_not_promoted",
        "gate": "G2' (crash_mechanism_protocol.md, 2026-10-13)",
        "constants": {"outbreak_cpue": OUTBREAK_CPUE, "separation": SEPARATION, "trough_bound": "Garwood 95% upper"},
        "trough_gated_reefs": sorted(trough_reefs),
        "observed_troughs": troughs,
        "passing": result[result.passes_g2prime][key].to_dict(orient="records"),
        "inputs": {f"{r.name}/{n}": sha256(r / n) for r in rescores for n in (
            "metadata.toml", "detectability_scores.csv", "detectability_series.csv")} |
                  {"reef_cots.csv": sha256(OBS)},
        "script_sha256": sha256(Path(__file__)),
        "outputs": {"gate_g2prime.csv": sha256(files["gate_g2prime.csv"])},
    }
    files["gate_g2prime.json"].write_text(json.dumps(report, indent=2, sort_keys=True, default=str) + "\n")
    pd.set_option("display.width", 260)
    pd.set_option("display.max_columns", 40)
    shown = result[result.a_half.isin([0.0, 10.0])]
    print(json.dumps(troughs, indent=1))
    v2_cols = ["summed_loss_smoothed", "summed_loss_v2", "matched_v2"] if has_v2 else []
    v2_cols += ["summed_loss_v3"] if "loss_v3" in scores.columns else []
    print(shown[key + v2_cols + ["matched_raw", "matched_smoothed", "worst_separation_ratio",
                       "liz_trough_model_mean_cpue", "macg_trough_model_mean_cpue", "nd_trough_model_mean_cpue",
                       "c1_two_waves", "c2_peaks_matched_raw_and_smoothed", "c3_distinct_waves",
                       "c4_trough_consistent", "passes_g2prime", "passes_old_gate_raw_count"]]
          .round(3).to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("rescores", type=Path, nargs="+")
    parser.add_argument("--trough-reefs", nargs="+", default=None,
                        help="model reef names whose trough is gated (default: all two-wave reefs)")
    args = parser.parse_args()
    main(args.output.resolve(), [r.resolve() for r in args.rescores],
         None if args.trough_reefs is None else set(args.trough_reefs))
