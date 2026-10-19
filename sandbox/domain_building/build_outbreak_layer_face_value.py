"""Outbreak-probability layer for option (b): the diver/cull correction enters external supply at
face value, not boosted by the external multiplier (lead's decision, 2026-10-15, sixteenth round;
crash_mechanism_protocol.md, "Preregistered face-value correction (option b)").

External supply is linear in the source outbreak probability P, and ADRIA multiplies it by Ext (x5).
With
    P* = P_orig + (P_corr - P_orig) / Ext
the supply under Ext is Ext * P_orig + (P_corr - P_orig): the manta-based predictions keep the
x5 boost and the correction counts once. P* lies in [0, 1] because P_corr in [0, 1].

--weight w (added 2026-10-16, eighteenth round; default 1 reproduces the face-value layer) gives
    P_w = P_orig + w (P_corr - P_orig) / Ext,   supply under Ext = Ext * P_orig + w (P_corr - P_orig),
i.e. the correction counts w times (w = Ext reproduces the fully boosted corrected layer). Lead's
reason for 1 < w < Ext: cull CPUE reflects densities before culling, so the full conversion may be
too strong. P_w stays in [0, 1] for 0 <= w <= Ext.

Inputs: sandbox/data/gbrPredsAdj_20262408.csv + gbrPreds2026Adj_20262408.csv (original) and
sandbox/data/derived/outbreak_layer_corrected_2018_2026/gbrPredsCorrected_1991_2026.csv.
Output (licensed inputs, git-ignored): sandbox/data/derived/outbreak_layer_face_value_x5_2018_2026/.

--corrected (added 2026-10-19, twenty-fifth round) selects the corrected layer CSV (default: the binary
layer above; e.g. the graded layer from build_corrected_outbreak_layer.py --graded).

    python build_outbreak_layer_face_value.py [--ext 5] [--weight 1] [--corrected <csv>] [--output <dir>]
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "sandbox" / "data"
ORIGINAL = [DATA / "gbrPredsAdj_20262408.csv", DATA / "gbrPreds2026Adj_20262408.csv"]
CORRECTED = DATA / "derived" / "outbreak_layer_corrected_2018_2026" / "gbrPredsCorrected_1991_2026.csv"
DEFAULT_OUT = DATA / "derived" / "outbreak_layer_face_value_x5_2018_2026"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main(ext: float, out: Path, weight: float = 1.0, corrected: Path = CORRECTED) -> None:
    if not 0.0 <= weight <= ext:
        raise ValueError("weight must lie in [0, ext]")
    if out.exists():
        raise FileExistsError(f"Refusing to overwrite {out}")
    orig = pd.concat([pd.read_csv(p)[["year", "reefName", "outbrProb", "x", "y"]] for p in ORIGINAL], ignore_index=True)
    corr = pd.read_csv(corrected)[["year", "reefName", "outbrProb"]]
    for name, d in (("original", orig), ("corrected", corr)):
        if d.duplicated(["year", "reefName"]).any():
            raise ValueError(f"Duplicate reef-years in the {name} layer")
    m = orig.merge(corr, on=["year", "reefName"], suffixes=("_orig", "_corr"), validate="one_to_one", how="outer",
                   indicator=True)
    if (m._merge != "both").any():
        raise ValueError("Original and corrected layers cover different reef-years")
    m["outbrProb"] = m.outbrProb_orig + weight * (m.outbrProb_corr - m.outbrProb_orig) / ext
    if not np.all((m.outbrProb >= 0.0) & (m.outbrProb <= 1.0)):
        raise ValueError("P* outside [0, 1]")
    unchanged = m.year < 2018
    if float((m.outbrProb[unchanged] - m.outbrProb_orig[unchanged]).abs().max()) != 0.0:
        raise ValueError("Years before 2018 changed")
    out.mkdir(parents=True)
    path = out / "gbrPredsFaceValue_1991_2026.csv"
    m.sort_values(["year", "reefName"])[["year", "reefName", "outbrProb", "x", "y"]].to_csv(path, index=False)
    summary = m.groupby("year")[["outbrProb_orig", "outbrProb_corr", "outbrProb"]].mean().loc[2017:]
    meta = {"status": "derived_input_not_promoted", "ext": ext, "weight": weight,
            "formula": "P* = P_orig + weight (P_corr - P_orig) / ext",
            "gbr_mean_probability_by_year": json.loads(summary.round(4).to_json(orient="index")),
            "inputs": {str(p.relative_to(ROOT)): sha256(p) for p in [*ORIGINAL, corrected]},
            "script_sha256": sha256(Path(__file__)), "outputs": {path.name: sha256(path)}}
    (out / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(summary.round(4).to_string())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ext", type=float, default=5.0)
    parser.add_argument("--weight", type=float, default=1.0)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--corrected", type=Path, default=CORRECTED)
    a = parser.parse_args()
    main(a.ext, a.output.resolve(), a.weight, a.corrected.resolve())
