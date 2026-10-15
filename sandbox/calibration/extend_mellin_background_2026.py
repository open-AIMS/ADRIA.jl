"""Extend the Mellin et al. (2019) background coral mortality table (site x year, 1991-2024;
build_mellin_background_mortality.jl) to 2026 for the Lizard v0.4 domain. Labelled stand-in
(lead's decision, 2026-10-14/15): 2025 and 2026 take each site's 2015-2024 mean.

    python extend_mellin_background_2026.py [--output <dir>]
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "sandbox" / "data" / "derived" / "mellin_background_mortality_1991_2024" / "background_mortality.csv"
DEFAULT_OUT = ROOT / "sandbox" / "data" / "derived" / "mellin_background_mortality_1991_2026"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(out: Path) -> None:
    if out.exists():
        raise FileExistsError(f"Refusing to overwrite {out}")
    table = pd.read_csv(SOURCE, dtype={"site_id": str})
    years = [str(y) for y in range(2015, 2025)]
    if not set(years) <= set(table.columns):
        raise ValueError("Source lacks 2015-2024 columns")
    stand_in = table[years].mean(axis=1)
    for year in ("2025", "2026"):
        table[year] = stand_in
    out.mkdir(parents=True)
    path = out / "background_mortality.csv"
    table.to_csv(path, index=False)
    meta = {"status": "derived_input_not_promoted", "source": str(SOURCE.relative_to(ROOT)), "source_sha256": sha256(SOURCE),
            "extension": "2025 and 2026 = site mean of 2015-2024 (labelled stand-in)",
            "script_sha256": sha256(Path(__file__)), "outputs": {"background_mortality.csv": sha256(path)}}
    (out / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps(meta, indent=1))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    main(parser.parse_args().output.resolve())
