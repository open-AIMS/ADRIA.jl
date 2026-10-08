"""Validate one immutable Owen flux audit against the frozen paired gate."""

from __future__ import annotations

import argparse
import hashlib
import tomllib
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
GATE = ROOT / "sandbox/calibration/runs/20261001T113401_lizard_domain_gate/trajectories.csv"


def digest(path: Path) -> str:
    sha = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            sha.update(block)
    return sha.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-id", required=True)
    args = parser.parse_args()
    if not args.run_id.replace("_", "").replace("-", "").isalnum():
        raise ValueError("Invalid run ID")
    run = ROOT / "sandbox/calibration/runs" / args.run_id
    with (run / "metadata.toml").open("rb") as stream:
        metadata = tomllib.load(stream)
    if metadata["hydrodynamic_mode"] != "sample" or metadata["seed"] != 20260930:
        raise AssertionError("Audit is not the sampled-year paired-gate replay")
    if max(metadata["max_absolute_residual_ha"].values()) >= 1e-7:
        raise AssertionError("Flux accounting residual exceeds tolerance")
    for name, expected in metadata["output_sha256"].items():
        if digest(run / name) != expected:
            raise AssertionError(f"Output hash changed: {name}")
    for name, expected in metadata["source_sha256"].items():
        if digest(Path(name)) != expected:
            raise AssertionError(f"Input hash changed: {name}")
    for name, expected in metadata["code_sha256"].items():
        if digest(Path(name)) != expected:
            raise AssertionError(f"Code hash changed: {name}")

    source = pd.read_csv(run / "source_survival_checks.csv")
    if len(source) != 113 or source.reef_id.nunique() != 113:
        raise AssertionError("Unexpected local source coverage")
    np.testing.assert_allclose(source.rme_2017, source.applied_2018, rtol=1e-10, atol=1e-14)

    audit = pd.read_csv(run / "target_reef_fluxes.csv")
    gate = pd.read_csv(GATE)
    gate = gate[(gate.seed == 20260930) & (gate.treatment == "v2_owen_boundary")]
    merged = audit.merge(gate, on=["reef_name", "year"], validate="one_to_one")
    if len(merged) != 156 or merged.reef_name.nunique() != 4:
        raise AssertionError("Incomplete four-reef replay")
    np.testing.assert_array_equal(merged.forcing_year_x, merged.forcing_year_y)
    np.testing.assert_allclose(merged.adult_stage_ha, merged.adults_ha, rtol=0, atol=1e-12)
    np.testing.assert_allclose(merged.coral_cover_x, merged.coral_cover_y, rtol=0, atol=1e-12)
    print(f"PASS {args.run_id}: 113 sources, 156 reef-years, archived control and hashes")


if __name__ == "__main__":
    main()
