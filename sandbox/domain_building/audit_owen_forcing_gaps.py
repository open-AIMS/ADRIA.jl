"""Quantify unresolved Owen source and survival coverage without imputing it.

This is a read-only audit of the source matrices and RME ancillary inputs.  It
does not alter the experimental connectivity artifact or infer missing reef
areas/survival.  Results are written to a new, never-overwritten run directory.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

sys.dont_write_bytecode = True
from build_lizard_cots_connectivity_owen import (
    DATA_ROOT,
    DEFAULT_DOMAIN,
    DEFAULT_SURVIVAL,
    NEW_TO_OLD_ALIASES,
    REEF_SHAPE,
    REEF_SHAPE_SOURCES,
    RME_ID_LIST,
    SOURCE_ROOT,
    matrix_year,
    reef_id_from_name,
)

REPO_ROOT = Path(__file__).resolve().parents[2]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    parser.add_argument(
        "--run-id",
        default=datetime.now().strftime("%Y%m%dT%H%M%S") + "_owen_coverage_audit",
    )
    args = parser.parse_args()
    if not args.run_id.replace("_", "").replace("-", "").isalnum():
        raise ValueError("Run ID may contain only letters, numbers, hyphens, and underscores")
    out_dir = REPO_ROOT / "sandbox" / "calibration" / "runs" / args.run_id
    if out_dir.exists():
        raise FileExistsError(f"Refusing to overwrite {out_dir}")

    reef_shape = gpd.read_file(REEF_SHAPE)
    global_ids = [reef_id_from_name(name) for name in reef_shape["reefName"]]
    if len(global_ids) != 3861 or len(set(global_ids)) != len(global_ids):
        raise ValueError("Unexpected or duplicate global reef IDs")
    lookup = {reef_id: index for index, reef_id in enumerate(global_ids)}
    local_ids = sorted(set(pd.read_csv(args.domain / "site_to_reef.csv")["reef_id"]), key=lookup.get)
    if len(local_ids) != 113:
        raise ValueError("Expected 113 Lizard parent reefs")
    local_indices = np.asarray([lookup[reef_id] for reef_id in local_ids])

    rme_ids = pd.read_csv(RME_ID_LIST, comment="#", header=None, dtype={0: str})
    rme_area_ids = set(rme_ids.iloc[:, 0].astype(str))
    survival = pd.read_csv(DEFAULT_SURVIVAL, dtype={"id": str}).set_index("id")
    historical_years = sorted(int(col) for col in survival.columns)
    if historical_years != list(range(2010, 2018)):
        raise ValueError(f"Unexpected survival years: {historical_years}")
    mapped_ids = [reef_id if reef_id in rme_area_ids else NEW_TO_OLD_ALIASES.get(reef_id) for reef_id in global_ids]
    valid = np.asarray([
        old_id in rme_area_ids and old_id in survival.index for old_id in mapped_ids
    ])
    local_mask = np.zeros(len(global_ids), dtype=bool)
    local_mask[local_indices] = True
    if not valid[local_indices].all():
        raise ValueError("A local reef lacks survival or area")
    missing_indices = np.flatnonzero(~valid & ~local_mask)
    mapped_external_indices = np.flatnonzero(valid & ~local_mask)
    if len(missing_indices) != 56:
        raise ValueError(f"Expected 56 unmatched external sources, found {len(missing_indices)}")
    survival_matrix = survival.loc[[mapped_ids[index] for index in mapped_external_indices],
                                   [str(year) for year in historical_years]].to_numpy(dtype=float)
    if not np.isfinite(survival_matrix).all() or np.any((survival_matrix < 0) | (survival_matrix > 1)):
        raise ValueError("Invalid historical larval survival")

    matrix_paths = sorted(SOURCE_ROOT.glob("conMatCots*.csv"), key=matrix_year)
    if [matrix_year(path) for path in matrix_paths] != list(range(2018, 2023)):
        raise ValueError("Expected hydrodynamic years 2018–2022")
    sink_rows: list[dict[str, object]] = []
    missing_rows: list[dict[str, object]] = []
    nonfinite_rows: list[dict[str, object]] = []
    for path in matrix_paths:
        year = matrix_year(path)
        matrix = np.loadtxt(path, delimiter=",", dtype=np.float64)
        if matrix.shape != (len(global_ids), len(global_ids)):
            raise ValueError(f"Bad matrix dimensions: {path.name}: {matrix.shape}")
        nonfinite = np.sum(~np.isfinite(matrix), axis=1)
        if np.any((nonfinite > 0) & (nonfinite < len(global_ids))):
            raise ValueError(f"Partly nonfinite source rows in {path.name}")
        dead_rows = np.flatnonzero(nonfinite == len(global_ids))
        for index in dead_rows:
            nonfinite_rows.append({"hydrodynamic_year": year, "source_reef_id": global_ids[index],
                                   "is_unmatched_source": bool(index in missing_indices)})
        matrix[dead_rows, :] = 0
        if not np.isfinite(matrix).all() or np.any((matrix < 0) | (matrix > 1)):
            raise ValueError(f"Invalid connectivity probabilities in {path.name}")

        local_columns = matrix[:, local_indices]
        missing = local_columns[missing_indices, :]
        mapped = local_columns[mapped_external_indices, :]
        latest = survival_matrix[:, -1]
        historical_mean = survival_matrix.mean(axis=1)
        historical_min = survival_matrix.min(axis=1)
        historical_max = survival_matrix.max(axis=1)
        for j, reef_id in enumerate(local_ids):
            missing_raw = float(missing[:, j].sum())
            mapped_raw = float(mapped[:, j].sum())
            sink_rows.append({
                "hydrodynamic_year": year,
                "sink_reef_id": reef_id,
                "missing_source_raw_connectivity": missing_raw,
                "mapped_external_raw_connectivity": mapped_raw,
                "missing_fraction_of_external_raw_connectivity":
                    missing_raw / (missing_raw + mapped_raw) if missing_raw + mapped_raw > 0 else 0.0,
                "mapped_external_survival_weighted_2017": float(mapped[:, j] @ latest),
                "mapped_external_survival_weighted_2010_2017_mean": float(mapped[:, j] @ historical_mean),
                "mapped_external_survival_weighted_sourcewise_min": float(mapped[:, j] @ historical_min),
                "mapped_external_survival_weighted_sourcewise_max": float(mapped[:, j] @ historical_max),
            })
        for k, index in enumerate(missing_indices):
            missing_rows.append({
                "hydrodynamic_year": year,
                "source_reef_id": global_ids[index],
                "raw_connectivity_to_113_lizard_reefs": float(missing[k, :].sum()),
                "number_of_lizard_sinks_reached": int(np.count_nonzero(missing[k, :])),
            })
        print(f"Audited {path.name}", flush=True)

    sink_frame = pd.DataFrame(sink_rows)
    missing_frame = pd.DataFrame(missing_rows)
    nonfinite_frame = pd.DataFrame(nonfinite_rows)
    out_dir.mkdir(parents=True)
    outputs = {
        "per_sink": out_dir / "owen_source_coverage_per_sink.csv",
        "unmatched_sources": out_dir / "owen_unmatched_source_connectivity.csv",
        "nonfinite_rows": out_dir / "owen_nonfinite_source_rows.csv",
    }
    sink_frame.to_csv(outputs["per_sink"], index=False)
    missing_frame.to_csv(outputs["unmatched_sources"], index=False)
    nonfinite_frame.to_csv(outputs["nonfinite_rows"], index=False)
    metadata = {
        "purpose": "coverage and historical-survival sensitivity audit; no imputation",
        "source_rows_are_larval_sources": True,
        "raw_connectivity_units": "dimensionless Owen matrix weights; not recruits",
        "survival_weighted_units": "dimensionless connectivity times RME survival",
        "unmatched_external_source_count": len(missing_indices),
        "matched_external_source_count": len(mapped_external_indices),
        "survival_years": historical_years,
        "missing_fraction_mean": float(sink_frame["missing_fraction_of_external_raw_connectivity"].mean()),
        "missing_fraction_max": float(sink_frame["missing_fraction_of_external_raw_connectivity"].max()),
        "missing_source_nonzero_rows": int((missing_frame["raw_connectivity_to_113_lizard_reefs"] > 0).sum()),
        "builder_sha256": sha256(Path(__file__)),
        "sources": {str(path.relative_to(REPO_ROOT)): sha256(path) for path in
                    [*REEF_SHAPE_SOURCES, RME_ID_LIST, DEFAULT_SURVIVAL,
                     args.domain / "site_to_reef.csv", *matrix_paths]},
        "outputs": {key: sha256(path) for key, path in outputs.items()},
        "interpretation_limit": (
            "Missing sources have no evidence-backed habitable area or survival. "
            "Raw connectivity fractions cannot be converted to recruitment fractions. "
            "Historical survival extrema are sensitivities, not 2018–2022 estimates."
        ),
    }
    (out_dir / "metadata.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: metadata[key] for key in (
        "unmatched_external_source_count", "missing_fraction_mean",
        "missing_fraction_max", "missing_source_nonzero_rows")}, indent=2))
    print(f"Wrote {out_dir}")


if __name__ == "__main__":
    main()
