"""Build the experimental Owen 2018-2023 COTS forcing for Lizard v0.2.

The dataset is written below ``cots_connectivity/datasets`` and is never selected
unless ``ADRIA_COTS_CONNECTIVITY_DATASET=owen_global_2018_2023`` is set.  Default
policies fail closed.  Provisional policies must be supplied explicitly and are
recorded in provenance so they cannot be mistaken for evidence-backed defaults.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_ROOT = REPO_ROOT / "sandbox" / "data"
DEFAULT_DOMAIN = DATA_ROOT / "Lizard_Historical_v0.2"
SOURCE_ROOT = DATA_ROOT / "OwenHydro" / "connectivityMatsYearlyGlobal"
REEF_SHAPE = SOURCE_ROOT / "reefShapefiles" / "gbrShapeLL.shp"
REEF_SHAPE_SOURCES = sorted(REEF_SHAPE.parent.glob("gbrShapeLL.*"))
PREDICTIONS = DATA_ROOT / "gbrPredsAdj_20262408.csv"
RME_ROOT = DATA_ROOT / "rme_ml_2025_06_05" / "data_files"
RME_ID_LIST = RME_ROOT / "id" / "id_list_2024_12_01.csv"
DEFAULT_SURVIVAL = RME_ROOT / "water_csv" / "COTS_LARVAL_REDUCTION_q3baseline.csv"
DATASET_NAME = "owen_global_2018_2023"
OUTBREAK_CPUE_THRESHOLD = 0.22

NEW_TO_OLD_ALIASES = {
    "11-325": "10-441",
    "11-244e": "11-288",
    "11-244f": "11-303",
    "11-244g": "11-310",
    "11-244h": "11-311",
    "20198": "20-198",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def source_record(path: Path) -> dict[str, object]:
    return {
        "path": path.relative_to(REPO_ROOT).as_posix(),
        "sha256": sha256_file(path),
        "bytes": path.stat().st_size,
    }


def reef_id_from_name(name: str) -> str:
    match = re.search(r"\(([^()]*)\)\s*$", str(name))
    if match is None:
        raise ValueError(f"Cannot extract reef ID from {name!r}")
    return match.group(1).strip()


def matrix_year(path: Path) -> int:
    match = re.fullmatch(r"conMatCots(\d{4})-\d{2}\.csv", path.name)
    if match is None:
        raise ValueError(f"Cannot extract hydrodynamic year from {path.name}")
    return int(match.group(1))


def write_matrix_csv(path: Path, matrix: np.ndarray, labels: list[str], first_column: str) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow([first_column, *labels])
        for label, row in zip(labels, matrix):
            writer.writerow([label, *(format(float(value), ".17g") for value in row)])


def write_static_site_matrix(
    path: Path,
    reef_matrix: np.ndarray,
    site_ids: list[str],
    site_to_reef: np.ndarray,
) -> float:
    reef_site_counts = np.bincount(site_to_reef, minlength=reef_matrix.shape[0])
    if np.any(reef_site_counts == 0):
        raise ValueError("Every selected reef must contain at least one V2 site")
    sink_divisors = reef_site_counts[site_to_reef]
    max_mass_error = 0.0
    for source_reef in range(reef_matrix.shape[0]):
        site_row = reef_matrix[source_reef, site_to_reef] / sink_divisors
        recovered = np.bincount(
            site_to_reef, weights=site_row, minlength=reef_matrix.shape[0]
        )
        max_mass_error = max(
            max_mass_error,
            float(np.max(np.abs(recovered - reef_matrix[source_reef, :]))),
        )
    if max_mass_error > 1e-12:
        raise ValueError(f"Reef-to-site expansion failed mass check: {max_mass_error}")
    temporary_path = path.with_suffix(".csv.tmp")
    with temporary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["Source", *site_ids])
        for site_index, (site_id, source_reef) in enumerate(zip(site_ids, site_to_reef), start=1):
            row = reef_matrix[source_reef, site_to_reef] / sink_divisors
            writer.writerow([site_id, *(format(float(value), ".17g") for value in row)])
            if site_index % 250 == 0:
                print(f"  expanded {site_index}/{len(site_ids)} COTS site rows", flush=True)
    temporary_path.replace(path)
    return max_mass_error


def quote_toml(value: str) -> str:
    return json.dumps(value, ensure_ascii=False)


def write_forcing_toml(
    path: Path,
    years: list[int],
    matrix_files: list[str],
    policies: dict[str, str],
) -> None:
    lines = [
        "schema_version = 3",
        f"dataset = {quote_toml(DATASET_NAME)}",
        'promotion_status = "experimental_opt_in"',
        'orientation = "reef connectivity rows are sources and columns are sinks"',
        'state_density_units = "COTS ha^-1"',
        'external_source_units = "pre-dispersal recruits ha^-1 yr^-1"',
        f"years = [{', '.join(str(year) for year in years)}]",
        'site_map = "site_reef_map.csv"',
        'source_survival = "source_larval_survival.csv"',
        'external_supply = "external_supply_coefficients.csv"',
        'external_outbreak_boundary = "external_outbreak_boundary.csv"',
        'outbreak_prediction_mapping = "outbreak_prediction_mapping.csv"',
        f"matrix_files = [{', '.join(quote_toml(name) for name in matrix_files)}]",
        f"nonfinite_source_policy = {quote_toml(policies['nonfinite_source_policy'])}",
        f"missing_source_policy = {quote_toml(policies['missing_source_policy'])}",
        f"survival_year_policy = {quote_toml(policies['survival_year_policy'])}",
        'outbreak_boundary_units = "dimensionless multiplier of pre-pelagic outbreak production density in COTS ha^-1 yr^-1"',
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def load_prediction_matrix(global_ids: list[str]) -> tuple[np.ndarray, list[int], pd.DataFrame]:
    predictions = pd.read_csv(PREDICTIONS)
    predictions["reef_id"] = predictions["reefName"].map(reef_id_from_name)
    if bool(predictions.duplicated(["reef_id", "year"]).any()):
        raise ValueError("Outbreak predictions contain duplicate reef-year keys")
    years = sorted(predictions["year"].astype(int).unique().tolist())
    pivot = predictions.pivot(index="reef_id", columns="year", values="outbrProb")
    missing_ids = sorted(set(global_ids) - set(pivot.index.astype(str)))
    extra_ids = sorted(set(pivot.index.astype(str)) - set(global_ids))
    if missing_ids or extra_ids:
        raise ValueError(
            f"Prediction/global reef mismatch: missing={missing_ids[:10]}, extra={extra_ids[:10]}"
        )
    matrix = pivot.loc[global_ids, years].to_numpy(dtype=np.float64)
    if not np.all(np.isfinite(matrix)) or np.any((matrix < 0.0) | (matrix > 1.0)):
        raise ValueError("Outbreak probabilities must be complete and in [0, 1]")
    mapping = pd.DataFrame(
        {
            "connectivity_reef_id": global_ids,
            "prediction_reef_id": global_ids,
            "mapping_status": "exact",
        }
    )
    return matrix, years, mapping


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    parser.add_argument(
        "--nonfinite-source-policy",
        choices=("error", "zero_full_rows"),
        default="error",
    )
    parser.add_argument(
        "--missing-source-policy",
        choices=("error", "exclude"),
        default="error",
    )
    parser.add_argument(
        "--survival-year-policy",
        choices=("error", "latest_available"),
        default="error",
    )
    arguments = parser.parse_args()
    domain_root = arguments.domain.resolve()
    site_map_path = domain_root / "site_to_reef.csv"
    if not site_map_path.is_file():
        raise ValueError(f"Missing V2 site map: {site_map_path}")

    reef_shape = gpd.read_file(REEF_SHAPE)
    if len(reef_shape) != 3_861:
        raise ValueError(f"Expected 3861 global reefs, found {len(reef_shape)}")
    global_ids = [reef_id_from_name(name) for name in reef_shape["reefName"]]
    if len(set(global_ids)) != len(global_ids):
        raise ValueError("Global reef-shape order contains duplicate reef IDs")
    global_lookup = {reef_id: index for index, reef_id in enumerate(global_ids)}

    sites = pd.read_csv(site_map_path, dtype=str)
    sites["indexV2"] = sites["indexV2"].astype(int)
    if not np.array_equal(sites["indexV2"].to_numpy(), np.arange(1, len(sites) + 1)):
        raise ValueError("V2 site map is not in strict indexV2 order")
    site_ids = sites["site_id"].astype(str).tolist()
    local_global_indices = sorted({global_lookup[reef_id] for reef_id in sites["reef_id"]})
    if len(local_global_indices) != 113:
        raise ValueError(f"Expected 113 V2 parent reefs, found {len(local_global_indices)}")
    local_ids = [global_ids[index] for index in local_global_indices]
    local_lookup = {reef_id: index for index, reef_id in enumerate(local_ids)}
    site_to_reef = np.asarray([local_lookup[reef_id] for reef_id in sites["reef_id"]], dtype=int)

    id_list = pd.read_csv(RME_ID_LIST, comment="#", header=None, dtype={0: str})
    old_area_ha = {
        str(row.iloc[0]): float(row.iloc[1]) * 100.0 * (1.0 - float(row.iloc[2]))
        for _, row in id_list.iterrows()
    }
    survival = pd.read_csv(DEFAULT_SURVIVAL, dtype={"id": str})
    survival_lookup = {reef_id: row for row, reef_id in enumerate(survival["id"].astype(str))}

    mapped_old_ids = [
        reef_id if reef_id in old_area_ha else NEW_TO_OLD_ALIASES.get(reef_id)
        for reef_id in global_ids
    ]
    missing_global_ids = [
        reef_id for reef_id, old_id in zip(global_ids, mapped_old_ids)
        if old_id is None or old_id not in survival_lookup or old_id not in old_area_ha
    ]
    if missing_global_ids and arguments.missing_source_policy == "error":
        raise ValueError(
            f"{len(missing_global_ids)} Owen reefs lack RME survival/area data; "
            "rerun only with an explicitly justified --missing-source-policy exclude"
        )
    if any(reef_id in missing_global_ids for reef_id in local_ids):
        raise ValueError("A local V2 reef lacks RME survival/area data")

    outbreak_probability, prediction_years, prediction_mapping = load_prediction_matrix(global_ids)
    prediction_mapping["rme_source_id"] = [old_id or "" for old_id in mapped_old_ids]
    prediction_mapping["rme_mapping_status"] = [
        "missing_excluded"
        if reef_id in missing_global_ids
        else ("documented_alias" if reef_id in NEW_TO_OLD_ALIASES else "exact")
        for reef_id in global_ids
    ]

    matrix_paths = sorted(SOURCE_ROOT.glob("conMatCots*.csv"), key=matrix_year)
    years = [matrix_year(path) for path in matrix_paths]
    available_survival_years = sorted(int(column) for column in survival.columns if column != "id")
    survival_year_by_hydrodynamic_year: dict[int, int] = {}
    for year in years:
        if str(year) in survival.columns:
            survival_year_by_hydrodynamic_year[year] = year
        elif arguments.survival_year_policy == "latest_available":
            survival_year_by_hydrodynamic_year[year] = max(available_survival_years)
        else:
            raise ValueError(
                f"No COTS larval-survival column for hydrodynamic year {year}; "
                "rerun only with an explicitly justified --survival-year-policy latest_available"
            )

    dataset_root = domain_root / "cots_connectivity" / "datasets" / DATASET_NAME
    forcing_root = dataset_root / "annual_forcing"
    dataset_root.mkdir(parents=True, exist_ok=True)
    forcing_root.mkdir(parents=True, exist_ok=True)

    local_areas = np.asarray([old_area_ha[mapped_old_ids[index]] for index in local_global_indices])
    site_reef_map = pd.DataFrame(
        {
            "site_id": site_ids,
            "reef_index": site_to_reef + 1,
            "reef_id": [local_ids[index] for index in site_to_reef],
            "reef_habitable_area_ha": local_areas[site_to_reef],
        }
    )
    site_reef_map.to_csv(forcing_root / "site_reef_map.csv", index=False, lineterminator="\n")

    local_survival_columns: dict[int, np.ndarray] = {}
    external_supply_columns: dict[int, np.ndarray] = {}
    internal_matrices: list[np.ndarray] = []
    boundary_rows: list[pd.DataFrame] = []
    nonfinite_report: list[dict[str, object]] = []
    matrix_output_names: list[str] = []
    global_local_mask = np.zeros(len(global_ids), dtype=bool)
    global_local_mask[local_global_indices] = True
    mapped_mask = np.asarray([reef_id not in missing_global_ids for reef_id in global_ids])
    external_mask = ~global_local_mask & mapped_mask

    for source_path, hydrodynamic_year in zip(matrix_paths, years):
        print(f"Reading {source_path.name}", flush=True)
        matrix = np.loadtxt(source_path, delimiter=",", dtype=np.float64)
        if matrix.shape != (len(global_ids), len(global_ids)):
            raise ValueError(f"{source_path.name} has shape {matrix.shape}")
        nonfinite_count = np.sum(~np.isfinite(matrix), axis=1)
        mixed_rows = np.flatnonzero((nonfinite_count > 0) & (nonfinite_count < len(global_ids)))
        if len(mixed_rows):
            raise ValueError(f"{source_path.name} has partially non-finite rows: {mixed_rows[:10] + 1}")
        full_nonfinite_rows = np.flatnonzero(nonfinite_count == len(global_ids))
        if len(full_nonfinite_rows) and arguments.nonfinite_source_policy == "error":
            raise ValueError(
                f"{source_path.name} has {len(full_nonfinite_rows)} wholly non-finite source rows; "
                "rerun only with --nonfinite-source-policy zero_full_rows"
            )
        if len(full_nonfinite_rows):
            matrix[full_nonfinite_rows, :] = 0.0
        if not np.all(np.isfinite(matrix)) or np.any((matrix < 0.0) | (matrix > 1.0)):
            raise ValueError(f"{source_path.name} contains invalid connectivity probabilities")
        nonfinite_report.append(
            {
                "matrix": source_path.name,
                "full_nonfinite_source_count": int(len(full_nonfinite_rows)),
                "full_nonfinite_source_ids": [global_ids[index] for index in full_nonfinite_rows],
            }
        )

        survival_year = survival_year_by_hydrodynamic_year[hydrodynamic_year]
        full_survival = np.zeros(len(global_ids), dtype=np.float64)
        full_area = np.zeros(len(global_ids), dtype=np.float64)
        for index, old_id in enumerate(mapped_old_ids):
            if not mapped_mask[index]:
                continue
            full_survival[index] = float(survival.loc[survival_lookup[old_id], str(survival_year)])
            full_area[index] = old_area_ha[old_id]
        if np.any((full_survival < 0.0) | (full_survival > 1.0)):
            raise ValueError("Larval-survival modifiers must be in [0, 1]")

        internal = matrix[np.ix_(local_global_indices, local_global_indices)]
        internal_matrices.append(internal)
        local_survival_columns[hydrodynamic_year] = full_survival[local_global_indices]

        external_connection = matrix[np.ix_(external_mask, local_global_indices)]
        external_survival = full_survival[external_mask]
        external_supply_columns[hydrodynamic_year] = (
            external_connection * external_survival[:, None]
        ).sum(axis=0)
        pressure = external_connection * (
            full_area[external_mask, None] / local_areas[None, :]
        )
        potential = pressure.T @ outbreak_probability[external_mask, :]
        pelagic = (pressure * external_survival[:, None]).T @ outbreak_probability[external_mask, :]
        if not np.all(np.isfinite(potential)) or not np.all(np.isfinite(pelagic)):
            raise ValueError("External outbreak boundary contains non-finite values")
        boundary_rows.append(
            pd.DataFrame(
                {
                    "reef_id": np.repeat(local_ids, len(prediction_years)),
                    "calendar_year": np.tile(prediction_years, len(local_ids)),
                    "hydrodynamic_year": hydrodynamic_year,
                    "probability_area_connectivity": potential.reshape(-1),
                    "survival_weighted_probability_area_connectivity": pelagic.reshape(-1),
                }
            )
        )

        matrix_name = f"reef_connectivity_{hydrodynamic_year}.csv"
        write_matrix_csv(forcing_root / matrix_name, internal, local_ids, "Source_Reef")
        matrix_output_names.append(matrix_name)
        del matrix

    mean_internal = np.mean(np.stack(internal_matrices, axis=0), axis=0)
    max_site_expansion_mass_error = write_static_site_matrix(
        dataset_root / "Lizard_COTS_Connectivity.csv",
        mean_internal,
        site_ids,
        site_to_reef,
    )

    survival_output = pd.DataFrame({"reef_id": local_ids})
    external_output = pd.DataFrame({"reef_id": local_ids})
    for year in years:
        survival_output[str(year)] = local_survival_columns[year]
        external_output[str(year)] = external_supply_columns[year]
    survival_output.to_csv(
        forcing_root / "source_larval_survival.csv", index=False, lineterminator="\n"
    )
    external_output.to_csv(
        forcing_root / "external_supply_coefficients.csv", index=False, lineterminator="\n"
    )
    prediction_mapping.to_csv(
        forcing_root / "outbreak_prediction_mapping.csv", index=False, lineterminator="\n"
    )
    pd.concat(boundary_rows, ignore_index=True).to_csv(
        forcing_root / "external_outbreak_boundary.csv", index=False, lineterminator="\n"
    )

    policies = {
        "nonfinite_source_policy": arguments.nonfinite_source_policy,
        "missing_source_policy": arguments.missing_source_policy,
        "survival_year_policy": arguments.survival_year_policy,
    }
    write_forcing_toml(
        forcing_root / "provenance.toml", years, matrix_output_names, policies
    )

    output_paths = sorted(
        path for path in dataset_root.rglob("*")
        if path.is_file() and path != dataset_root / "provenance.json"
    )
    provenance = {
        "schema_version": 1,
        "dataset": DATASET_NAME,
        "builder": Path(__file__).relative_to(REPO_ROOT).as_posix(),
        "builder_sha256": sha256_file(Path(__file__)),
        "promotion_status": "experimental_opt_in",
        "selection": {
            "ADRIA_COTS_CONNECTIVITY_MODE": "cots",
            "ADRIA_COTS_CONNECTIVITY_DATASET": DATASET_NAME,
        },
        "orientation": "rows are larval sources; columns are settlement sinks",
        "hydrodynamic_years": years,
        "survival_year_by_hydrodynamic_year": survival_year_by_hydrodynamic_year,
        "policies": policies,
        "missing_rme_survival_or_area_source_count": len(missing_global_ids),
        "missing_rme_survival_or_area_source_ids": missing_global_ids,
        "nonfinite_sources": nonfinite_report,
        "local_reef_count": len(local_ids),
        "site_count": len(site_ids),
        "maximum_reef_to_site_mass_error": max_site_expansion_mass_error,
        "outbreak_probability": {
            "calendar_years": prediction_years,
            "event_definition": "probability COTS CPUE exceeds 0.22 COTS per tow",
            "cpue_threshold_cots_per_tow": OUTBREAK_CPUE_THRESHOLD,
        },
        "known_limitations": [
            "2018-2022 hydrodynamic matrices have no same-year survival columns in the supplied RME water-quality file",
            "56 Owen source reefs have no RME habitable-area or larval-survival match",
            "wholly non-finite source rows are present in the first four matrices",
            "the dataset is not a default and has not passed peak-timing/amplitude promotion gates",
        ],
        "sources": [
            *(source_record(path) for path in REEF_SHAPE_SOURCES),
            source_record(PREDICTIONS),
            source_record(RME_ID_LIST),
            source_record(DEFAULT_SURVIVAL),
            *(source_record(path) for path in matrix_paths),
        ],
        "outputs": [
            {
                "path": path.relative_to(dataset_root).as_posix(),
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
            }
            for path in output_paths
        ],
    }
    (dataset_root / "provenance.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(
        f"Built experimental {DATASET_NAME}: sites={len(site_ids)}, reefs={len(local_ids)}, "
        f"years={years}, missing_sources={len(missing_global_ids)}"
    )


if __name__ == "__main__":
    main()
