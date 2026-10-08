"""Build the versioned Lizard_Historical_v0.2 domain inputs.

This builder deliberately leaves Lizard_Historical_v0.1 untouched.  It treats the
V2 ``indexV2`` field as the authoritative site order, creates stable site labels,
repairs invalid geometries deterministically, and labels the otherwise-unlabelled
distance workbook in that exact order.

The source workbook is streamed from its XML representation because loading the
8.4 million-cell sheet with a general spreadsheet object model is unnecessarily
slow and memory intensive.  No source file is modified.
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import hashlib
import json
import math
import re
import shutil
import zipfile
from pathlib import Path
from typing import Iterable
from xml.etree import ElementTree as ET

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from scipy.io import netcdf_file
from scipy.stats import lognorm
from shapely import make_valid
from shapely.geometry import GeometryCollection, MultiPolygon, Polygon
from shapely.ops import unary_union


REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_ROOT = REPO_ROOT / "sandbox" / "data"
V2_VECTOR_ROOT = DATA_ROOT / "OwenHydro" / "farNorthGbrSubreefs" / "farNorthGbrSubreefsV2"
V2_VECTOR = V2_VECTOR_ROOT / "farNorthGbrSubreefsV2LL.shp"
CONNECTIVITY_WORKBOOK = DATA_ROOT / "OwenHydro" / "conMatFngbrDist 1.xlsx"
RME_ROOT = DATA_ROOT / "rme_ml_2025_06_05" / "data_files"
RME_ID_LIST = RME_ROOT / "id" / "id_list_2024_12_01.csv"
RME_INITIAL_ROOT = RME_ROOT / "initial_csv"
HISTORICAL_DHW = DATA_ROOT / "dhw_historical.csv"
COTS_RASTER = DATA_ROOT / "COTS_prob_0.02_cpue_year2025_clean.tif"
GLOBAL_REEF_SHAPE_ROOT = DATA_ROOT / "OwenHydro" / "connectivityMatsYearlyGlobal" / "reefShapefiles"
GLOBAL_REEF_SHAPE = GLOBAL_REEF_SHAPE_ROOT / "gbrShapeLL.shp"
DEFAULT_OUTPUT = DATA_ROOT / "Lizard_Historical_v0.2"

SITE_CRS = "EPSG:4326"
AREA_CRS = "EPSG:32755"
SITE_ID_PREFIX = "FNG_V2_"
EXPECTED_SITE_COUNT = 2_895
EXPECTED_REEF_COUNT = 113
HISTORICAL_YEARS = np.arange(1985, 2025, dtype=np.int32)


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


def polygonal_geometry(geometry):
    """Return a valid polygonal geometry without changing an already-valid input."""
    repaired = geometry if geometry.is_valid else make_valid(geometry)
    if isinstance(repaired, (Polygon, MultiPolygon)):
        return repaired
    if isinstance(repaired, GeometryCollection):
        polygons = [part for part in repaired.geoms if isinstance(part, (Polygon, MultiPolygon))]
        if polygons:
            merged = unary_union(polygons)
            if isinstance(merged, (Polygon, MultiPolygon)):
                return merged
    raise ValueError(f"Geometry repair did not yield a polygon: {repaired.geom_type}")


def sample_cots_density(points: gpd.GeoSeries, raster_path: Path) -> tuple[np.ndarray, dict[str, object]]:
    with rasterio.open(raster_path) as dataset:
        raster_points = points.to_crs(dataset.crs)
        coordinates = [(point.x, point.y) for point in raster_points]
        values = np.fromiter(
            (float(sample[0]) for sample in dataset.sample(coordinates, indexes=1)),
            dtype=np.float64,
            count=len(coordinates),
        )
        invalid = ~np.isfinite(values)
        if dataset.nodata is not None:
            invalid |= np.isclose(values, float(dataset.nodata))
        invalid_count = int(invalid.sum())
        values[invalid] = 0.0
        if np.any(values < 0.0):
            raise ValueError("COTS raster produced negative centroid values")
        details = {
            "method": "nearest raster cell sampled at projected polygon centroid",
            "raster_crs": str(dataset.crs),
            "nodata_or_nonfinite_replaced_with_zero": invalid_count,
            "minimum": float(values.min()),
            "maximum": float(values.max()),
        }
        return values, details


def build_spatial(output_root: Path, skip_cots_density: bool) -> tuple[gpd.GeoDataFrame, dict[str, object]]:
    source = gpd.read_file(V2_VECTOR)
    if source.crs is None or source.crs.to_epsg() != 4326:
        raise ValueError(f"Expected V2 polygons in EPSG:4326, found {source.crs}")
    if len(source) != EXPECTED_SITE_COUNT:
        raise ValueError(f"Expected {EXPECTED_SITE_COUNT} V2 polygons, found {len(source)}")
    if "indexV2" not in source or source["indexV2"].isna().any():
        raise ValueError("V2 polygons require a complete indexV2 field")

    sites = source.sort_values("indexV2", kind="stable").reset_index(drop=True).copy()
    expected_order = np.arange(1, EXPECTED_SITE_COUNT + 1)
    if not np.array_equal(sites["indexV2"].to_numpy(dtype=np.int64), expected_order):
        raise ValueError("indexV2 must be unique and exactly equal to 1:2895")
    if sites["UNIQUE_ID"].isna().any() or sites["GBRMPA_ID"].isna().any():
        raise ValueError("Every V2 site must have UNIQUE_ID and GBRMPA_ID parent identifiers")
    if sites["GBRMPA_ID"].astype(str).nunique() != EXPECTED_REEF_COUNT:
        raise ValueError(f"Expected {EXPECTED_REEF_COUNT} parent reefs")

    original_valid = sites.geometry.is_valid.to_numpy()
    sites.geometry = sites.geometry.map(polygonal_geometry)
    if not bool(sites.geometry.is_valid.all()):
        raise ValueError("Some geometries remain invalid after deterministic make_valid repair")

    projected = sites.to_crs(AREA_CRS)
    area_m2 = projected.geometry.area.to_numpy(dtype=np.float64)
    if not np.all(np.isfinite(area_m2)) or np.any(area_m2 <= 0.0):
        raise ValueError("Projected site areas must be positive and finite")
    projected_centroids = projected.geometry.centroid
    lonlat_centroids = gpd.GeoSeries(projected_centroids, crs=AREA_CRS).to_crs(SITE_CRS)

    source_site_id = sites["site_id"].astype(str).copy()
    canonical_site_id = [f"{SITE_ID_PREFIX}{idx:04d}" for idx in expected_order]
    if len(set(canonical_site_id)) != EXPECTED_SITE_COUNT:
        raise AssertionError("Canonical V2 site IDs are not unique")

    if "fid" in sites:
        sites = sites.rename(columns={"fid": "source_fid"})
    sites = sites.rename(columns={"area": "source_area_native"})
    sites["source_site_id"] = source_site_id
    sites["site_id"] = canonical_site_id
    sites["UNIQUE_ID"] = sites["UNIQUE_ID"].astype(str)
    sites["GBRMPA_ID"] = sites["GBRMPA_ID"].astype(str)
    reef_names = gpd.read_file(GLOBAL_REEF_SHAPE, columns=["LABEL_ID", "reefName"])
    if bool(reef_names["LABEL_ID"].astype(str).duplicated().any()):
        raise ValueError("Global reef shapefile contains duplicate LABEL_ID values")
    reef_name_lookup = dict(
        zip(reef_names["LABEL_ID"].astype(str), reef_names["reefName"].astype(str))
    )
    missing_reef_names = sorted(set(sites["GBRMPA_ID"]) - set(reef_name_lookup))
    if missing_reef_names:
        raise ValueError(f"V2 reefs missing global reef names: {missing_reef_names}")
    sites["reef_name"] = sites["GBRMPA_ID"].map(reef_name_lookup)
    sites["area"] = area_m2
    sites["x_coord"] = lonlat_centroids.x.to_numpy(dtype=np.float64)
    sites["y_coord"] = lonlat_centroids.y.to_numpy(dtype=np.float64)
    sites["k"] = 0.5
    sites["depth_med"] = 5.0
    sites["recs_per_m2"] = 10.0
    sites["geometry_repaired"] = (~original_valid).astype(np.int8)

    cots_details: dict[str, object]
    if skip_cots_density:
        sites["cots_density"] = 0.0
        cots_details = {"method": "skipped by command-line option; all values set to zero"}
    else:
        values, cots_details = sample_cots_density(lonlat_centroids, COTS_RASTER)
        sites["cots_density"] = values

    preferred = [
        "site_id",
        "source_site_id",
        "indexV2",
        "indexV1",
        "UNIQUE_ID",
        "GBRMPA_ID",
        "reef_name",
        "area",
        "source_area_native",
        "x_coord",
        "y_coord",
        "k",
        "depth_med",
        "recs_per_m2",
        "cots_density",
        "geometry_repaired",
    ]
    remaining = [column for column in sites.columns if column not in preferred + ["geometry"]]
    sites = sites[preferred + remaining + ["geometry"]]

    spatial_dir = output_root / "spatial"
    spatial_dir.mkdir(parents=True, exist_ok=True)
    spatial_path = spatial_dir / "lizard_cluster.gpkg"
    if spatial_path.exists():
        spatial_path.unlink()
    sites.to_file(spatial_path, layer="lizard_cluster", driver="GPKG", index=False)

    site_map = pd.DataFrame(
        {
            "site_id": canonical_site_id,
            "source_site_id": source_site_id,
            "indexV2": expected_order,
            "UNIQUE_ID": sites["UNIQUE_ID"],
            "reef_id": sites["GBRMPA_ID"],
            "reef_name": sites["reef_name"],
            "area_m2": area_m2,
        }
    )
    site_map.to_csv(output_root / "site_to_reef.csv", index=False, lineterminator="\n")

    details = {
        "site_count": int(len(sites)),
        "parent_reef_count": int(sites["GBRMPA_ID"].nunique()),
        "site_id_definition": "FNG_V2_ followed by zero-padded indexV2",
        "order_definition": "ascending indexV2; asserted equal to 1:2895",
        "source_site_id_unique_count": int(source_site_id.nunique()),
        "invalid_geometry_count": int((~original_valid).sum()),
        "geometry_repair": "Shapely make_valid; retain or union polygonal components only",
        "area_units": "m^2",
        "area_crs": AREA_CRS,
        "total_area_m2": float(area_m2.sum()),
        "output_crs": SITE_CRS,
        "cots_density": cots_details,
    }
    return sites, details


def workbook_sheet_path(archive: zipfile.ZipFile, requested_name: str) -> str:
    ns = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    rel_ns = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
    workbook = ET.fromstring(archive.read("xl/workbook.xml"))
    rel_id = None
    for sheet in workbook.findall("m:sheets/m:sheet", ns):
        if sheet.attrib.get("name") == requested_name:
            rel_id = sheet.attrib[f"{{{rel_ns}}}id"]
            break
    if rel_id is None:
        raise ValueError(f"Workbook does not contain sheet {requested_name!r}")

    package_ns = {"r": "http://schemas.openxmlformats.org/package/2006/relationships"}
    rels = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
    for rel in rels.findall("r:Relationship", package_ns):
        if rel.attrib["Id"] == rel_id:
            target = rel.attrib["Target"].replace(chr(92), "/")
            return target if target.startswith("xl/") else f"xl/{target.lstrip('/')}"
    raise ValueError(f"Cannot resolve worksheet relationship {rel_id}")


CELL_COLUMN = re.compile(r"([A-Z]+)")


def excel_column_index(cell_reference: str) -> int:
    match = CELL_COLUMN.match(cell_reference)
    if match is None:
        raise ValueError(f"Invalid Excel cell reference {cell_reference!r}")
    value = 0
    for character in match.group(1):
        value = value * 26 + ord(character) - ord("A") + 1
    return value - 1


def worksheet_rows(stream, expected_columns: int) -> Iterable[list[str]]:
    cell_tag = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}c"
    value_tag = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}v"
    row_tag = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}row"
    for _, element in ET.iterparse(stream, events=("end",)):
        if element.tag != row_tag:
            continue
        row = [""] * expected_columns
        for cell in element.findall(cell_tag):
            if cell.attrib.get("t") not in (None, "n"):
                raise ValueError(f"Unexpected non-numeric cell type {cell.attrib.get('t')!r}")
            value_node = cell.find(value_tag)
            if value_node is None or value_node.text is None:
                raise ValueError(f"Blank cell encountered at {cell.attrib.get('r')}")
            column = excel_column_index(cell.attrib["r"])
            if column >= expected_columns:
                raise ValueError(f"Unexpected column at {cell.attrib['r']}")
            row[column] = value_node.text
        if any(value == "" for value in row):
            raise ValueError(f"Worksheet row {element.attrib.get('r')} is incomplete")
        yield row
        element.clear()


def build_connectivity(output_root: Path, site_ids: list[str]) -> dict[str, object]:
    expected = len(site_ids)
    connectivity_dir = output_root / "connectivity"
    connectivity_dir.mkdir(parents=True, exist_ok=True)
    output_path = connectivity_dir / "Lizard_Connectivity.csv"
    temporary_output = output_path.with_suffix(".csv.tmp")

    matrix_path = output_root / ".connectivity_matrix.float64.tmp"
    with contextlib.nullcontext():
        matrix = np.memmap(
            matrix_path,
            dtype=np.float64,
            mode="w+",
            shape=(expected, expected),
        )
        with zipfile.ZipFile(CONNECTIVITY_WORKBOOK) as archive:
            sheet_path = workbook_sheet_path(archive, "in")
            with archive.open(sheet_path) as worksheet, temporary_output.open(
                "w", newline="", encoding="utf-8"
            ) as output_handle:
                writer = csv.writer(output_handle, lineterminator="\n")
                writer.writerow(["Source", *site_ids])
                rows_seen = 0
                for row_index, raw_values in enumerate(worksheet_rows(worksheet, expected)):
                    if row_index >= expected:
                        raise ValueError(f"Workbook has more than {expected} data rows")
                    numeric = np.asarray(raw_values, dtype=np.float64)
                    if not np.all(np.isfinite(numeric)):
                        raise ValueError(f"Connectivity row {row_index + 1} contains non-finite values")
                    if np.any(numeric < 0.0):
                        raise ValueError(f"Connectivity row {row_index + 1} contains negative values")
                    matrix[row_index, :] = numeric
                    writer.writerow([site_ids[row_index], *raw_values])
                    rows_seen += 1
                    if rows_seen % 250 == 0:
                        print(f"  converted {rows_seen}/{expected} connectivity rows", flush=True)
        if rows_seen != expected:
            raise ValueError(f"Workbook has {rows_seen} rows; expected {expected}")
        matrix.flush()

        max_asymmetry = 0.0
        for start in range(0, expected, 256):
            stop = min(start + 256, expected)
            max_asymmetry = max(
                max_asymmetry,
                float(np.max(np.abs(matrix[start:stop, :] - matrix[:, start:stop].T))),
            )
        if max_asymmetry > 1e-12:
            raise ValueError(f"Interim matrix is not symmetric; max difference {max_asymmetry}")
        row_sums = np.asarray(matrix.sum(axis=1))
        details = {
            "dimensions": [expected, expected],
            "orientation": "rows are sources; columns are sinks",
            "label_and_order_source": "ascending V2 indexV2 (1:2895)",
            "minimum": float(matrix.min()),
            "maximum": float(matrix.max()),
            "row_sum_minimum": float(row_sums.min()),
            "row_sum_maximum": float(row_sums.max()),
            "maximum_absolute_asymmetry": max_asymmetry,
            "workbook_sheet": "in",
        }
        del matrix
        matrix_path.unlink()

    shutil.move(temporary_output, output_path)
    return details


def rme_initial_cover_by_reef() -> tuple[dict[str, np.ndarray], list[Path]]:
    files = [RME_INITIAL_ROOT / f"coral_sp{species}_2023.csv" for species in range(2, 7)]
    frames = [pd.read_csv(path, comment="#", header=None, dtype={0: str}) for path in files]
    reference_ids = frames[0].iloc[:, 0].astype(str).to_numpy()
    if any(not np.array_equal(frame.iloc[:, 0].astype(str).to_numpy(), reference_ids) for frame in frames[1:]):
        raise ValueError("RME coral initial-cover files do not share an identical reef order")
    group_cover = np.column_stack(
        [frame.iloc[:, 1:].astype(float).mean(axis=1).to_numpy() / 100.0 for frame in frames]
    )

    id_list = pd.read_csv(RME_ID_LIST, comment="#", header=None, dtype={0: str})
    rme_k = dict(zip(id_list.iloc[:, 0].astype(str), 1.0 - id_list.iloc[:, 2].astype(float)))

    diameter_edges_cm = np.asarray(
        [
            [2.5, 7.5, 12.5, 25.0, 50.0, 80.0, 120.0, 160.0],
            [2.5, 7.5, 12.5, 20.0, 30.0, 60.0, 100.0, 150.0],
            [2.5, 7.5, 12.5, 20.0, 30.0, 40.0, 50.0, 60.0],
            [2.5, 5.0, 7.5, 10.0, 20.0, 40.0, 50.0, 100.0],
            [2.5, 5.0, 7.5, 10.0, 20.0, 40.0, 50.0, 100.0],
        ]
    )
    edge_area_cm2 = math.pi * (diameter_edges_cm / 2.0) ** 2
    cdf = lognorm.cdf(edge_area_cm2, s=math.log(4.0), scale=700.0)
    weights = np.diff(cdf, axis=1)
    weights /= weights.sum(axis=1, keepdims=True)

    output: dict[str, np.ndarray] = {}
    for row, reef_id in enumerate(reference_ids):
        if reef_id not in rme_k or not (0.0 < rme_k[reef_id] <= 1.0):
            raise ValueError(f"Invalid RME habitable fraction for {reef_id}")
        # This reproduces the existing Lizard builder's RME-relative state semantics.
        output[reef_id] = (group_cover[row, :, None] * weights / rme_k[reef_id]).reshape(35)
    return output, files


def write_initial_cover(output_root: Path, reef_ids: list[str]) -> dict[str, object]:
    reef_cover, _ = rme_initial_cover_by_reef()
    missing = sorted(set(reef_ids) - set(reef_cover))
    if missing:
        raise ValueError(f"V2 reefs absent from RME initial-cover files: {missing}")
    cover = np.column_stack([reef_cover[reef_id] for reef_id in reef_ids]).astype(np.float64)
    if cover.shape != (35, len(reef_ids)) or np.any(cover < 0.0) or np.any(cover.sum(axis=0) > 1.0):
        raise ValueError("Derived initial coral cover is invalid")

    path = output_root / "initial_cover.nc"
    with netcdf_file(path, "w") as dataset:
        dataset.createDimension("species", cover.shape[0])
        dataset.createDimension("location", cover.shape[1])
        dataset.createVariable("species", "i4", ("species",))[:] = np.arange(1, cover.shape[0] + 1)
        dataset.createVariable("location", "i4", ("location",))[:] = np.arange(1, cover.shape[1] + 1)
        # NetCDF.jl/YAXArrays exposes dimensions in reverse storage order for this
        # classic NetCDF file, so store location x species to load as species x site.
        variable = dataset.createVariable("covers", "f8", ("location", "species"))
        variable[:] = cover.T
        variable.units = "fraction of ReefMod habitable area"
        dataset.history = "Derived by build_lizard_domain_v2.py from mean ReefMod 2023 repeats"
    return {
        "dimensions": list(cover.shape),
        "minimum_total_cover": float(cover.sum(axis=0).min()),
        "maximum_total_cover": float(cover.sum(axis=0).max()),
        "semantics": "RME absolute cover divided by each source reef's RME habitable fraction, matching the v0.1 extraction path",
        "known_limitation": "spatial k remains 0.5; reconciling cover to the Lizard k definition is a separate mechanism change",
    }


def normalized_unique_id(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    if re.fullmatch(r"\d+\.0", text):
        return text[:-2]
    return text


def write_historical_dhw(output_root: Path, unique_ids: list[str]) -> dict[str, object]:
    source = pd.read_csv(HISTORICAL_DHW, dtype={"RME_UNIQUE_ID": str})
    source["RME_UNIQUE_ID"] = source["RME_UNIQUE_ID"].map(normalized_unique_id)
    source["timestep"] = source["timestep"].astype(int)
    source["dhw"] = source["dhw"].astype(float)
    duplicate_mask = source.duplicated(["RME_UNIQUE_ID", "timestep"], keep=False)
    duplicate_row_count = int(duplicate_mask.sum())
    if duplicate_row_count:
        duplicate_groups = source.loc[duplicate_mask].groupby(
            ["RME_UNIQUE_ID", "timestep"], sort=False
        )["dhw"]
        conflicting = duplicate_groups.nunique().gt(1)
        if bool(conflicting.any()):
            bad_keys = list(conflicting[conflicting].index[:10])
            raise ValueError(f"Historical DHW has conflicting duplicate keys: {bad_keys}")
        source = source.drop_duplicates(["RME_UNIQUE_ID", "timestep"], keep="first")
    lookup = source.set_index(["RME_UNIQUE_ID", "timestep"])["dhw"]
    missing = [
        (reef_id, int(year))
        for reef_id in sorted(set(unique_ids))
        for year in HISTORICAL_YEARS
        if (reef_id, int(year)) not in lookup.index
    ]
    if missing:
        raise ValueError(f"Historical DHW is missing V2 reef-years; first entries: {missing[:10]}")

    dhw = np.empty((len(HISTORICAL_YEARS), len(unique_ids), 2), dtype=np.float32)
    for site_index, reef_id in enumerate(unique_ids):
        values = np.asarray([lookup.loc[(reef_id, int(year))] for year in HISTORICAL_YEARS], dtype=np.float32)
        dhw[:, site_index, 0] = values
        dhw[:, site_index, 1] = values
    if not np.all(np.isfinite(dhw)) or np.any(dhw < 0.0):
        raise ValueError("Derived historical DHW cube is invalid")

    dhw_dir = output_root / "DHWs"
    dhw_dir.mkdir(parents=True, exist_ok=True)
    for name in ("dhw_RCPhistorical.nc", "dhw_RCP45.nc"):
        path = dhw_dir / name
        with netcdf_file(path, "w") as dataset:
            dataset.createDimension("time", dhw.shape[0])
            dataset.createDimension("location", dhw.shape[1])
            dataset.createDimension("scenario", dhw.shape[2])
            dataset.createVariable("time", "i4", ("time",))[:] = HISTORICAL_YEARS
            dataset.createVariable("location", "i4", ("location",))[:] = np.arange(1, dhw.shape[1] + 1)
            dataset.createVariable("scenario", "i4", ("scenario",))[:] = np.arange(1, dhw.shape[2] + 1)
            # Store in reverse order so YAXArrays loads time x location x scenario,
            # matching the historical Lizard loader contract.
            variable = dataset.createVariable("dhw", "f4", ("scenario", "location", "time"))
            variable[:] = np.transpose(dhw, (2, 1, 0))
            variable.units = "degree heating weeks"
            dataset.history = "Historical 1985-2024 reef values expanded to V2 sites"
    return {
        "dimensions": list(dhw.shape),
        "calendar_years": [int(HISTORICAL_YEARS[0]), int(HISTORICAL_YEARS[-1])],
        "minimum": float(dhw.min()),
        "maximum": float(dhw.max()),
        "nonzero_count": int(np.count_nonzero(dhw)),
        "identical_duplicate_source_rows_removed": duplicate_row_count // 2,
        "RCP45_alias": "contains the same historical cube for compatibility; it is not a future RCP projection",
    }


def write_datapackage(output_root: Path) -> None:
    metadata = {
        "name": "Lizard_Historical",
        "version": "0.8.0",
        "domain_data_version": "0.2.0",
        "description": "Versioned Far North GBR/Lizard domain using V2 subreef polygons and indexV2-labelled interim site connectivity",
        "spatial": "spatial/lizard_cluster.gpkg",
        "connectivity": "connectivity/Lizard_Connectivity.csv",
        "initial_cover": "initial_cover.nc",
        "dhw": "DHWs/dhw_RCPhistorical.nc",
    }
    (output_root / "datapackage.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--skip-cots-density",
        action="store_true",
        help="Set site COTS density to zero instead of sampling the source raster",
    )
    arguments = parser.parse_args()
    output_root = arguments.output.resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    for directory in ("waves", "cyclones", "cots_connectivity"):
        (output_root / directory).mkdir(exist_ok=True)

    print(f"Building spatial layer in {output_root}", flush=True)
    sites, spatial_details = build_spatial(output_root, arguments.skip_cots_density)
    site_ids = sites["site_id"].astype(str).tolist()
    reef_ids = sites["GBRMPA_ID"].astype(str).tolist()
    unique_ids = sites["UNIQUE_ID"].astype(str).tolist()

    print("Streaming and labelling interim site connectivity", flush=True)
    connectivity_details = build_connectivity(output_root, site_ids)
    print("Deriving initial coral cover from ReefMod source tables", flush=True)
    cover_details = write_initial_cover(output_root, reef_ids)
    print("Deriving corrected historical DHW cubes", flush=True)
    dhw_details = write_historical_dhw(output_root, unique_ids)
    write_datapackage(output_root)

    vector_sources = sorted(
        path for path in V2_VECTOR_ROOT.glob("farNorthGbrSubreefsV2LL.*") if path.is_file()
    )
    reef_name_sources = sorted(
        path for path in GLOBAL_REEF_SHAPE_ROOT.glob("gbrShapeLL.*") if path.is_file()
    )
    coral_sources = [RME_INITIAL_ROOT / f"coral_sp{species}_2023.csv" for species in range(2, 7)]
    sources = vector_sources + reef_name_sources + [
        CONNECTIVITY_WORKBOOK,
        RME_ID_LIST,
        *coral_sources,
        HISTORICAL_DHW,
    ]
    if not arguments.skip_cots_density:
        sources.append(COTS_RASTER)

    output_paths = [
        output_root / "spatial" / "lizard_cluster.gpkg",
        output_root / "connectivity" / "Lizard_Connectivity.csv",
        output_root / "site_to_reef.csv",
        output_root / "initial_cover.nc",
        output_root / "DHWs" / "dhw_RCPhistorical.nc",
        output_root / "DHWs" / "dhw_RCP45.nc",
        output_root / "datapackage.json",
    ]
    provenance = {
        "schema_version": 1,
        "domain_version": "Lizard_Historical_v0.2",
        "builder": Path(__file__).relative_to(REPO_ROOT).as_posix(),
        "builder_sha256": sha256_file(Path(__file__)),
        "baseline_preserved": "sandbox/data/Lizard_Historical_v0.1 is not modified",
        "spatial": spatial_details,
        "coral_connectivity": connectivity_details,
        "initial_cover": cover_details,
        "dhw": dhw_details,
        "sources": [source_record(path) for path in sources],
        "outputs": [
            {
                "path": path.relative_to(output_root).as_posix(),
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
            }
            for path in output_paths
        ],
    }
    provenance_path = output_root / "provenance.json"
    provenance_path.write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"Built {output_root}")
    print(json.dumps({"spatial": spatial_details, "connectivity": connectivity_details, "dhw": dhw_details}, indent=2))


if __name__ == "__main__":
    main()
