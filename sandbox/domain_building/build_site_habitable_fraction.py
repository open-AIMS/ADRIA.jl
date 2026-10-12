"""Derive site-level habitable fractions for the Lizard V2 domain from the
GBR10 2025 benthic habitat map (INTERNAL ONLY source; outputs stay in the
git-ignored sandbox/data tree).

Two habitable-fraction definitions are written:

* ``k_benthic`` (k_total): Rock (13) + Coral/Algae (15) cells divided by all
  cells inside the site polygon; sand (11), rubble (12), unlabelled and unmapped
  cells count as non-habitable.
* ``k_benthic_mapped``: the same habitable cells divided by *mapped* cells,
  assuming unmapped area (typically beyond the optically mapped depth) has the
  mapped composition. Sites with less than ``MIN_MAPPED`` mapped area use the
  area-weighted mapped composition of their reef, then the RME reef value
  (1 - sand proportion) if the reef itself is unmapped.

* ``k_geomorphic_infill``: benthic-mapped cells keep their benthic class; each
  benthic-unmapped cell contributes the habitable share ``p(g)`` of its GBR10
  2025 geomorphic zone ``g``. ``p(g)`` is estimated over all domain-site cells
  where both layers are mapped (habitable / benthic-mapped cells in zone ``g``).
  Cells unmapped in both layers, and zones with fewer than ``MIN_ZONE_CELLS``
  benthic-mapped cells, contribute zero.

Per-class fractions and a site x benthic x geomorphic cell-count table are
written so alternative definitions can be computed without re-reading the
rasters.

Run from the ADRIA.jl repository root (requires rasterio with GDAL >= 3.7 for
FileGDB rasters):

    python sandbox/domain_building/build_site_habitable_fraction.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.features import geometry_mask
from rasterio.windows import Window, from_bounds

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SITES = REPO_ROOT / "sandbox/data/Lizard_Historical_v0.2/spatial/lizard_cluster.gpkg"
DEFAULT_BENTHIC = Path(
    r"C:\Users\smatthew\Documents\GitKraken\SpatialData\GBR-HabitatMapping_INTERNALONLY"
    r"\gbr10_benthic_2025\gbr10_benthic_2025.gdb")
DEFAULT_GEOMORPHIC = Path(
    r"C:\Users\smatthew\Documents\GitKraken\SpatialData\GBR-HabitatMapping_INTERNALONLY"
    r"\gbr10_geomorphic_2025\gbr10_geomorphic_2025.gdb")
DEFAULT_OUTPUT = REPO_ROOT / "sandbox/data/Lizard_Historical_v0.2/derived/site_habitable_fraction_benthic_geomorphic2025.csv"
GEOMORPHIC_CLASSES = {2: "deep_water", 11: "shallow_lagoon", 12: "deep_lagoon",
                      13: "inner_reef_flat", 14: "outer_reef_flat", 15: "reef_crest",
                      21: "sheltered_reef_slope", 22: "reef_slope", 23: "plateau",
                      24: "back_reef_slope"}
MIN_ZONE_CELLS = 1000
CLASSES = {11: "sand", 12: "rubble", 13: "rock", 15: "coral_algae"}
HABITABLE = (13, 15)
MIN_MAPPED = 0.2
RME_ID_LIST = REPO_ROOT / "sandbox/data/rme_ml_2025_06_05/data_files/id/id_list_2024_12_01.csv"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_tree(root: Path) -> str:
    """Order-stable digest over every file name and content in a directory."""
    digest = hashlib.sha256()
    for path in sorted(p for p in root.rglob("*") if p.is_file()):
        digest.update(str(path.relative_to(root)).encode())
        digest.update(sha256_file(path).encode())
    return digest.hexdigest()


def main(sites_path: Path, benthic_path: Path, geomorphic_path: Path, output: Path) -> None:
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite {output}")
    sites = gpd.read_file(sites_path)
    if sites.crs is None or sites.crs.to_epsg() != 4326:
        raise ValueError("Site polygons must be EPSG:4326")
    if sites.site_id.duplicated().any():
        raise ValueError("Duplicate site IDs")
    rows = []
    crosstab = []
    with rasterio.open(benthic_path) as raster, rasterio.open(geomorphic_path) as geomorphic:
        if raster.crs.to_epsg() != 4326 or raster.count != 1:
            raise ValueError("Unexpected benthic raster CRS or band count")
        if geomorphic.transform != raster.transform or geomorphic.shape != raster.shape:
            raise ValueError("Geomorphic and benthic rasters are not on the same grid")
        nodata = raster.nodata
        for site in sites.itertuples(index=False):
            exact = from_bounds(*site.geometry.bounds, transform=raster.transform)
            col0, row0 = int(np.floor(exact.col_off)), int(np.floor(exact.row_off))
            window = Window(col0, row0,
                            int(np.ceil(exact.col_off + exact.width)) - col0,
                            int(np.ceil(exact.row_off + exact.height)) - row0)
            values = raster.read(1, window=window, boundless=True, fill_value=0)
            zones = geomorphic.read(1, window=window, boundless=True, fill_value=0)
            transform = raster.window_transform(window)
            inside = ~geometry_mask([site.geometry], out_shape=values.shape,
                                    transform=transform, all_touched=False)
            n_inside = int(inside.sum())
            if n_inside == 0:
                raise ValueError(f"No raster cells inside site {site.site_id}")
            classes = values[inside]
            site_zones = zones[inside]
            pairs, pair_counts = np.unique(
                np.stack([classes.astype(np.int64), site_zones.astype(np.int64)]),
                axis=1, return_counts=True)
            for (b, g), n in zip(pairs.T, pair_counts):
                crosstab.append((site.site_id, int(b), int(g), int(n)))
            counts = {name: int((classes == code).sum()) for code, name in CLASSES.items()}
            mapped_other = int(((classes != 0) & ~np.isin(classes, list(CLASSES)) &
                                (classes != (nodata if nodata is not None else -1))).sum())
            unmapped = n_inside - sum(counts.values()) - mapped_other
            row = {"site_id": site.site_id, "GBRMPA_ID": site.GBRMPA_ID,
                   "reef_name": site.reef_name, "n_cells": n_inside}
            for name, count in counts.items():
                row[f"frac_{name}"] = count / n_inside
            row["frac_other_class"] = mapped_other / n_inside
            row["frac_unmapped"] = unmapped / n_inside
            row["k_benthic"] = sum(counts[CLASSES[c]] for c in HABITABLE) / n_inside
            row["k_benthic_with_rubble"] = row["k_benthic"] + row["frac_rubble"]
            rows.append(row)
    table = pd.DataFrame(rows)
    table["area_m2"] = sites.set_index("site_id").loc[table.site_id, "area"].to_numpy()
    table["frac_mapped"] = 1.0 - table.frac_unmapped
    habitable_area = (table.frac_rock + table.frac_coral_algae) * table.area_m2
    mapped_area = table.frac_mapped * table.area_m2
    reef_mapped_k = (habitable_area.groupby(table.GBRMPA_ID).sum() /
                     mapped_area.groupby(table.GBRMPA_ID).sum())
    ids = pd.read_csv(RME_ID_LIST, comment="#", header=None, dtype={0: str})
    rme_k = dict(zip(ids[0].astype(str), 1.0 - ids[2].astype(float)))
    own = habitable_area / mapped_area.where(table.frac_mapped >= MIN_MAPPED)
    reef_fallback = table.GBRMPA_ID.map(reef_mapped_k)
    rme_fallback = table.GBRMPA_ID.map(rme_k)
    table["k_benthic_mapped"] = own.fillna(reef_fallback).fillna(rme_fallback)
    table["k_benthic_mapped_source"] = np.where(own.notna(), "site_mapped",
        np.where(reef_fallback.notna(), "reef_mapped", "rme_reef"))
    if table.k_benthic_mapped.isna().any() or not table.k_benthic_mapped.between(-1e-9, 1 + 1e-9).all():
        raise ValueError("Mapped habitable fraction missing or outside [0, 1]")
    table["k_benthic_mapped"] = table.k_benthic_mapped.clip(0.0, 1.0)  # float round-off only
    cells = pd.DataFrame(crosstab, columns=["site_id", "benthic", "geomorphic", "n_cells"])
    benthic_mapped = cells.benthic != 0
    habitable_cell = cells.benthic.isin(HABITABLE)
    zone = cells[benthic_mapped].groupby("geomorphic").n_cells.sum()
    zone_habitable = cells[benthic_mapped & habitable_cell].groupby("geomorphic").n_cells.sum()
    zone_share = (zone_habitable.reindex(zone.index, fill_value=0) / zone).where(zone >= MIN_ZONE_CELLS)
    zone_share = zone_share.drop(index=0, errors="ignore")  # geomorphic 0 = unmapped
    unmapped_cells = cells[~benthic_mapped].copy()
    unmapped_cells["infill"] = unmapped_cells.n_cells * unmapped_cells.geomorphic.map(zone_share).fillna(0.0)
    infill = unmapped_cells.groupby("site_id").infill.sum()
    infilled_cells = unmapped_cells[unmapped_cells.geomorphic.map(zone_share).notna()].groupby("site_id").n_cells.sum()
    habitable_cells = cells[habitable_cell].groupby("site_id").n_cells.sum()
    n_cells = table.set_index("site_id").n_cells
    table["k_geomorphic_infill"] = ((habitable_cells.reindex(n_cells.index, fill_value=0) +
                                     infill.reindex(n_cells.index, fill_value=0.0)) / n_cells).to_numpy()
    table["frac_geomorphic_infilled"] = (infilled_cells.reindex(n_cells.index, fill_value=0) / n_cells).to_numpy()
    table["frac_unknown_both"] = table.frac_unmapped - table.frac_geomorphic_infilled
    if not table.k_geomorphic_infill.between(0.0, 1.0 + 1e-9).all():
        raise ValueError("Geomorphic-infilled habitable fraction outside [0, 1]")
    table["k_geomorphic_infill"] = table.k_geomorphic_infill.clip(0.0, 1.0)
    class_columns = [f"frac_{name}" for name in CLASSES.values()] + ["frac_other_class", "frac_unmapped"]
    if not np.allclose(table[class_columns].sum(axis=1), 1.0):
        raise ValueError("Class fractions do not sum to one")
    output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output, index=False)
    cells_path = output.with_name(output.stem + "_cells.csv")
    cells.to_csv(cells_path, index=False)
    provenance = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": str(Path(__file__).relative_to(REPO_ROOT)),
        "script_sha256": sha256_file(Path(__file__)),
        "sites": str(sites_path), "sites_sha256": sha256_file(sites_path),
        "benthic_source": str(benthic_path),
        "benthic_source_tree_sha256": sha256_tree(benthic_path),
        "benthic_classes": {str(k): v for k, v in CLASSES.items()},
        "habitable_classes": list(HABITABLE),
        "definition": "k_benthic = (rock + coral/algae cells) / cells inside site polygon; "
                      "sand, rubble, unlabelled and unmapped cells are non-habitable",
        "definition_mapped": f"k_benthic_mapped = (rock + coral/algae) / mapped cells for sites "
                             f">= {MIN_MAPPED} mapped; else reef area-weighted mapped composition; "
                             "else RME 1 - sand_proportion",
        "rme_id_list_sha256": sha256_file(RME_ID_LIST),
        "geomorphic_source": str(geomorphic_path),
        "geomorphic_source_tree_sha256": sha256_tree(geomorphic_path),
        "geomorphic_classes": {str(k): v for k, v in GEOMORPHIC_CLASSES.items()},
        "definition_geomorphic_infill": "benthic-mapped cells keep their class; benthic-unmapped cells "
                                        "contribute p(zone) = habitable / benthic-mapped cells of that "
                                        f"geomorphic zone over domain sites (zones with < {MIN_ZONE_CELLS} "
                                        "mapped cells, or unmapped in both layers, contribute 0)",
        "zone_habitable_share": {str(int(k)): (None if pd.isna(v) else float(v)) for k, v in zone_share.items()},
        "zone_benthic_mapped_cells": {str(int(k)): int(v) for k, v in zone.items()},
        "k_geomorphic_infill_summary": table.k_geomorphic_infill.describe().to_dict(),
        "cells_output": str(cells_path), "cells_output_sha256": sha256_file(cells_path),
        "k_benthic_mapped_summary": table.k_benthic_mapped.describe().to_dict(),
        "k_benthic_mapped_sources": table.k_benthic_mapped_source.value_counts().to_dict(),
        "cell_rule": "cell centre inside polygon (rasterio geometry_mask, all_touched=False); "
                     "EPSG:4326 cells treated as equal area within a site",
        "raster_resolution_deg": list(rasterio.open(benthic_path).res),
        "access": "derived from INTERNAL ONLY data; do not publish or commit outputs",
        "n_sites": int(len(table)),
        "k_benthic_summary": table.k_benthic.describe().to_dict(),
        "output": str(output), "output_sha256": sha256_file(output),
    }
    output.with_suffix(".provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(table[["k_benthic", "k_benthic_mapped", "k_geomorphic_infill",
                 "frac_geomorphic_infilled", "frac_unknown_both"]].describe().round(3).to_string())
    print("zone habitable share:", {GEOMORPHIC_CLASSES.get(int(k), int(k)): round(float(v), 3)
                                    for k, v in zone_share.items() if pd.notna(v)})
    print(table.k_benthic_mapped_source.value_counts().to_string())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sites", type=Path, default=DEFAULT_SITES)
    parser.add_argument("--benthic", type=Path, default=DEFAULT_BENTHIC)
    parser.add_argument("--geomorphic", type=Path, default=DEFAULT_GEOMORPHIC)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    main(args.sites, args.benthic, args.geomorphic, args.output)
