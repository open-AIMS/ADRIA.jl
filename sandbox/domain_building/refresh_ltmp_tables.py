"""Refresh the AIMS LTMP tables used by calibration (sandbox/data/reef_cots.csv, reef_manta.csv)
from a dated COTS-Coral-OutcomesV2 snapshot (scripts/01_download_aims_ltmp.R, which uses
R/aims_api.R). Lead's decision, 2026-10-16, nineteenth round.

- reef_cots: survey-level COTS counts and tows (AIMS `cots-by-domain`). `year` is the AIMS
  reporting year (surveys from about September count toward the next year).
- reef_manta: AIMS's MODELLED manta hard-coral estimates per reef x reporting year.
Both are written in the existing column layout (reef_name = domain_name; snapshot-only columns
query_domain and count_group are dropped). The current tables are first copied to
sandbox/data/archive_ltmp/<old sha256[:12]>/ (never overwritten).

Checks (fail closed):
- no duplicate survey keys; reef names match the query domain;
- every existing COTS survey (survey_id, domain_name) keeps identical year, tows, cots and
  cotsptow; removed surveys are classed as renamed (same survey_id under another reef name) or
  dropped by AIMS, reported, and fail the refresh if any lies in the domain's sectors (CL, CA);
- the change report lists added surveys and coral revisions.

    python refresh_ltmp_tables.py --snapshot <COTS-Coral-OutcomesV2/data/raw/aims-ltmp/YYYY-MM-DD> [--photo]

--photo (added 2026-10-19) instead refreshes sandbox/data/reef_photo_transect.csv from reef_photo.parquet
(AIMS modelled photo-transect cover; HARD CORAL rows with purpose COMPOSITION carry taxa in
reefpage_category, GROUP_LEVEL rows the totals), same layout, previous copy archived; report in
ltmp_refresh_<date>_photo/.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "sandbox" / "data"
TABLES = {"reef_cots": DATA / "reef_cots.csv", "reef_manta": DATA / "reef_manta.csv"}
ARCHIVE = DATA / "archive_ltmp"
REPORT_ROOT = DATA / "derived"
COTS_KEY = ["survey_id", "domain_name"]
DOMAIN_SECTORS = ("CL", "CA")  # Cooktown-Lizard and Cairns: the Lizard domain's AIMS sectors
MANTA_KEY = ["domain_name", "report_year", "depth", "series", "variable", "project_code"]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


PHOTO = DATA / "reef_photo_transect.csv"
PHOTO_KEY = ["domain_name", "report_year", "depth", "series", "variable", "purpose", "reefpage_category", "project_code"]


def refresh_photo(snapshot: Path) -> None:
    report_dir = REPORT_ROOT / f"ltmp_refresh_{snapshot.name}_photo"
    if report_dir.exists():
        raise FileExistsError(f"Refusing to overwrite {report_dir}")
    old = pd.read_csv(PHOTO, low_memory=False)
    new = pd.read_parquet(snapshot / "reef_photo.parquet")
    if (new.domain_name != new.query_domain).any():
        raise ValueError("reef_photo: domain_name differs from query_domain")
    photo = new.drop(columns=[c for c in ("query_domain", "count_group") if c in new]).assign(reef_name=lambda d: d.domain_name)
    photo = photo[list(old.columns)]
    for col, dtype in old.dtypes.items():
        if pd.api.types.is_integer_dtype(dtype):
            photo[col] = pd.to_numeric(photo[col]).round().astype("Int64")
    o = old.assign(depth=old.depth.astype(float), reefpage_category=old.reefpage_category.fillna("_"))
    n = photo.assign(depth=photo.depth.astype(float), reefpage_category=photo.reefpage_category.fillna("_"))
    if n.duplicated(PHOTO_KEY).any():
        raise ValueError("Duplicate photo-transect keys")
    b = o.merge(n, on=PHOTO_KEY, how="outer", suffixes=("_old", "_new"), indicator=True)
    both = b[b._merge == "both"]
    revision = (both.median_new - both.median_old).abs()
    report_dir.mkdir(parents=True)
    digest = sha256(PHOTO)
    dest = ARCHIVE / digest[:12]
    dest.mkdir(parents=True, exist_ok=True)
    if not (dest / PHOTO.name).exists():
        shutil.copy2(PHOTO, dest / PHOTO.name)
    if sha256(dest / PHOTO.name) != digest:
        raise ValueError("Archive copy of reef_photo_transect.csv does not match")
    photo.to_csv(PHOTO, index=False)
    meta = {"snapshot": str(snapshot), "snapshot_input": sha256(snapshot / "reef_photo.parquet"),
            "archived_previous": {"sha256": digest, "archive": str((dest / PHOTO.name).relative_to(ROOT))},
            "rows_old": len(old), "rows_new": len(photo), "matched": len(both),
            "added": int((b._merge == "right_only").sum()), "removed": int((b._merge == "left_only").sum()),
            "added_by_year": {str(int(k)): int(v) for k, v in b[b._merge == "right_only"].report_year.value_counts().sort_index().items()},
            "median_revised_rows": int((revision > 1e-6).sum()), "median_revision_max_abs": float(revision.max()),
            "output": sha256(PHOTO), "script_sha256": sha256(Path(__file__))}
    (report_dir / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps({k: meta[k] for k in ("rows_old", "rows_new", "matched", "added", "removed", "added_by_year",
                                            "median_revised_rows", "median_revision_max_abs")}, indent=1))


def main(snapshot: Path) -> None:
    report_dir = REPORT_ROOT / f"ltmp_refresh_{snapshot.name}"
    if report_dir.exists():
        raise FileExistsError(f"Refusing to overwrite {report_dir}")
    old = {k: pd.read_csv(p) for k, p in TABLES.items()}
    new_cots = pd.read_parquet(snapshot / "reef_cots.parquet")
    new_manta = pd.read_parquet(snapshot / "reef_manta.parquet")
    for name, d in (("reef_cots", new_cots), ("reef_manta", new_manta)):
        if (d.domain_name != d.query_domain).any():
            raise ValueError(f"{name}: domain_name differs from query_domain")

    # COTS: same layout as before.
    cots = new_cots.drop(columns="query_domain").assign(reef_name=lambda d: d.domain_name)
    for col in ("year", "tows", "cots"):
        if not np.allclose(cots[col], cots[col].round()):
            raise ValueError(f"reef_cots.{col} is not integer-valued")
        cots[col] = cots[col].round().astype(int)
    cots = cots[list(old["reef_cots"].columns)]
    if cots.duplicated(COTS_KEY).any():
        raise ValueError("Duplicate COTS survey keys")
    both = old["reef_cots"].merge(cots, on=COTS_KEY, how="outer", suffixes=("_old", "_new"), indicator=True)
    matched = both[both._merge == "both"]
    changed = {c: int((np.abs(matched[f"{c}_old"].astype(float) - matched[f"{c}_new"].astype(float)) > 1e-9).sum())
               for c in ("year", "tows", "cots", "cotsptow")}
    if any(changed.values()):
        raise ValueError(f"Existing COTS surveys changed: {changed}")
    removed = both[both._merge == "left_only"].copy()
    renamed_to = cots.set_index("survey_id").domain_name
    removed["status"] = np.where(removed.survey_id.isin(renamed_to.index), "renamed", "dropped_by_aims")
    removed["renamed_to"] = removed.survey_id.map(renamed_to)
    if removed.a_sector_old.isin(DOMAIN_SECTORS).any():
        raise ValueError("COTS surveys removed in the domain's sectors")
    added = both[both._merge == "right_only"]

    # Manta coral: same layout (reef_name added, count_group dropped).
    manta = new_manta.drop(columns=["query_domain", "count_group"]).assign(reef_name=lambda d: d.domain_name)
    manta = manta[list(old["reef_manta"].columns)]
    for col, dtype in old["reef_manta"].dtypes.items():
        if pd.api.types.is_integer_dtype(dtype):
            manta[col] = pd.to_numeric(manta[col]).round().astype("Int64")
    mo = old["reef_manta"].assign(depth=lambda d: d.depth.astype(float))
    mb = mo.merge(manta.assign(depth=lambda d: d.depth.astype(float)), on=MANTA_KEY, how="outer",
                  suffixes=("_old", "_new"), indicator=True)
    mm = mb[mb._merge == "both"]
    revision = (mm.median_new - mm.median_old).abs()

    report_dir.mkdir(parents=True)
    archived = {}
    for name, path in TABLES.items():
        digest = sha256(path)
        dest = ARCHIVE / digest[:12]
        dest.mkdir(parents=True, exist_ok=True)
        target = dest / path.name
        if not target.exists():
            shutil.copy2(path, target)
        if sha256(target) != digest:
            raise ValueError(f"Archive copy of {path.name} does not match")
        archived[name] = {"sha256": digest, "archive": str(target.relative_to(ROOT))}
    cots.to_csv(TABLES["reef_cots"], index=False)
    manta.to_csv(TABLES["reef_manta"], index=False)
    added[["domain_name", "a_sector_new", "survey_id", "year_new", "survey_date_new", "tows_new", "cots_new", "cotsptow_new"]] \
        .rename(columns=lambda c: c.replace("_new", "")).to_csv(report_dir / "added_cots_surveys.csv", index=False)
    removed[["domain_name", "a_sector_old", "survey_id", "status", "renamed_to", "year_old", "tows_old", "cots_old"]].to_csv(report_dir / "removed_cots_surveys.csv", index=False)
    meta = {"snapshot": str(snapshot), "snapshot_manifest_sha256": sha256(snapshot / "manifest.csv"),
            "snapshot_inputs": {f: sha256(snapshot / f) for f in ("reef_cots.parquet", "reef_manta.parquet")},
            "archived_previous": archived,
            "cots": {"rows_old": len(old["reef_cots"]), "rows_new": len(cots), "added": len(added), "removed": len(removed),
                     "removed_reefs": sorted(set(removed.domain_name)), "removed_status": removed.status.value_counts().to_dict(),
                     "matched_unchanged": len(matched)},
            "manta": {"rows_old": len(old["reef_manta"]), "rows_new": len(manta), "matched": len(mm),
                      "added": int((mb._merge == "right_only").sum()), "removed": int((mb._merge == "left_only").sum()),
                      "median_revised_rows": int((revision > 1e-6).sum()), "median_revision_max_abs": float(revision.max())},
            "outputs": {p.name: sha256(p) for p in TABLES.values()}, "script_sha256": sha256(Path(__file__))}
    (report_dir / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps({k: meta[k] for k in ("cots", "manta")}, indent=1))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snapshot", type=Path, required=True)
    parser.add_argument("--photo", action="store_true")
    a = parser.parse_args()
    (refresh_photo if a.photo else main)(a.snapshot.resolve())
