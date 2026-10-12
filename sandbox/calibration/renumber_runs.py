"""Give run directories without a time-of-day stamp a per-date sequence number,
`YYYYMMDDT_NN_<name>`, ordered by directory creation time, so file explorers list
them in the order they were tried.

Run contents are never modified (they are the historical record). References in
the calibration docs, scripts and the sibling COTSMod README are rewritten, and
every rename is appended to `runs/RUN_RENAMES.csv` (old name, new name, creation
time). Use `--apply` to rename; the default is a dry run.

    python sandbox/calibration/renumber_runs.py [--apply] [--extra-root DIR ...]
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import os
import re
from pathlib import Path

CAL = Path(__file__).resolve().parent
ROOT = CAL.parents[1]
RUNS = CAL / "runs"
LOG = RUNS / "RUN_RENAMES.csv"
UNSTAMPED = re.compile(r"^(\d{8})T(?!\d{6})(?!_\d{2}_)_?(.+)$")
TEXT_SUFFIXES = {".md", ".py", ".jl", ".toml", ".txt"}
SKIP_DIRS = {RUNS, ROOT / "sandbox" / "data", ROOT / ".git"}


def planned_renames() -> list[tuple[str, str, float]]:
    candidates = []
    for path in RUNS.iterdir():
        match = UNSTAMPED.match(path.name) if path.is_dir() else None
        if match:
            candidates.append((path.stat().st_ctime, path.name, match.group(1), match.group(2)))
    existing = {}
    for path in RUNS.iterdir():
        numbered = re.match(r"^(\d{8})T_(\d{2})_", path.name)
        if path.is_dir() and numbered:
            existing[numbered.group(1)] = max(existing.get(numbered.group(1), 0), int(numbered.group(2)))
    plan = []
    for created, name, date, rest in sorted(candidates):
        existing[date] = existing.get(date, 0) + 1
        plan.append((name, f"{date}T_{existing[date]:02d}_{rest}", created))
    return plan


def text_files(roots: list[Path]):
    for root in roots:
        for dirpath, dirnames, filenames in os.walk(root):
            here = Path(dirpath)
            dirnames[:] = [d for d in dirnames if here / d not in SKIP_DIRS and not d.startswith(".")]
            for filename in filenames:
                path = here / filename
                if path.suffix in TEXT_SUFFIXES:
                    yield path


def main(apply: bool, extra_roots: list[Path]) -> None:
    plan = planned_renames()
    if not plan:
        print("Nothing to renumber")
        return
    for old, new, _ in plan:
        print(f"{old} -> {new}")
    # Longest names first so a name that prefixes another is never replaced inside it.
    patterns = [(re.compile(r"(?<![A-Za-z0-9_])" + re.escape(old) + r"(?![A-Za-z0-9_])"), new)
                for old, new, _ in sorted(plan, key=lambda item: -len(item[0]))]
    edits = {}
    for path in text_files([ROOT, *extra_roots]):
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue
        updated = text
        for pattern, new in patterns:
            updated = pattern.sub(new, updated)
        if updated != text:
            edits[path] = updated
    for path in sorted(edits):
        print("references:", path)
    if not apply:
        print("Dry run; pass --apply to rename")
        return
    for old, new, _ in plan:
        if (RUNS / new).exists():
            raise FileExistsError(f"Refusing to overwrite {new}")
    new_log = not LOG.exists()
    with LOG.open("a", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        if new_log:
            writer.writerow(["old_name", "new_name", "created_local", "renamed_utc"])
        for old, new, created in plan:
            (RUNS / old).rename(RUNS / new)
            writer.writerow([old, new, dt.datetime.fromtimestamp(created).isoformat(timespec="seconds"),
                             dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")])
    for path, text in edits.items():
        path.write_text(text, encoding="utf-8")
    print(f"Renamed {len(plan)} runs; updated {len(edits)} files")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--apply", action="store_true")
    parser.add_argument("--extra-root", type=Path, action="append", default=[])
    arguments = parser.parse_args()
    main(arguments.apply, arguments.extra_root)
