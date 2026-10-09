"""Index every calibration run directory without interpreting its ecology.

The curated decisions live in MODEL_LOG.md. This deterministic inventory keeps
unreviewed, incomplete, failed and superseded run directories discoverable.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
import tomllib
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNS = HERE / "runs"
LOG = HERE / "MODEL_LOG.md"
OUT = HERE / "MODEL_RUN_INVENTORY.csv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def failure_count(path: Path) -> str:
    if not path.is_file():
        return ""
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return str(sum(1 for _ in csv.DictReader(stream)))


def main() -> None:
    reviewed = set(re.findall(r"runs/([A-Za-z0-9_.-]+)", LOG.read_text()))
    rows = []
    for run in sorted(path for path in RUNS.iterdir()
                      if path.is_dir() and not path.name.startswith(".")):
        metadata = next((run / name for name in ("metadata.toml", "metadata.json")
                         if (run / name).is_file()), None)
        status = "missing_metadata"
        parse_error = ""
        if metadata:
            try:
                if metadata.suffix == ".toml":
                    with metadata.open("rb") as stream:
                        record = tomllib.load(stream)
                else:
                    record = json.loads(metadata.read_text(encoding="utf-8"))
                status = str(record.get("status", "unspecified"))
            except (OSError, ValueError, tomllib.TOMLDecodeError) as error:
                status = "metadata_parse_failed"
                parse_error = type(error).__name__
        files = {path.name for path in run.iterdir() if path.is_file()}
        rows.append({
            "run_id": run.name,
            "decision_log_linked": str(run.name in reviewed).lower(),
            "metadata_file": metadata.name if metadata else "",
            "metadata_sha256": sha256(metadata) if metadata else "",
            "metadata_status": status,
            "metadata_parse_error": parse_error,
            "scores_csv": str("scores.csv" in files).lower(),
            "trajectories_csv": str("trajectories.csv" in files).lower(),
            "failures_csv_rows": failure_count(run / "failures.csv"),
            "file_count": len(files),
        })
    with OUT.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Indexed {len(rows)} run directories; {sum(row['decision_log_linked'] == 'true' for row in rows)} linked in model log: {OUT}")


if __name__ == "__main__":
    main()
