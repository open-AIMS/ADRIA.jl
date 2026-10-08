"""Quantify what LTMP manta/photo coral can support for Lizard calibration.

This is a survey-method and support diagnostic, not a fitted observation operator.
The model is area-weighted whole-reef cover; both LTMP series are 9 m surveys.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd

from plot_lizard_ltmp_coral import GATE_RUN, REEF_ALIASES, ROOT, load_ltmp_coral


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
    output = ROOT / "sandbox/calibration/runs" / args.run_id
    if output.exists():
        raise FileExistsError(output)
    output.mkdir(parents=True)
    sources = {
        "manta": ROOT / "sandbox/data/reef_manta.csv",
        "photo": ROOT / "sandbox/data/reef_photo_transect.csv",
        "model": ROOT / "sandbox/calibration/runs" / GATE_RUN / "trajectories.csv",
        "site_map": ROOT / "sandbox/data/Lizard_Historical_v0.2/site_to_reef.csv",
        "domain_sites": ROOT / "sandbox/data/Lizard_Historical_v0.2/spatial/lizard_cluster.gpkg",
    }
    domain_sites = gpd.read_file(sources["domain_sites"])
    site_map = pd.read_csv(sources["site_map"])
    if (domain_sites.site_id.nunique() != len(domain_sites) or
            set(domain_sites.site_id) != set(site_map.site_id) or
            not np.allclose(domain_sites.depth_med, 5.0, rtol=0, atol=0)):
        raise ValueError("V2 polygon/depth support differs from the documented placeholder")
    coral = load_ltmp_coral(sources["manta"], sources["photo"])
    coral.to_csv(output / "mapped_ltmp_coral.csv", index=False)
    manta = coral[coral.method == "manta"]
    photo = coral[coral.method == "photo_transect"]
    paired = manta.merge(photo, on=["model_reef_name", "year"], suffixes=("_manta", "_photo"))
    paired["photo_minus_manta_cover"] = paired.median_photo - paired.median_manta
    paired["survey_date_gap_days"] = (
        pd.to_datetime(paired.survey_date_photo) - pd.to_datetime(paired.survey_date_manta)
    ).dt.days
    paired.to_csv(output / "paired_ltmp_methods.csv", index=False)

    model = pd.read_csv(sources["model"])
    selected = model[(model.treatment == "v2_owen_boundary") & (model.seed == 20260930)]
    if selected.groupby(["reef_name", "year"]).size().ne(1).any():
        raise ValueError("Duplicate model reef-year rows")
    matched = coral.merge(selected[["reef_name", "year", "coral_cover"]],
                          left_on=["model_reef_name", "year"],
                          right_on=["reef_name", "year"], validate="many_to_one")
    matched["whole_reef_model_minus_9m_survey"] = matched.coral_cover - matched["median"]
    matched.to_csv(output / "survey_model_diagnostic.csv", index=False)

    summaries = []
    for reef in REEF_ALIASES:
        reef_manta = manta[manta.model_reef_name == reef]
        reef_photo = photo[photo.model_reef_name == reef]
        reef_pairs = paired[paired.model_reef_name == reef]
        diffs = reef_pairs.photo_minus_manta_cover
        corr = (reef_pairs.median_manta.corr(reef_pairs.median_photo)
                if len(reef_pairs) >= 3 else np.nan)
        summaries.append({
            "reef_name": reef,
            "manta_n_years": len(reef_manta),
            "manta_first_year": int(reef_manta.year.min()),
            "manta_last_year": int(reef_manta.year.max()),
            "photo_n_years": len(reef_photo),
            "paired_n_years": len(reef_pairs),
            "paired_first_year": int(reef_pairs.year.min()) if len(reef_pairs) else None,
            "paired_last_year": int(reef_pairs.year.max()) if len(reef_pairs) else None,
            "photo_minus_manta_mean": float(diffs.mean()) if len(diffs) else None,
            "photo_minus_manta_median": float(diffs.median()) if len(diffs) else None,
            "photo_manta_pearson": float(corr) if np.isfinite(corr) else None,
            "paired_abs_date_gap_median_days": float(reef_pairs.survey_date_gap_days.abs().median())
            if len(reef_pairs) else None,
            "mapped_depth_m": 9,
            "model_depth_med_m": 5,
        })
    pd.DataFrame(summaries).to_csv(output / "coral_mapping_support.csv", index=False)

    for reef in REEF_ALIASES:
        frame = matched[matched.model_reef_name == reef]
        if frame[frame.method == "manta"].empty:
            raise ValueError(f"No mapped manta/model overlap for {reef}")
    output_names = ["mapped_ltmp_coral.csv", "paired_ltmp_methods.csv",
                    "survey_model_diagnostic.csv", "coral_mapping_support.csv"]
    metadata = {
        "status": "diagnostic_only",
        "observation_contract": "LTMP reef-category 9 m medians; report-year match; model area-weighted whole-reef cover",
        "manta_role": "common four-reef survey diagnostic, not calibrated cover likelihood",
        "photo_role": "independent method cross-check on three reefs, not pooled with manta",
        "model_depth_note": f"All {len(domain_sites)} V2 sites have placeholder median depth 5 m; no depth-specific cover observation operator",
        "adria_revision": subprocess.check_output(
            ["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip(),
        "code_sha256": {
            str(path): digest(path) for path in (
                Path(__file__), ROOT / "sandbox/calibration/plot_lizard_ltmp_coral.py",
                ROOT / "sandbox/domain_building/build_lizard_domain_v2.py")
        },
        "source_paths": {name: str(path) for name, path in sources.items()},
        "source_sha256": {name: digest(path) for name, path in sources.items()},
        "output_sha256": {name: digest(output / name) for name in output_names},
    }
    (output / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(output)
    print(pd.DataFrame(summaries).to_string(index=False))


if __name__ == "__main__":
    main()
