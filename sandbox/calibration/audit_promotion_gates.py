"""Promotion gates 4a (held-out years) and 4b (leave-one-reef-out) for the current configuration
(crash_mechanism_protocol.md, "Preregistered promotion checks", 2026-10-19).

Decision families re-ranked (seed means over 5 seed pairs; lower is better):
- start   : lin / h10 / h30 (no cull), runs 20261015T_36 (original layer) and _37 (corrected layer);
            chosen h10, judged on each layer;
- weight_q: correction weight w {1, 2, 3} x Q {10, 15, 20, 25}, runs 20261016T_01-03; chosen (2, 20);
- ext_fec : Ext {5, 2.5} x local output {1, 2, 3}, run 20261015T_45; chosen (5, 1).

4a: fit window = v2 manta 1992-2012; check window = v2 manta 2013-2021 + diver term (rescore_holdout.jl).
4b: summed v3 (primary rows, A_half 0) with each calibration reef dropped in turn.

Pass (fixed in the preregistration; "within the seed spread" made operational here): the chosen
option ranks first, or its paired seed-by-seed difference from the first-ranked option is positive
in fewer than 4 of 5 seeds. 4a additionally requires that it does not rank last on the check window.

    python audit_promotion_gates.py <out_dir> <holdout_dir> [<holdout_dir> ...]
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from audit_1991_switching_repair import sha256  # noqa: E402

RUNS = Path(__file__).resolve().parent / "runs"
REEFS = ("Lizard Island Reef", "MacGillivray Reef", "North Direction Reef", "Eyrie Reef")
RESCORES = {"20261015T_36_init_cull_original_1991_v1": "20261015T_38_rescore_init_cull_original",
            "20261015T_37_init_cull_corrected_1991_v1": "20261015T_40_rescore_init_cull_corrected",
            "20261016T_01_pit_q_w1_1991_v1": "20261016T_04_rescore_pit_q_w1",
            "20261016T_02_pit_q_w2_1991_v1": "20261016T_06_rescore_pit_q_w2",
            "20261016T_03_pit_q_w3_1991_v1": "20261016T_08_rescore_pit_q_w3",
            "20261015T_45_local_ext_facevalue_1991_v1": "20261015T_46_rescore_local_ext"}


def option(run: str, arm: str) -> tuple[str, str, str]:
    """(family, group, option label)."""
    if "init_cull" in run:
        return "start", "original" if "original" in run else "corrected", arm.split("_")[4]
    if "pit_q" in run:
        w = re.search(r"pit_q_w(\d)", run).group(1)
        q = re.search(r"_q(\d+)_", arm).group(1)
        return "weight_q", "all", f"w{w}_q{q}"
    m = re.search(r"_e(\d\d)(?:_f(\d{3}))?_", arm)
    return "ext_fec", "all", f"e{m.group(1)}_f{m.group(2) or '100'}"


CHOSEN = {"start": "h10", "weight_q": "w2_q20", "ext_fec": "e50_f100"}


def judge(table: pd.DataFrame, value: str, check_last: str | None = None) -> list[dict]:
    """table: rows (family, group, option, seed, value...). Returns one record per (family, group)."""
    out = []
    for (fam, grp), g in table.groupby(["family", "group"]):
        means = g.groupby("option")[value].mean().sort_values()
        best, chosen = means.index[0], CHOSEN[fam]
        pivot = g.pivot_table(index="seed", columns="option", values=value)
        worse = int((pivot[chosen] - pivot[best] > 0).sum()) if best != chosen else 0
        rec = {"family": fam, "group": grp, "metric": value, "chosen": chosen, "best": best,
               "chosen_rank": int(list(means.index).index(chosen) + 1), "n_options": len(means),
               "chosen_mean": float(means[chosen]), "best_mean": float(means[best]),
               "seeds_chosen_worse_than_best": worse, "within_spread": best == chosen or worse < 4}
        if check_last:
            cm = g.groupby("option")[check_last].mean().sort_values()
            rec.update({"check_rank": int(list(cm.index).index(chosen) + 1), "check_best": cm.index[0],
                        "check_not_last": list(cm.index).index(chosen) + 1 < len(cm)})
        out.append(rec)
    return out


def main(out: Path, holdouts: list[Path]) -> None:
    if out.exists():
        raise FileExistsError(f"Refusing to overwrite {out}")
    # 4a.
    h = pd.concat([pd.read_csv(d / "holdout_scores.csv") for d in holdouts], ignore_index=True)
    h = h.drop_duplicates(["run", "treatment", "reef_name"])
    arm = h.groupby(["run", "treatment"]).agg(loss_fit=("loss_fit", "sum"), loss_check=("loss_check", "sum")).reset_index()
    arm[["family", "group", "option"]] = [option(r, a) for r, a in zip(arm.run, arm.treatment)]
    arm["seed"] = arm.treatment.str.extract(r"_s(\d)_").astype(int)
    gate4a = judge(arm, "loss_fit", check_last="loss_check")
    for rec in gate4a:
        rec["pass"] = rec["within_spread"] and rec["check_not_last"]
    # 4b.
    rows = []
    for run, rescore in RESCORES.items():
        s = pd.read_csv(RUNS / rescore / "detectability_scores.csv")
        s = s[(s.a_half == 0.0) & (s.score_v2_primary.astype(str).str.lower() == "true")]
        for (treat), g in s.groupby("treatment"):
            fam, grp, opt = option(run, treat)
            if fam == "start" and "cull210" in treat:
                continue
            v = g.set_index("reef_name").loss_v3
            row = {"run": run, "treatment": treat, "family": fam, "group": grp, "option": opt,
                   "seed": int(re.search(r"_s(\d)_", treat).group(1)), "v3_all": float(v.sum())}
            for reef in REEFS:
                row[f"v3_without_{reef.split()[0]}"] = float(v.drop(reef).sum())
            rows.append(row)
    loro = pd.DataFrame(rows)
    gate4b = []
    for col in ["v3_all"] + [f"v3_without_{r.split()[0]}" for r in REEFS]:
        for rec in judge(loro, col):
            rec["pass"] = rec["within_spread"]
            gate4b.append(rec)
    out.mkdir(parents=True)
    pd.DataFrame(gate4a).to_csv(out / "gate4a_holdout.csv", index=False)
    pd.DataFrame(gate4b).to_csv(out / "gate4b_loro.csv", index=False)
    arm.to_csv(out / "holdout_by_arm.csv", index=False)
    loro.to_csv(out / "loro_by_arm.csv", index=False)
    report = {"status": "promotion_gates_4a_4b", "chosen": CHOSEN,
              "pass_rule": "chosen ranks first, or paired difference to first is positive in < 4/5 seeds; 4a also not last on check window",
              "gate4a_pass": all(r["pass"] for r in gate4a), "gate4b_pass": all(r["pass"] for r in gate4b),
              "gate4a": gate4a, "gate4b": gate4b,
              "inputs": {**{f"{d.name}/holdout_scores.csv": sha256(d / "holdout_scores.csv") for d in holdouts},
                         **{f"{r}/detectability_scores.csv": sha256(RUNS / r / "detectability_scores.csv") for r in RESCORES.values()}},
              "script_sha256": sha256(Path(__file__))}
    (out / "gates_4a_4b.json").write_text(json.dumps(report, indent=2, default=str) + "\n")
    pd.set_option("display.width", 250)
    print(pd.DataFrame(gate4a).round(2).to_string(index=False))
    print(pd.DataFrame(gate4b).round(2).to_string(index=False))
    print("4a pass:", report["gate4a_pass"], "| 4b pass:", report["gate4b_pass"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out", type=Path)
    parser.add_argument("holdouts", type=Path, nargs="+")
    a = parser.parse_args()
    main(a.out.resolve(), [d.resolve() for d in a.holdouts])
