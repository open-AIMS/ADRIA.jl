from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


SIM_TO_EMP_MAP = {
    "Lizard Island Reef": "Lizard Isles",
    "MacGillivray Reef": "Macgillivray Reef",
    "North Direction Reef": "North Direction Island",
    "Eyrie Reef": "Eyrie Reef",
}

REEFS_TO_PLOT = [
    "Lizard Island Reef",
    "MacGillivray Reef",
    "North Direction Reef",
    "Eyrie Reef",
]

VARIANTS = [
    ("nopulse", "No pulse", "black"),
    ("pulse_start20_dur1", "Pulse start 20, dur 1", "red"),
    ("pulse_start20_dur2", "Pulse start 20, dur 2", "darkorange"),
    ("pulse_start20_dur3", "Pulse start 20, dur 3", "purple"),
]


def validation_metrics(sim_series, obs_series):
    paired = pd.concat([sim_series.rename("sim"), obs_series.rename("obs")], axis=1).dropna()
    n = len(paired)
    if n == 0:
        return {"n": 0, "pearson": np.nan, "spearman": np.nan, "rmse": np.nan, "pbias": np.nan}

    sim = paired["sim"].astype(float)
    obs = paired["obs"].astype(float)
    pearson = np.corrcoef(sim, obs)[0, 1] if n >= 2 and sim.std() > 0 and obs.std() > 0 else np.nan
    sim_rank = sim.rank(method="average")
    obs_rank = obs.rank(method="average")
    spearman = np.corrcoef(sim_rank, obs_rank)[0, 1] if n >= 2 and sim_rank.std() > 0 and obs_rank.std() > 0 else np.nan
    rmse = np.sqrt(np.mean((sim - obs) ** 2))
    pbias = 100.0 * (sim.sum() - obs.sum()) / obs.sum() if obs.sum() != 0 else np.nan
    return {"n": n, "pearson": pearson, "spearman": spearman, "rmse": rmse, "pbias": pbias}


def load_variant(tag):
    path = Path(f"sandbox/data/best_stochastic_{tag}_trajectories.csv")
    if not path.exists():
        return None
    return pd.read_csv(path)


emp_df = pd.read_csv("sandbox/data/reef_cots.csv")
variant_data = [(tag, label, color, load_variant(tag)) for tag, label, color in VARIANTS]
variant_data = [(tag, label, color, df) for tag, label, color, df in variant_data if df is not None]

if not variant_data:
    raise FileNotFoundError("No requested pulse scenario trajectory CSVs were found under sandbox/data.")

metrics_rows = []
fig, axes = plt.subplots(2, 2, figsize=(16, 11), sharex=True, sharey=True)

for ax, reef in zip(axes.flat, REEFS_TO_PLOT):
    emp_name = SIM_TO_EMP_MAP[reef]
    emp_reef = emp_df[emp_df["reef_name"] == emp_name].copy()
    if not emp_reef.empty and emp_reef["cotsptow"].max() > 0:
        emp_reef["obs_norm"] = emp_reef["cotsptow"] / emp_reef["cotsptow"].max()
        ax.scatter(
            emp_reef["year"],
            emp_reef["obs_norm"],
            color="darkred",
            s=42,
            edgecolor="black",
            linewidth=0.4,
            zorder=4,
            label="Observed COTS",
        )

    for tag, label, color, df in variant_data:
        reef_df = df[df["reef_name"] == reef]
        if reef_df.empty:
            continue

        median = reef_df.groupby("year")["sim_cots_norm"].median()
        ax.plot(median.index, median.values, color=color, linewidth=2.0, label=label)

        if not emp_reef.empty and emp_reef["cotsptow"].max() > 0:
            obs = emp_reef.groupby("year")["obs_norm"].mean()
            merged = pd.merge(
                median.rename("sim").reset_index(),
                obs.rename("obs").reset_index(),
                on="year",
                how="inner",
            )
            metrics = validation_metrics(merged["sim"], merged["obs"])
            metrics_rows.append(
                {
                    "variant": tag,
                    "reef_name": reef,
                    "n_obs": metrics["n"],
                    "pearson": metrics["pearson"],
                    "spearman": metrics["spearman"],
                    "rmse": metrics["rmse"],
                    "percent_bias": metrics["pbias"],
                    "overall_peak_year": int(median.idxmax()),
                    "post_2010_peak_year": int(median[median.index >= 2010].idxmax()),
                }
            )

    ax.set_title(reef, fontsize=13, fontweight="bold")
    ax.set_xlim(1985, 2024)
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, linestyle="--", alpha=0.45)

handles, labels = axes.flat[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center", ncol=5, frameon=False)
fig.supylabel("Normalized Adult COTS Density")
fig.supxlabel("Year")
fig.suptitle("Pulse Timing Scenario Comparison", fontsize=16, fontweight="bold", y=0.98)
plt.tight_layout(rect=(0, 0, 1, 0.93))

metrics_df = pd.DataFrame(metrics_rows)
metrics_df.to_csv("sandbox/data/pulse_start20_scenario_comparison_metrics.csv", index=False)
plt.savefig("sandbox/pulse_start20_scenario_comparison.png", dpi=300)
print("Saved comparison plot to sandbox/pulse_start20_scenario_comparison.png")
print("Saved comparison metrics to sandbox/data/pulse_start20_scenario_comparison_metrics.csv")

