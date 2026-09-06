import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# Setup Paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.normpath(os.path.join(SCRIPT_DIR, "..", ".."))
DATA_DIR = os.path.join(REPO_ROOT, "sandbox", "data")
PLOTS_DIR = os.path.join(SCRIPT_DIR, "plots")
os.makedirs(PLOTS_DIR, exist_ok=True)

# Set style
plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica']

# Load Data
sim_path = os.path.join(DATA_DIR, "best_calibrated_trajectories.csv")
if not os.path.exists(sim_path):
    sim_path = os.path.join(DATA_DIR, "best_bbo_simulated_trajectories.csv")

obs_path = os.path.join(DATA_DIR, "reef_cots.csv")
cand_path = os.path.join(DATA_DIR, "bbo_evaluated_candidates.csv")
best_path = os.path.join(DATA_DIR, "bbo_best_summary.csv")

if not os.path.exists(sim_path):
    raise FileNotFoundError(f"Simulated trajectories missing: {sim_path}. Run simulate_best_calibration.jl first.")


sim_df = pd.read_csv(sim_path)
obs_df = pd.read_csv(obs_path)
cand_df = pd.read_csv(cand_path) if os.path.exists(cand_path) else None
best_df = pd.read_csv(best_path) if os.path.exists(best_path) else None

target_reefs = [
    ("Lizard Island Reef", "Lizard Isles"),
    ("MacGillivray Reef", "Macgillivray Reef"),
    ("North Direction Reef", "North Direction Island"),
    ("Eyrie Reef", "Eyrie Reef")
]

# -------------------------------------------------------------
# Figure 1: Normalized COTS Calibration Trajectories vs Observations
# -------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(14, 10), sharex=True, sharey=True)
axes = axes.flatten()

color_sim = "#1f77b4" # Deep blue
color_obs = "#d62728" # Red

for i, (sim_name, obs_name) in enumerate(target_reefs):
    ax = axes[i]
    
    # Filter sim
    sub_sim = sim_df[sim_df['reef_name'] == sim_name].sort_values('year')
    # Filter obs
    sub_obs = obs_df[obs_df['reef_name'] == obs_name].sort_values('year')
    
    if not sub_obs.empty and sub_obs['cotsptow'].max() > 0:
        obs_max = sub_obs['cotsptow'].max()
        # Group by year for obs mean
        obs_yearly = sub_obs.groupby('year')['cotsptow'].mean().reset_index()
        obs_yearly['obs_norm'] = obs_yearly['cotsptow'] / obs_max
        ax.scatter(obs_yearly['year'], obs_yearly['obs_norm'], color=color_obs, s=45, zorder=4, label='Historical Obs (Normalized PTOW)')
        ax.plot(obs_yearly['year'], obs_yearly['obs_norm'], color=color_obs, linestyle='--', alpha=0.6, zorder=3)

    if not sub_sim.empty:
        ax.plot(sub_sim['year'], sub_sim['sim_cots_norm'], color=color_sim, linewidth=2.5, zorder=5, label='Calibrated ADRIA Simulation')
        ax.fill_between(sub_sim['year'], 0, sub_sim['sim_cots_norm'], color=color_sim, alpha=0.15, zorder=2)

    ax.set_title(f"{sim_name}", fontsize=13, fontweight='bold', pad=8)
    ax.set_ylim(-0.05, 1.1)
    ax.set_xlim(1984, 2025)
    ax.grid(True, linestyle=':', alpha=0.6)
    
    if i >= 2:
        ax.set_xlabel("Year", fontsize=11, fontweight='bold')
    if i % 2 == 0:
        ax.set_ylabel("Normalized Adult COTS Density", fontsize=11, fontweight='bold')
    
    if i == 0:
        ax.legend(loc='upper right', frameon=True, framealpha=0.9, fontsize=10)

fig.suptitle("BlackBoxOptim Calibrated COTS Outbreak Trajectories vs Lizard Island Observations", fontsize=15, fontweight='bold', y=0.98)
plt.tight_layout(rect=[0, 0, 1, 0.96])
fig1_path = os.path.join(PLOTS_DIR, "cots_calibration_trajectories.png")
plt.savefig(fig1_path, dpi=300)
plt.close()
print(f"Saved: {fig1_path}")


# -------------------------------------------------------------
# Figure 2: Dual-Axis COTS & Coral Cover Trajectories
# -------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(14, 10), sharex=True)
axes = axes.flatten()

color_cots = "#d62728"
color_coral = "#2ca02c" # Green

for i, (sim_name, _) in enumerate(target_reefs):
    ax1 = axes[i]
    sub_sim = sim_df[sim_df['reef_name'] == sim_name].sort_values('year')
    
    if not sub_sim.empty:
        # COTS Adult Raw Density
        line1 = ax1.plot(sub_sim['year'], sub_sim['sim_cots_adult'], color=color_cots, linewidth=2.2, label='Simulated Adult COTS (ind/ha)')
        ax1.set_ylabel("Adult COTS Density (ind/ha)", color=color_cots, fontsize=11, fontweight='bold')
        ax1.tick_params(axis='y', labelcolor=color_cots)
        
        # Total Coral Cover
        ax2 = ax1.twinx()
        line2 = ax2.plot(sub_sim['year'], sub_sim['sim_coral_cover'] * 100, color=color_coral, linewidth=2.2, linestyle='--', label='Total Coral Cover (%)')
        ax2.set_ylabel("Coral Cover (%)", color=color_coral, fontsize=11, fontweight='bold')
        ax2.tick_params(axis='y', labelcolor=color_coral)
        ax2.grid(False) # Turn off grid for second axis to avoid overlap

    ax1.set_title(f"{sim_name} - COTS & Coral Cover Dynamics", fontsize=12, fontweight='bold', pad=8)
    ax1.set_xlim(1984, 2025)
    ax1.grid(True, linestyle=':', alpha=0.5)

    if i >= 2:
        ax1.set_xlabel("Year", fontsize=11, fontweight='bold')

fig.suptitle("Simulated COTS Population Outbreaks & Coral Cover Responses (Lizard Island Cluster)", fontsize=15, fontweight='bold', y=0.98)
plt.tight_layout(rect=[0, 0, 1, 0.96])
fig2_path = os.path.join(PLOTS_DIR, "cots_coral_cover_trajectories.png")
plt.savefig(fig2_path, dpi=300)
plt.close()
print(f"Saved: {fig2_path}")


# -------------------------------------------------------------
# Figure 3: BlackBoxOptim Candidate Diagnostics & Penalty Breakdown
# -------------------------------------------------------------
if cand_df is not None and not cand_df.empty:
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
    
    # Panel 1: Loss progression
    ax = axes[0]
    ax.plot(cand_df['eval_id'], cand_df['loss'], marker='o', color='#1f77b4', linewidth=1.5, markersize=5, label='Total Combined Loss')
    ax.plot(cand_df['eval_id'], cand_df['cycle_loss'], marker='s', color='#ff7f0e', linestyle='--', linewidth=1.5, markersize=4, label='Cycle Loss')
    ax.set_title("A. Loss Score Progression Across Evaluated Candidates", fontsize=11, fontweight='bold')
    ax.set_xlabel("Evaluation Step ID", fontsize=10, fontweight='bold')
    ax.set_ylabel("Loss Score (Lower is Better)", fontsize=10, fontweight='bold')
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(True, linestyle=':', alpha=0.6)
    
    # Panel 2: Cycle Penalty Sub-component Breakdown for Top Candidates
    ax = axes[1]
    top_cands = cand_df.sort_values('loss').head(10)
    penalty_cols = [
        ('mean_peak_count_penalty', 'Peak Count'),
        ('mean_peak_timing_penalty', 'Peak Timing'),
        ('mean_period_penalty', 'Periodicity'),
        ('mean_lag_correlation_penalty', 'Lag Correlation'),
        ('mean_amplitude_penalty', 'Amplitude')
    ]
    
    bottom = np.zeros(len(top_cands))
    x_labels = [f"#{row.eval_id}" for _, row in top_cands.iterrows()]
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    
    for idx, (p_col, label) in enumerate(penalty_cols):
        if p_col in top_cands.columns:
            vals = top_cands[p_col].values
            ax.bar(x_labels, vals, bottom=bottom, label=label, color=colors[idx], alpha=0.85)
            bottom += vals

    ax.set_title("B. Cycle Metric Penalty Breakdown (Top 10 Candidates)", fontsize=11, fontweight='bold')
    ax.set_xlabel("Candidate Eval ID", fontsize=10, fontweight='bold')
    ax.set_ylabel("Penalty Component Score", fontsize=10, fontweight='bold')
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(True, linestyle=':', alpha=0.5)

    # Panel 3: Parameter Space Exploration (a_F vs a_S vs Loss)
    ax = axes[2]
    scatter = ax.scatter(cand_df['a_F'], cand_df['a_S'], c=cand_df['loss'], cmap='viridis_r', s=60, edgecolors='k', linewidth=0.5)
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label("Loss Score", fontsize=10, fontweight='bold')
    
    if best_df is not None and not best_df.empty:
        best_row = best_df.iloc[0]
        ax.scatter(best_row['a_F'], best_row['a_S'], color='red', marker='*', s=200, label=f"Best Candidate (Loss={best_row['loss']:.2f})", zorder=5)

    ax.set_title("C. Parameter Search Space (a_F vs a_S)", fontsize=11, fontweight='bold')
    ax.set_xlabel("Feeding Rate Attenuation Factor (a_F)", fontsize=10, fontweight='bold')
    ax.set_ylabel("Survival Attenuation Factor (a_S)", fontsize=10, fontweight='bold')
    ax.legend(loc='upper left', fontsize=9)
    ax.grid(True, linestyle=':', alpha=0.6)

    fig.suptitle("BlackBoxOptim Driver Diagnostic & Convergence Analysis", fontsize=14, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    fig3_path = os.path.join(PLOTS_DIR, "bbo_optimization_diagnostics.png")
    plt.savefig(fig3_path, dpi=300)
    plt.close()
    print(f"Saved: {fig3_path}")

print("=== Calibration plotting completed successfully ===")
