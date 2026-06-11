import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load simulated top 5 ensemble
sim_df = pd.read_csv("sandbox/data/top5_trajectories.csv")

# Load empirical data
emp_df = pd.read_csv("sandbox/data/reef_cots.csv")
emp_coral_df = pd.read_csv("sandbox/data/reef_manta.csv")

sim_to_emp_map = {
    "Lizard Island Reef": "Lizard Isles",
    "MacGillivray Reef": "Macgillivray Reef",
    "North Direction Reef": "North Direction Island",
    "Eyrie Reef": "Eyrie Reef"
}

reefs_to_plot = ["Lizard Island Reef", "MacGillivray Reef", "North Direction Reef", "Eyrie Reef"]

plt.figure(figsize=(16, 12))

for i, reef in enumerate(reefs_to_plot, 1):
    ax1 = plt.subplot(2, 2, i)
    ax2 = ax1.twinx()
    
    # 1. Plot simulated ensemble trajectories
    reef_sim = sim_df[sim_df['reef_name'] == reef]
    
    for run_id in reef_sim['run_id'].unique():
        run_data = reef_sim[reef_sim['run_id'] == run_id]
        ax1.plot(run_data['year'], run_data['sim_cots_norm'], color='red', alpha=0.3, linewidth=1.5)
        ax2.plot(run_data['year'], run_data['sim_coral'], color='blue', alpha=0.3, linewidth=1.5)
        
    # Calculate and plot the ensemble mean
    if not reef_sim.empty:
        mean_sim = reef_sim.groupby('year').mean(numeric_only=True).reset_index()
        ax1.plot(mean_sim['year'], mean_sim['sim_cots_norm'], color='red', linewidth=3, label='Sim COTS (Mean)')
        ax2.plot(mean_sim['year'], mean_sim['sim_coral'], color='blue', linewidth=3, label='Sim Coral (Mean)')
    
    # 2. Plot empirical observations
    emp_name = sim_to_emp_map[reef]
    
    # Emp COTS
    emp_reef = emp_df[emp_df['reef_name'] == emp_name]
    if not emp_reef.empty:
        max_cots = emp_reef['cotsptow'].max()
        if max_cots > 0:
            emp_reef_norm = emp_reef['cotsptow'] / max_cots
            ax1.scatter(emp_reef['year'], emp_reef_norm, label='Obs COTS', color='darkred', marker='o', s=60, zorder=5)
            
    # Emp Coral
    emp_coral_reef = emp_coral_df[emp_coral_df['reef_name'] == emp_name]
    if not emp_coral_reef.empty:
        ax2.scatter(emp_coral_reef['report_year'], emp_coral_reef['mean'], label='Obs Coral', color='darkblue', marker='s', s=60, zorder=5)
        
    # Formatting
    ax1.set_title(reef, fontsize=14, fontweight='bold')
    ax1.set_xlabel('Year', fontsize=12)
    ax1.set_ylabel('Normalized Adult COTS Density', color='red', fontsize=12)
    ax2.set_ylabel('Total Coral Cover (Proportion)', color='blue', fontsize=12)
    
    ax1.set_xlim(1985, 2024)
    ax1.set_ylim(-0.05, 1.05)
    ax2.set_ylim(-0.05, 1.05)
    
    ax1.grid(True, linestyle='--', alpha=0.6)
    
    # Combine legends from both axes
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    # Avoid duplicate labels
    by_label = dict(zip(labels1 + labels2, lines1 + lines2))
    ax1.legend(by_label.values(), by_label.keys(), loc='upper left', fontsize=10)

plt.tight_layout()
plt.savefig("sandbox/ensemble_calibration_plot.png", dpi=300)
print("Ensemble plot saved to sandbox/ensemble_calibration_plot.png")
