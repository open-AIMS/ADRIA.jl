import pandas as pd
import matplotlib.pyplot as plt
import os

sim_df = pd.read_csv("sandbox/data/calibration_results_sim.csv")
# Read raw empirical data instead of pre-filtered to avoid missing mapped names
emp_df = pd.read_csv("sandbox/data/reef_cots.csv")
emp_coral_df = pd.read_csv("sandbox/data/reef_manta.csv")

# Map simulation reef names to empirical AIMS reef names
sim_to_emp_map = {
    "Lizard Island Reef": "Lizard Isles",
    "MacGillivray Reef": "Macgillivray Reef",
    "North Direction Reef": "North Direction Island",
    "Eyrie Reef": "Eyrie Reef"
}

reefs_to_plot = ["Lizard Island Reef", "MacGillivray Reef", "North Direction Reef", "Eyrie Reef"]

plt.figure(figsize=(16, 12))

for i, reef in enumerate(reefs_to_plot, 1):
    plt.subplot(2, 2, i)
    
    # Create twin axes
    ax1 = plt.gca()
    ax2 = ax1.twinx()
    
    # Sim COTS
    sim_reef = sim_df[sim_df['reef_name'] == reef]
    ax1.plot(sim_reef['year'], sim_reef['sim_cots_norm'], label='Simulated COTS (Norm)', color='red', linewidth=2)
    
    # Emp COTS
    emp_name = sim_to_emp_map[reef]
    emp_reef = emp_df[emp_df['reef_name'] == emp_name]
    
    if not emp_reef.empty:
        # Normalize empirical COTS for this specific reef
        max_cots = emp_reef['cotsptow'].max()
        if max_cots > 0:
            emp_reef_norm = emp_reef['cotsptow'] / max_cots
            emp_reef_sorted = emp_reef.sort_values(by='year')
            emp_reef_norm_sorted = emp_reef_norm.loc[emp_reef_sorted.index]
            ax1.plot(emp_reef_sorted['year'], emp_reef_norm_sorted, label='Observed COTS (Norm)', color='lightcoral', marker='o', linestyle='--', linewidth=1.5, markerfacecolor='darkred', markeredgecolor='darkred', zorder=5)
    
    # Sim Coral
    ax2.plot(sim_reef['year'], sim_reef['sim_coral'], label='Simulated Coral Cover', color='blue', linewidth=2)
    
    # Emp Coral
    emp_coral_reef = emp_coral_df[emp_coral_df['reef_name'] == emp_name]
    if not emp_coral_reef.empty:
        emp_coral_sorted = emp_coral_reef.sort_values(by='report_year')
        ax2.plot(emp_coral_sorted['report_year'], emp_coral_sorted['mean'], label='Observed Coral Cover', color='cornflowerblue', marker='s', linestyle='--', linewidth=1.5, markerfacecolor='darkblue', markeredgecolor='darkblue', zorder=5)
    
    plt.title(reef, fontsize=14, fontweight='bold')
    ax1.set_xlabel('Year', fontsize=12)
    ax1.set_ylabel('Normalized COTS Density', color='darkred', fontsize=12)
    ax2.set_ylabel('Coral Cover (Proportion)', color='darkblue', fontsize=12)
    
    ax1.set_xlim(1985, 2024)
    ax1.set_ylim(0, 1.1)
    ax2.set_ylim(0, 1.0)
    
    # Combine legends
    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()
    ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc='upper right')
    ax1.grid(True, linestyle='--', alpha=0.5)

plt.tight_layout()
plt.savefig("sandbox/calibration_plot.png", dpi=300)
print("Plot saved to sandbox/calibration_plot.png")
