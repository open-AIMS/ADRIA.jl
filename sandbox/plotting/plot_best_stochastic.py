import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# Load stochastic simulation runs (fixed best COTS parameters)
sim_df = pd.read_csv("sandbox/data/best_stochastic_trajectories.csv")
site_path = Path("sandbox/data/best_stochastic_site_trajectories.csv")
site_df = pd.read_csv(site_path) if site_path.exists() else pd.DataFrame()

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
n_runs = sim_df['sim_id'].nunique()

for i, reef in enumerate(reefs_to_plot, 1):
    ax1 = plt.subplot(2, 2, i)
    ax2 = ax1.twinx()
    
    reef_sim = sim_df[sim_df['reef_name'] == reef]
    
    if not reef_sim.empty:
        # 1. Plot individual stochastic simulation runs in light alpha to show envelope
        for sim_id in reef_sim['sim_id'].unique()[:15]:  # Plot up to 15 individual traces
            run_data = reef_sim[reef_sim['sim_id'] == sim_id]
            ax1.plot(run_data['year'], run_data['sim_cots_norm'], color='red', alpha=0.12, linewidth=1.0)
            ax2.plot(run_data['year'], run_data['sim_coral'], color='blue', alpha=0.12, linewidth=1.0)
            
        # 2. Calculate median and 10th-90th percentile uncertainty bounds across runs.
        summary_cots = reef_sim.groupby('year')['sim_cots_norm'].agg(['median', lambda x: np.percentile(x, 10), lambda x: np.percentile(x, 90)]).reset_index()
        summary_cots.columns = ['year', 'median', 'p10', 'p90']
        
        summary_coral = reef_sim.groupby('year')['sim_coral'].agg(['median', lambda x: np.percentile(x, 10), lambda x: np.percentile(x, 90)]).reset_index()
        summary_coral.columns = ['year', 'median', 'p10', 'p90']
        
        # Plot COTS Median + Shading
        ax1.plot(summary_cots['year'], summary_cots['median'], color='red', linewidth=3.0, label='Sim COTS (Median, Best Fit)')
        ax1.fill_between(summary_cots['year'], summary_cots['p10'], summary_cots['p90'], color='red', alpha=0.2, label='COTS 10th-90th %ile Range')
        
        # Plot Coral Median + Shading
        ax2.plot(summary_coral['year'], summary_coral['median'], color='blue', linewidth=3.0, label='Sim Coral Cover (Median)')
        ax2.fill_between(summary_coral['year'], summary_coral['p10'], summary_coral['p90'], color='blue', alpha=0.15, label='Coral 10th-90th %ile Range')

        # 3. Plot within-reef site spread when site-level diagnostics are available.
        if not site_df.empty:
            reef_sites = site_df[site_df['reef_name'] == reef].copy()
            if not reef_sites.empty:
                max_site_cots = reef_sites['sim_cots_adult'].max()
                if max_site_cots > 0:
                    reef_sites['sim_cots_norm'] = reef_sites['sim_cots_adult'] / max_site_cots
                    site_cots = reef_sites.groupby('year')['sim_cots_norm'].agg([
                        lambda x: np.percentile(x, 10),
                        lambda x: np.percentile(x, 90)
                    ]).reset_index()
                    site_cots.columns = ['year', 'p10', 'p90']
                    ax1.fill_between(
                        site_cots['year'],
                        site_cots['p10'],
                        site_cots['p90'],
                        color='orange',
                        alpha=0.16,
                        label='COTS Site 10th-90th %ile Range'
                    )

                site_coral = reef_sites.groupby('year')['sim_coral'].agg([
                    lambda x: np.percentile(x, 10),
                    lambda x: np.percentile(x, 90)
                ]).reset_index()
                site_coral.columns = ['year', 'p10', 'p90']
                ax2.fill_between(
                    site_coral['year'],
                    site_coral['p10'],
                    site_coral['p90'],
                    color='cyan',
                    alpha=0.12,
                    label='Coral Site 10th-90th %ile Range'
                )
    
    # 4. Plot empirical observations
    emp_name = sim_to_emp_map[reef]
    
    # Emp COTS
    emp_reef = emp_df[emp_df['reef_name'] == emp_name]
    if not emp_reef.empty:
        max_cots = emp_reef['cotsptow'].max()
        if max_cots > 0:
            emp_reef_norm = emp_reef['cotsptow'] / max_cots
            ax1.scatter(emp_reef['year'], emp_reef_norm, label='Obs COTS (cots/tow norm)', color='darkred', marker='o', s=65, zorder=5, edgecolor='black', linewidth=0.5)
            
    # Emp Coral
    emp_coral_reef = emp_coral_df[emp_coral_df['reef_name'] == emp_name]
    if not emp_coral_reef.empty:
        ax2.scatter(emp_coral_reef['report_year'], emp_coral_reef['mean'], label='Obs Coral Cover', color='darkblue', marker='s', s=65, zorder=5, edgecolor='black', linewidth=0.5)
        
    # Formatting
    ax1.set_title(f"{reef}\n(Best Parameter Candidate: run_id=2, Stochastic Runs={n_runs})", fontsize=13, fontweight='bold')
    ax1.set_xlabel('Year', fontsize=12)
    ax1.set_ylabel('Normalized Adult COTS Density', color='red', fontsize=12)
    ax2.set_ylabel('Total Coral Cover (Proportion)', color='blue', fontsize=12)
    
    ax1.set_xlim(1985, 2024)
    ax1.set_ylim(-0.05, 1.05)
    ax2.set_ylim(-0.05, 1.05)
    
    ax1.grid(True, linestyle='--', alpha=0.6)
    
    # Combine legends
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    by_label = dict(zip(labels1 + labels2, lines1 + lines2))
    ax1.legend(by_label.values(), by_label.keys(), loc='upper left', fontsize=9, framealpha=0.85)

plt.tight_layout()
plt.savefig("sandbox/best_stochastic_plot.png", dpi=300)
print("Stochastic best-fit plot saved to sandbox/best_stochastic_plot.png")
