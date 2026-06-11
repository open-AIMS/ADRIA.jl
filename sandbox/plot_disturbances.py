import pandas as pd
import matplotlib.pyplot as plt

# Load data
dhw_df = pd.read_csv("sandbox/data/dhw_historical.csv")
site_df = pd.read_csv("sandbox/data/Lizard_Historical_v0.1/site_to_reef.csv")

# Clean reef names in site mapping
site_df['reef_name_clean'] = site_df['reef_name'].apply(lambda x: x.split(" (")[0])

# Merge DHW with site mapping to get reef names
merged = pd.merge(dhw_df, site_df, left_on="RME_UNIQUE_ID", right_on="UNIQUE_ID")

reefs_to_plot = ["Lizard Island Reef", "MacGillivray Reef", "North Direction Reef", "Eyrie Reef"]

plt.figure(figsize=(16, 12))

for i, reef in enumerate(reefs_to_plot, 1):
    plt.subplot(2, 2, i)
    
    # Filter for this reef
    reef_data = merged[merged['reef_name_clean'] == reef]
    
    if not reef_data.empty:
        # Group by year (timestep) and calculate mean DHW across all sites in the reef
        yearly_dhw = reef_data.groupby('timestep')['dhw'].mean().reset_index()
        
        # Plot DHW
        plt.bar(yearly_dhw['timestep'], yearly_dhw['dhw'], color='orange', alpha=0.7, label='DHW (Degree Heating Weeks)')
        
        # Plot zero for Cyclones
        plt.plot(yearly_dhw['timestep'], [0]*len(yearly_dhw), color='black', linestyle='--', linewidth=2, label='Cyclone Tracks (Unforced)')
        
        plt.title(reef, fontsize=14, fontweight='bold')
        plt.xlabel('Year', fontsize=12)
        plt.ylabel('DHW (°C-weeks)', fontsize=12)
        plt.xlim(1985, 2024)
        plt.ylim(0, max(10, yearly_dhw['dhw'].max() * 1.2)) # Ensure we can see peaks up to at least 10
        plt.legend(loc='upper right')
        plt.grid(True, linestyle='--', alpha=0.5, axis='y')

plt.tight_layout()
plt.savefig("sandbox/disturbance_plot.png", dpi=300)
print("Disturbance plot saved to sandbox/disturbance_plot.png")
