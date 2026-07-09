import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def validation_metrics(sim_series, obs_series):
    paired = pd.concat([sim_series.rename('sim'), obs_series.rename('obs')], axis=1).dropna()
    n = len(paired)
    if n == 0:
        return {'n': 0, 'pearson': np.nan, 'spearman': np.nan, 'rmse': np.nan, 'pbias': np.nan}

    sim = paired['sim'].astype(float)
    obs = paired['obs'].astype(float)
    pearson = np.corrcoef(sim, obs)[0, 1] if n >= 2 and sim.std() > 0 and obs.std() > 0 else np.nan
    sim_rank = sim.rank(method='average')
    obs_rank = obs.rank(method='average')
    spearman = np.corrcoef(sim_rank, obs_rank)[0, 1] if n >= 2 and sim_rank.std() > 0 and obs_rank.std() > 0 else np.nan
    rmse = np.sqrt(np.mean((sim - obs) ** 2))
    obs_sum = obs.sum()
    pbias = 100.0 * (sim.sum() - obs_sum) / obs_sum if obs_sum != 0 else np.nan
    return {'n': n, 'pearson': pearson, 'spearman': spearman, 'rmse': rmse, 'pbias': pbias}


def fmt_metric(value, digits=2, suffix=''):
    if pd.isna(value):
        return 'NA'
    return f'{value:.{digits}f}{suffix}'


def fmt_metrics(metrics):
    return (
        f"COTS fit (n={metrics['n']}): "
        f"r={fmt_metric(metrics['pearson'])} | "
        f"rho={fmt_metric(metrics['spearman'])} | "
        f"RMSE={fmt_metric(metrics['rmse'])} | "
        f"PBIAS={fmt_metric(metrics['pbias'], 1, '%')}"
    )


def summarize_validation(reef, metric_name, sim_by_year, obs_by_year):
    merged = pd.merge(sim_by_year, obs_by_year, on='year', how='inner')
    metrics = validation_metrics(merged['sim'], merged['obs'])
    return {
        'reef_name': reef,
        'metric_target': metric_name,
        'n_obs': metrics['n'],
        'pearson': metrics['pearson'],
        'spearman': metrics['spearman'],
        'rmse': metrics['rmse'],
        'percent_bias': metrics['pbias']
    }, metrics


sim_df = pd.read_csv('sandbox/data/calibration_results_sim.csv')
emp_df = pd.read_csv('sandbox/data/reef_cots.csv')
emp_coral_df = pd.read_csv('sandbox/data/reef_manta.csv')

# Map simulation reef names to empirical AIMS reef names
sim_to_emp_map = {
    'Lizard Island Reef': 'Lizard Isles',
    'MacGillivray Reef': 'Macgillivray Reef',
    'North Direction Reef': 'North Direction Island',
    'Eyrie Reef': 'Eyrie Reef'
}

reefs_to_plot = ['Lizard Island Reef', 'MacGillivray Reef', 'North Direction Reef', 'Eyrie Reef']
validation_rows = []

plt.figure(figsize=(16, 13))

for i, reef in enumerate(reefs_to_plot, 1):
    ax1 = plt.subplot(2, 2, i)
    ax2 = ax1.twinx()
    cots_metric_text = 'COTS fit: not enough matched observations'

    # Sim COTS
    sim_reef = sim_df[sim_df['reef_name'] == reef]
    if not sim_reef.empty:
        ax1.plot(sim_reef['year'], sim_reef['sim_cots_norm'], label='Simulated COTS (Norm)', color='red', linewidth=2)

    # Emp COTS
    emp_name = sim_to_emp_map[reef]
    emp_reef = emp_df[emp_df['reef_name'] == emp_name].copy()

    if not emp_reef.empty:
        max_cots = emp_reef['cotsptow'].max()
        if max_cots > 0:
            emp_reef['obs_norm'] = emp_reef['cotsptow'] / max_cots
            emp_reef_sorted = emp_reef.sort_values(by='year')
            ax1.plot(
                emp_reef_sorted['year'],
                emp_reef_sorted['obs_norm'],
                label='Observed COTS (Norm)',
                color='lightcoral',
                marker='o',
                linestyle='--',
                linewidth=1.5,
                markerfacecolor='darkred',
                markeredgecolor='darkred',
                zorder=5
            )

            if not sim_reef.empty:
                sim_cots_metric = sim_reef[['year', 'sim_cots_norm']].rename(columns={'sim_cots_norm': 'sim'})
                obs_cots_metric = emp_reef.groupby('year', as_index=False)['obs_norm'].mean().rename(columns={'obs_norm': 'obs'})
                row, metrics = summarize_validation(reef, 'cots_norm', sim_cots_metric, obs_cots_metric)
                validation_rows.append(row)
                cots_metric_text = fmt_metrics(metrics)

    # Sim Coral
    if not sim_reef.empty:
        ax2.plot(sim_reef['year'], sim_reef['sim_coral'], label='Simulated Coral Cover', color='blue', linewidth=2)

    # Emp Coral
    emp_coral_reef = emp_coral_df[emp_coral_df['reef_name'] == emp_name].copy()
    if not emp_coral_reef.empty:
        emp_coral_sorted = emp_coral_reef.sort_values(by='report_year')
        ax2.plot(
            emp_coral_sorted['report_year'],
            emp_coral_sorted['mean'],
            label='Observed Coral Cover',
            color='cornflowerblue',
            marker='s',
            linestyle='--',
            linewidth=1.5,
            markerfacecolor='darkblue',
            markeredgecolor='darkblue',
            zorder=5
        )

        if not sim_reef.empty:
            sim_coral_metric = sim_reef[['year', 'sim_coral']].rename(columns={'sim_coral': 'sim'})
            obs_coral_metric = emp_coral_reef.groupby('report_year', as_index=False)['mean'].mean().rename(columns={'report_year': 'year', 'mean': 'obs'})
            row, _ = summarize_validation(reef, 'coral_cover', sim_coral_metric, obs_coral_metric)
            validation_rows.append(row)

    ax1.set_title(f'{reef}\n{cots_metric_text}', fontsize=12, fontweight='bold')
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

if validation_rows:
    pd.DataFrame(validation_rows).to_csv('sandbox/data/calibration_validation_metrics.csv', index=False)

plt.tight_layout()
plt.savefig('sandbox/calibration_plot.png', dpi=300)
print('Plot saved to sandbox/calibration_plot.png')
if validation_rows:
    print('Validation metrics saved to sandbox/data/calibration_validation_metrics.csv')