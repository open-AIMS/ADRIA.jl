import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from pathlib import Path


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


# Load stochastic simulation runs (fixed best COTS parameters)
output_tag = os.environ.get('COTS_OUTPUT_TAG', '').strip()
output_prefix = 'best_stochastic' if output_tag == '' else f'best_stochastic_{output_tag}'
sim_df = pd.read_csv(f'sandbox/data/{output_prefix}_trajectories.csv')
site_path = Path(f'sandbox/data/{output_prefix}_site_trajectories.csv')
site_df = pd.read_csv(site_path) if site_path.exists() else pd.DataFrame()
meta_path = Path(f'sandbox/data/{output_prefix}_metadata.csv')
meta_df = pd.read_csv(meta_path) if meta_path.exists() else pd.DataFrame()

# Load empirical data
emp_df = pd.read_csv('sandbox/data/reef_cots.csv')
emp_coral_df = pd.read_csv('sandbox/data/reef_manta.csv')

sim_to_emp_map = {
    'Lizard Island Reef': 'Lizard Isles',
    'MacGillivray Reef': 'Macgillivray Reef',
    'North Direction Reef': 'North Direction Island',
    'Eyrie Reef': 'Eyrie Reef'
}

reefs_to_plot = ['Lizard Island Reef', 'MacGillivray Reef', 'North Direction Reef', 'Eyrie Reef']
validation_rows = []

plt.figure(figsize=(16, 13))
n_runs = sim_df['sim_id'].nunique()
stochastic_mode = meta_df['stochastic_mode'].iloc[0] if not meta_df.empty and 'stochastic_mode' in meta_df else 'unknown'

for i, reef in enumerate(reefs_to_plot, 1):
    ax1 = plt.subplot(2, 2, i)
    ax2 = ax1.twinx()
    cots_ylim_top = 1.05
    cots_metric_text = 'COTS fit: not enough matched observations'

    reef_sim = sim_df[sim_df['reef_name'] == reef]

    if not reef_sim.empty:
        reef_cots_scale = reef_sim['sim_cots_adult'].max() if 'sim_cots_adult' in reef_sim else 0.0

        # 1. Plot individual stochastic simulation runs in light alpha to show envelope
        for sim_id in reef_sim['sim_id'].unique()[:15]:
            run_data = reef_sim[reef_sim['sim_id'] == sim_id]
            ax1.plot(run_data['year'], run_data['sim_cots_norm'], color='red', alpha=0.12, linewidth=1.0)
            ax2.plot(run_data['year'], run_data['sim_coral'], color='blue', alpha=0.12, linewidth=1.0)

        # 2. Calculate median and 10th-90th percentile uncertainty bounds across runs.
        summary_cots = reef_sim.groupby('year')['sim_cots_norm'].agg([
            'median',
            lambda x: np.percentile(x, 10),
            lambda x: np.percentile(x, 90)
        ]).reset_index()
        summary_cots.columns = ['year', 'median', 'p10', 'p90']
        cots_ylim_top = max(cots_ylim_top, summary_cots['p90'].max() * 1.05)

        summary_coral = reef_sim.groupby('year')['sim_coral'].agg([
            'median',
            lambda x: np.percentile(x, 10),
            lambda x: np.percentile(x, 90)
        ]).reset_index()
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
                if reef_cots_scale > 0:
                    reef_sites['sim_cots_norm'] = reef_sites['sim_cots_adult'] / reef_cots_scale
                    site_cots = reef_sites.groupby('year')['sim_cots_norm'].agg([
                        lambda x: np.percentile(x, 10),
                        lambda x: np.percentile(x, 90)
                    ]).reset_index()
                    site_cots.columns = ['year', 'p10', 'p90']
                    cots_ylim_top = max(cots_ylim_top, site_cots['p90'].max() * 1.05)
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

    # 4. Plot empirical observations and calculate matched-year validation metrics.
    emp_name = sim_to_emp_map[reef]

    emp_reef = emp_df[emp_df['reef_name'] == emp_name].copy()
    if not emp_reef.empty:
        max_cots = emp_reef['cotsptow'].max()
        if max_cots > 0:
            emp_reef['obs_norm'] = emp_reef['cotsptow'] / max_cots
            ax1.scatter(emp_reef['year'], emp_reef['obs_norm'], label='Obs COTS (cots/tow norm)', color='darkred', marker='o', s=65, zorder=5, edgecolor='black', linewidth=0.5)

            if not reef_sim.empty:
                sim_cots_metric = summary_cots[['year', 'median']].rename(columns={'median': 'sim'})
                obs_cots_metric = emp_reef.groupby('year', as_index=False)['obs_norm'].mean().rename(columns={'obs_norm': 'obs'})
                row, metrics = summarize_validation(reef, 'cots_norm', sim_cots_metric, obs_cots_metric)
                validation_rows.append(row)
                cots_metric_text = fmt_metrics(metrics)

    emp_coral_reef = emp_coral_df[emp_coral_df['reef_name'] == emp_name].copy()
    if not emp_coral_reef.empty:
        ax2.scatter(emp_coral_reef['report_year'], emp_coral_reef['mean'], label='Obs Coral Cover', color='darkblue', marker='s', s=65, zorder=5, edgecolor='black', linewidth=0.5)

        if not reef_sim.empty:
            sim_coral_metric = summary_coral[['year', 'median']].rename(columns={'median': 'sim'})
            obs_coral_metric = emp_coral_reef.groupby('report_year', as_index=False)['mean'].mean().rename(columns={'report_year': 'year', 'mean': 'obs'})
            row, _ = summarize_validation(reef, 'coral_cover', sim_coral_metric, obs_coral_metric)
            validation_rows.append(row)

    # Formatting
    ax1.set_title(
        f'{reef}\nBest Candidate: run_id=2, Runs={n_runs}, Mode={stochastic_mode}\n{cots_metric_text}',
        fontsize=11,
        fontweight='bold'
    )
    ax1.set_xlabel('Year', fontsize=12)
    ax1.set_ylabel('Normalized Adult COTS Density', color='red', fontsize=12)
    ax2.set_ylabel('Total Coral Cover (Proportion)', color='blue', fontsize=12)

    ax1.set_xlim(1985, 2024)
    ax1.set_ylim(-0.05, cots_ylim_top)
    ax2.set_ylim(-0.05, 1.05)

    ax1.grid(True, linestyle='--', alpha=0.6)

    # Combine legends
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    by_label = dict(zip(labels1 + labels2, lines1 + lines2))
    ax1.legend(by_label.values(), by_label.keys(), loc='upper left', fontsize=9, framealpha=0.85)

if validation_rows:
    pd.DataFrame(validation_rows).to_csv(f'sandbox/data/{output_prefix}_validation_metrics.csv', index=False)

plt.tight_layout()
plt.savefig(f'sandbox/{output_prefix}_plot.png', dpi=300)
print(f'Stochastic best-fit plot saved to sandbox/{output_prefix}_plot.png')
if validation_rows:
    print(f'Validation metrics saved to sandbox/data/{output_prefix}_validation_metrics.csv')