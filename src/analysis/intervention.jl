using ADRIA.decision
using ADRIA: relative_leftover_space, connectivity_strength
"""
	intervention_frequency(rs::ResultSet, scen_indices::NamedTuple, log_type::Symbol)::NamedDimsArray

Count number of times a location of selected for intervention
Count frequency of seeded sites for scenarios satisfying a condition.

# Arguments
- 'rs' : ResultSet
- `scen_indices` : rcp_id => scenario id that satisfy a condition of interest.
- 'log_type` : the intervention log to use in calculating frequencies (one of :seed, :shade or :fog).

# Returns
NamedDimsArray(:locations, :rcps)

# Example
```julia
using ADRIA, Statistics

rs = ADRIA.load_results("some result set")

tac = ADRIA.metrics.scenario_total_cover(rs)
rsv = ADRIA.metrics.scenario_rsv(rs)

# Create matrix of mean scenario outcomes
mean_tac = vec(mean(tac, dims=1))
mean_sv = vec(mean(rsv, dims=1))
y = hcat(mean_tac, mean_sv)

# Find all pareto optimal scenarios where all metrics >= 0.9
rule_func = x -> all(x .>= 0.9)
robust_scens = ADRIA.analysis.find_robust(rs, y, rule_func, [45, 60])

# Retrieve seeding intervention frequency for robust scenarios
robust_selection_frequencies = ADRIA.analysis.intervention_frequency(rs, robust_scens, :seed)
"""
function intervention_frequency(
	rs::ResultSet, scen_indices::NamedTuple, log_type::Symbol
)::NamedDimsArray
	log_type ∈ [:seed, :shade, :fog] || ArgumentError("Unsupported log")

	# Get requested log
	interv_log = getfield(rs, Symbol("$(log_type)_log"))
	rcps = collect(Symbol.(keys(scen_indices)))
	n_locs = n_locations(rs)

	interv_freq = NamedDimsArray(
		zeros(n_locs, length(rcps)); locations = rs.site_ids, rcps = rcps
	)
	for rcp in rcps
		# Select scenarios satisfying condition and tally selection for each location
		logged_data = dropdims(
			sum(interv_log[scenarios = scen_indices[rcp]]; dims = :coral_id);
			dims = :coral_id,
		)
		interv_freq(; rcps = rcp) .= vec(
			dropdims(
				sum(logged_data .> 0; dims = [:timesteps, :scenarios]); dims = :timesteps
			),
		)
	end

	return interv_freq
end

"""
	decision_matrices(rs::ResultSet, criteria_weights::DataFrameRow;
		loc_coral_cover = rs.site_max_coral_cover::Vector{Float64},
		area_to_seed::Float64 = 962.11)
		
Calculates a decision matrix for a specified intervention, using a scenario specification 
	and Domain alone. These can be visualised spatially using `viz.decision_matrices`.

# Arguments
- `rs` : ADRIA ResultSet
- `criteria_weights` :  A row of a scenario dataframe, containing intervention criteria weights.
- `int_type` : Intervention type (e.g. :seed or :fog)
- `loc_coral_cover` : Relative coral cover to site k area (dims: nspecies*nsites), default 
	is max cover over scenarios in rs.
- `area_to_seed` : Area of corals seeded in the scenario considered

# Returns
Selection score
"""
function decision_matrices(
	rs::ResultSet,
	criteria_weights::DataFrameRow;
	loc_coral_cover = rs.site_max_coral_cover::Vector{Float64},
	area_to_seed::Float64 = 962.11,
	RCP::String = "45",
)
	site_ids = collect(1:length(rs.site_data.site_id))
	leftover_space = relative_leftover_space(loc_coral_cover) .* site_k_area(rs)

	wave_stress = ADRIA.decision.summary_stat_env(
		rs.wave_stats[RCP](; stat = "mean"), (:wave_scenario)
	)
	heat_stress = ADRIA.decision.summary_stat_env(
		rs.dhw_stats[RCP](; stat = "mean"), (:dhw_scenario)
	)
	TP_data = rs.connectivity_data[RCP]
	connectivity_data = connectivity_strength(
		TP_data .* site_k_area(rs), vec(loc_coral_cover), TP_data
	)
	strong_pred = connectivity_data.strongest_predecessor

	predec = ADRIA.decision.priority_predecessor_criteria(
		strong_pred, vec(rs.sim_constants["priority_sites"]), length(site_ids)
	)
	zones_crit = ADRIA.decision.zones_criteria(
		vec(rs.sim_constants["priority_zones"]), rs.site_data.zone_type, strong_pred,
		site_ids,
	)

	A, filtered_sites = ADRIA.decision.create_decision_matrix(
		site_ids,
		connectivity_data.in_conn,
		connectivity_data.out_conn,
		leftover_space[site_ids],
		wave_stress[site_ids],
		heat_stress[site_ids],
		rs.site_data.depth_med[site_ids],
		predec,
		zones_crit,
		criteria_weights.deployed_coral_risk_tol,
	)

	S, wse = ADRIA.decision.create_seed_matrix(
		A,
		criteria_weights.seed_in_connectivity,
		criteria_weights.seed_out_connectivity,
		criteria_weights.seed_wave_stress,
		criteria_weights.seed_heat_stress,
		criteria_weights.seed_priority,
		criteria_weights.seed_zone,
		criteria_weights.coral_cover_low,
		criteria_weights.seed_depth,
	)

	SE = NamedDimsArray(
		S[:, 2:end];
		locations = rs.site_data.site_id,
		criteria = [
			:in_connectivity,
			:out_connectivity,
			:wave_stress,
			:heat_stress,
			:priority_predecessor,
			:zones,
			:leftover_space,
			:depth,
		],
	)

	S, wsh = ADRIA.decision.create_fog_matrix(
		A,
		site_k_area(rs)[site_ids][filtered_sites],
		criteria_weights.fog_connectivity,
		criteria_weights.fog_wave_stress,
		criteria_weights.fog_heat_stress,
		criteria_weights.fog_priority,
		criteria_weights.fog_zone,
		criteria_weights.coral_cover_high,
	)
	SH = NamedDimsArray(
		S[:, 2:end];
		locations = rs.site_data.site_id,
		criteria = [
			:in_connectivity,
			:out_connectivity,
			:wave_stress,
			:heat_stress,
			:priority_predecessor,
			:zones,
			:coral_area,
		],
	)

	return SE, wse, SH, wsh
end
