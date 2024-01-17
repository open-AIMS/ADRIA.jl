using ADRIA.decision
using ADRIA: relative_leftover_space, connectivity_strength, model_spec, criteria_params,
    component_params

"""
    decision_matrices(rs::ResultSet, criteria_row::DataFrameRow;
        loc_coral_cover = rs.site_max_coral_cover::Vector{Float64}, RCP::String = "45")

Calculates a decision matrix for a specified intervention, using a scenario specification 
	and Domain alone. These can be visualised spatially using `viz.decision_matrices`.

# Arguments
- `rs` : ADRIA ResultSet
- `criteria_row` :  A row of a scenario dataframe, containing intervention criteria weights.
- `loc_coral_cover` : Relative coral cover to site k area (dims: nspecies*nsites), default 
is max cover over scenarios in rs.
-  `RCP` : RCP scenario to use (waves and dhws), default is "45".

# Returns
Selection score
"""
function decision_matrices(
    rs::ResultSet,
    criteria_row::DataFrameRow;
    loc_coral_cover = rs.site_max_coral_cover::Vector{Float64},
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

    TP_data = rs.connectivity_data[RCP]  # connectivity matrix
    connectivity_data = connectivity_strength(
        TP_data .* site_k_area(rs), vec(loc_coral_cover), TP_data
    )
    # Strongest larval source location for each location
    strong_pred = connectivity_data.strongest_predecessor

    # Criteria for strongest larval sources to priority locations
    priority_source_criteria = ADRIA.decision.priority_predecessor_criteria(
        strong_pred, vec(rs.sim_constants["priority_sites"]), length(site_ids)
    )
    # Criteria for strongest larval sources to/members of priority zones
    zones_criteria = ADRIA.decision.zones_criteria(
        vec(rs.sim_constants["priority_zones"]), rs.site_data.zone_type, strong_pred,
        site_ids,
    )

    spec = model_spec(rs)
    crit = component_params(spec, CriteriaWeights)
    weights_seed_crit = criteria_params(crit, "(:seed, :weight)")
    weights_fog_crit = criteria_params(crit, "(:fog, :weight)")

    A, filtered_sites = ADRIA.decision.create_decision_matrix(
        site_ids,
        connectivity_data.in_conn,
        connectivity_data.out_conn,
        leftover_space[site_ids],
        wave_stress[site_ids],
        heat_stress[site_ids],
        rs.site_data.depth_med[site_ids],
        priority_source_criteria,
        zones_criteria,
        criteria_row.deployed_coral_risk_tol,
    )

    S, wse = ADRIA.decision.create_seed_matrix(
        A,
        criteria_row.seed_in_connectivity,
        criteria_row.seed_out_connectivity,
        criteria_row.seed_wave_stress,
        criteria_row.seed_heat_stress,
        criteria_row.seed_priority,
        criteria_row.seed_zone,
        criteria_row.seed_coral_cover_low,
        criteria_row.seed_depth;
        filter_space = -1.0,
    )
    SE = NamedDimsArray(
        S[:, 2:end];
        locations = rs.site_data.site_id[Int.(S[:, 1])],
        criteria = weights_seed_crit.fieldname,
    )

    S, wsh = ADRIA.decision.create_fog_matrix(
        A,
        site_k_area(rs)[site_ids][filtered_sites],
        criteria_row.fog_in_connectivity,
        criteria_row.fog_out_connectivity,
        criteria_row.fog_wave_stress,
        criteria_row.fog_heat_stress,
        criteria_row.fog_priority,
        criteria_row.fog_zone,
        criteria_row.fog_coral_cover_high,
    )
    SH = NamedDimsArray(
        S[:, 2:end];
        locations = rs.site_data.site_id[Int.(S[:, 1])],
        criteria = weights_fog_crit.fieldname,
    )

	return SE, wse, SH, wsh
end
