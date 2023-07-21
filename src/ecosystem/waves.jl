"""
    wave_damage!(Sw_t::Matrix{T}, wave_scen::Array{T, 3}, wave_mortality::Vector{T}, 
        n_species::Int64)::Array{Float64,3} where {T<:Float64}

Calculate wave damage for each species/group using 90th percentile wave mortality data.

# Arguments
- `Sw_t` : cache matrix
- `wave_scen` : wave stress trajectory for all time steps and locations
- `wave_mortality` : 90th percentile wave mortality data
- `n_species` : number of species represented

# Returns
Proportion of corals that survive wave stress
"""
function wave_damage!(Sw_t::Matrix{T}, wave_scen::Array{T,3}, wave_mortality::Vector{T},
    n_species::Int64)::Array{Float64,3} where {T<:Float64}

    for sp::Int64 in 1:n_species
        @views Sw_t[:, sp, :] .= wave_mortality[sp] .* wave_scen
    end

    clamp!(Sw_t, 0.0, 1.0)

    # Wave damage survival
    Sw_t .= 1.0 .- Sw_t

    return Sw_t
end
