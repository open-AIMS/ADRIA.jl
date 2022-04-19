"""
The scope that different coral groups and size classes have for 
producing larvae without consideration of environment.

Coral fecundity per coral area of the different size classes.  
When multiplied by the relative cover of each size class within taxa,
this produces an estimate of the relative fecundity of each coral group and size.  
Total relative fecundity of a group is then calculated as the sum of 
fecundities across size classes. 

Parameters
----------
fec_groups : Matrix[n_classes, n_sites], memory cache to place results into
fec_all : Matrix[n_taxa, n_sites], temporary cache to place intermediate fecundity values into
fec_params : Vector, coral fecundity parameters
Y_pstep : Matrix[n_taxa, n_sites], of values in previous time step
site_area : Vector[n_sites], of site areas

Returns
-------
Matrix[n_classes, n_sites] : fecundity per m2 of coral
"""
function fecundity_scope!(fec_groups::Array{Float64, 2}, fec_all::Array{Float64, 2}, fec_params::Array{Float64}, 
                          Y_pstep::Array{Float64, 2}, site_area::Array{Float64})::Nothing
    ngroups::Int64 = size(fec_groups, 1)   # number of coral groups: 6
    nclasses::Int64 = size(fec_params, 1)  # number of coral size classes: 36

    fec_all .= fec_params .* Y_pstep .* site_area;
    for (i, (s, e)) in enumerate(zip(1:ngroups:nclasses, ngroups:ngroups:nclasses+1))
        @views fec_groups[i, :] = sum(fec_all[s:e, :], dims=1)
    end

    # Above is equivalent to the below, but generic to any group/class size
    # @views fec_groups[1, :] = sum(fec_all[1:6, :], dims=1);   # Tabular Acropora enhanced
    # @views fec_groups[2, :] = sum(fec_all[7:12, :], dims=1);  # Tabular Acropora unenhanced
    # @views fec_groups[3, :] = sum(fec_all[13:18, :], dims=1); # Corymbose Acropora enhanced
    # @views fec_groups[4, :] = sum(fec_all[19:24, :], dims=1); # Corymbose Acropora unenhanced
    # @views fec_groups[5, :] = sum(fec_all[25:30, :], dims=1); # Small massives and encrusting
    # @views fec_groups[6, :] = sum(fec_all[31:36, :], dims=1); # Large massives

    return nothing
end