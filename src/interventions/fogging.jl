function fog_locations!(Yfog, locs, dhw_t, fogging)
    dhw_t[locs] .= dhw_t[locs] .* (1.0 .- fogging)
    return Yfog[locs] .= fogging
end
