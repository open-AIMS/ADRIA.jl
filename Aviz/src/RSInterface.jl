struct RSInterface
    ResultSet
end


function get_components(rs::RSInterface)
    return rs.ResultSet.model_spec[:components]
end

function get_bounds(rs::RSInterface)
    return rs.ResultSet.model_spec[:full_bounds]
end
