"""
    screen_scenarios(y::AbstractArray, rule)

Return the indices of rows (scenarios) from `y` whose every column value
passes `rule` after column-wise normalization.

# Arguments
- `y`: an AbstractArray with scenarios in rows and features/variables in columns.
- `rule`: a predicate applied elementwise to the column-normalized data.
  The predicate must accept a scalar and return a Bool (for example `x -> x > 0.1`).
  A row is selected only when `rule` is true for every element in that row.

# Returns
A vector of integer row indices indicating which scenarios passed the rule
on all columns.

# Notes
- Columns are normalized by `col_normalize(copy(y))` before applying `rule`.
- If `y` contains missing values or `rule` does not return Bool for scalars,
  the function may error.

# Example
```julia
y = [1.0 2.0; 3.0 4.0]
screen_scenarios(y, x -> x >= 0.25)  # returns indices of rows meeting the condition
```
"""
function screen_scenarios(y::AbstractArray, rule)
    y_star = col_normalize(copy(y))

    return findall(all.(eachrow(map(rule, y_star))))
end
