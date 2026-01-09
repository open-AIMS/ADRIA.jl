"""
    sample_balanced(d::Domain, n::Int64; sample_method=SobolSample(R=OwenScramble(base=2, pad=32)))::DataFrame

Generate equal samples across counterfactual, unguided, and guided scenarios.

# Arguments
- `d` : Domain.
- `n` : number of samples to generate *per intervention type*.
- `sample_method` : type of sampler to use.

# Returns
Scenario specification
"""
function sample_balanced(
    d::Domain, n::Int64; sample_method=SobolSample(; R=OwenScramble(; base=2, pad=32))
)::DataFrame
    cf = sample_cf(d, n; sample_method)
    ug = sample_unguided(d, n; sample_method)
    gd = sample_guided(d, n; sample_method)

    return vcat(cf, ug, gd)
end
