# Crown-of-thorns starfish (COTS) submodel: data-package inputs and the adapter to COTSMod.jl.
# The factor declarations (`cots_factors.jl`) are included earlier, with the other model
# factors; this file needs `Domain` and so is included after it.
include("cots_inputs.jl")
include("cots_dynamics.jl")
