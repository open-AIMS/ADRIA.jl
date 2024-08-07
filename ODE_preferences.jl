using Preferences, UUIDs
using OrdinaryDiffEq

set_preferences!(UUID("1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"), "PrecompileNonStiff" => true)
set_preferences!(UUID("1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"), "PrecompileStiff" => false)
set_preferences!(
    UUID("1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"), "PrecompileAutoSwitch" => false
)
set_preferences!(
    UUID("1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"), "PrecompileLowStorage" => false
)
set_preferences!(
    UUID("1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"), "PrecompileDefaultSpecialize" => true
)
set_preferences!(
    UUID("1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"), "PrecompileAutoSpecialize" => false
)
set_preferences!(
    UUID("1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"),
    "PrecompileFunctionWrapperSpecialize" => false
)
set_preferences!(
    UUID("1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"), "PrecompileNoSpecialize" => false
)
