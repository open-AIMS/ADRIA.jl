using Preferences, UUIDs
using OrdinaryDiffEq

ode_uuid = UUID("1dea7af3-3e70-54e6-95c3-0bf5283fa5ed")

set_preferences!(ode_uuid, "PrecompileNonStiff" => true)
set_preferences!(ode_uuid, "PrecompileStiff" => false)
set_preferences!(ode_uuid, "PrecompileAutoSwitch" => false)
set_preferences!(ode_uuid, "PrecompileLowStorage" => false)
set_preferences!(ode_uuid, "PrecompileDefaultSpecialize" => true)
set_preferences!(ode_uuid, "PrecompileAutoSpecialize" => false)
set_preferences!(ode_uuid, "PrecompileFunctionWrapperSpecialize" => false)
set_preferences!(ode_uuid, "PrecompileNoSpecialize" => false)
