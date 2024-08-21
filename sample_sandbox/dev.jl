using WGLMakie, GeoMakie, GraphMakie, Bonito
using ADRIA

Page(exportable=true, offline=true)


dom = ADRIA.load_domain("/data/input/Moore_2024-02-05_v050_rc2")

scens = ADRIA.sample(dom, 128)
rs = ADRIA.run_scenarios(dom, scens, "45")

s_tac = ADRIA.metrics.scenario_total_cover(rs)
f = ADRIA.viz.scenarios(rs, s_tac; fig_opts=Dict{Any,Any}(:size=>(1600, 800)))
save("/data/output/adria_scenario_outcomes.html", f)
# save("/data/output/adria_scenario_outcomes.png", f)