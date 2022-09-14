"""
┌──────────┐ ┌────────────────┐ ┌────┐ ┌──────────┐
│          │ │                │ │ H  │ │          │
│          │ │  Trajectories  │ │ i  │ │  Map     │
│          │ │                │ │ s  │ │          │
│          │ │                │ │ t  │ │          │
│          │ └────────────────┘ └────┘ └──────────┘
│ Controls │
│          │ ┌────────────────────────────────────┐
│          │ │                                    │
│          │ │      Parallel Coordinates          │
│          │ │                                    │
│          │ └────────────────────────────────────┘
│          │
│          │ ┌──────────────┐ ┌───────────────────┐
│          │ │              │ │                   │
│          │ │              │ │ PCP of Outcomes   │
│          │ │   Pairplot   │ │                   │
│          │ │              │ │                   │
│          │ │              │ │                   │
└──────────┘ └──────────────┘ └───────────────────┘
"""
function create_layout(; resolution=(1920, 1080))
    f = Figure(resolution=resolution)

    controls = f[1:3, 1] = GridLayout()
    main = f[1:3, 2:5] = GridLayout()

    spatial_temporal = main[1, 1:5] = GridLayout()
    temporal = Axis(
        spatial_temporal[1, 1:3],
        title="Scenario Trajectories",
        xlabel="Year",
        ylabel="Mean TAC (m²)"
    )
    # Show y-axis in millions
    temporal.ytickformat = xs -> ["$(x/1e6)M" for x in xs]

    scen_hist = Axis(
        spatial_temporal[1,4]
    )
    colgap!(spatial_temporal, 5)

    spatial = Axis(
        spatial_temporal[1, 5],
        title="Map",
        xlabel="Long",
        ylabel="Lat",
        # backgroundcolor=:lightskyblue1
    )

    interv_pcp = Axis(
        main[2, 1:5],
        title="Interventions"
    )

    outcomes = main[3, 1:5] = GridLayout()
    pairplot = outcomes[1:3,1:2] = GridLayout()
    outcome_pcp = Axis(
        outcomes[1:3, 3:5],
        title="Outcomes"
    )

    return (figure=f,
            controls=controls,
            trajectory=temporal, scen_hist=scen_hist, map=spatial,
            interv_pcp=interv_pcp, 
            pairplot=pairplot, outcome_pcp=outcome_pcp)
end


# """
# ┌─────────────┐  ┌─────────────────────┐ ┌────────┐
# │             │  │                     │ │        │
# │             │  │     Trajectories    │ │  Map   │
# │             │  │                     │ │        │
# │             │  └─────────────────────┘ └────────┘
# │             │
# │  Controls   │  ┌────────────────────────────────┐
# │             │  │                                │
# │             │  │      Parallel Coordinates      │
# │             │  │                                │
# │             │  └────────────────────────────────┘
# │             │
# │             │  ┌──────────────┐ ┌───────────────┐
# │             │  │              │ │    Violin     │
# │             │  │              │ └───────────────┘
# │             │  │   Pairplot   │
# │             │  │              │ ┌───────────────┐
# │             │  │              │ │   Histogram   │
# └─────────────┘  └──────────────┘ └───────────────┘
# """
# function create_layout(; resolution=(1600, 900))
#     f = Figure(resolution=resolution)

#     controls = f[1:3, 1] = GridLayout()
#     main = f[1:3, 2] = GridLayout()

#     spatial_temporal = main[1, 1:2] = GridLayout()
#     temporal = Axis(
#         spatial_temporal[1, 1:2],
#         title="Scenario Trajectories",
#         xlabel="Year",
#         ylabel="Mean TAC (m²)"
#     )
#     # Show y-axis in millions
#     temporal.ytickformat = xs -> ["$(x/1e6)M" for x in xs]

#     spatial = Axis(
#         spatial_temporal[1, 3],
#         title="Map",
#         xlabel="Long",
#         ylabel="Lat",
#         # backgroundcolor=:lightskyblue1
#     )

#     pcp = Axis(
#         main[2, 1:2],
#         title="Outcomes"
#     )

#     site_spread = main[3, 1:2] = GridLayout()
#     site_violin = Axis(
#         site_spread[1,2],
#         title="Spread"
#     )

#     site_hist = Axis(
#         site_spread[2:3,2]
#     )

#     pairplot = site_spread[1:3,1] = GridLayout()

#     return (figure=f, controls=controls,
#             trajectory=temporal, map=spatial, pcp=pcp, pairplot=pairplot,
#             violin=site_violin, hist=site_hist)
# end
