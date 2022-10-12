"""
┌──────────┐  ┌───────────────────────────────────────────┐ ┌──────┐
│          │  │                                           │ │      │
│          │  │                                           │ │  H   │
│          │  │                                           │ │  i   │
│          │  │   Trajectories with control sliders       │ │  s   │
│          │  │                                           │ │  t   │
│          │  │                                           │ │      │
│          │  └───────────────────────────────────────────┘ └──────┘
│   C      │
│   o      │  ┌───────────────────────┐ ┌─────────────────┐ ┌───────┐
│   n      │  │                       │ │                 │ │       │
│   t      │  │                       │ │                 │ │ Econ  │
│   r      │  │                       │ │                 │ │       │
│   o      │  │                       │ │  Outcome Matrix │ │       │
│   l      │  │                       │ │                 │ └───────┘
│   s      │  │    Map                │ │                 │
│          │  │                       │ │                 │ ┌───────┐
│          │  │                       │ │                 │ │       │
│          │  │                       │ └─────────────────┘ │ Econ  │
│          │  │                       │ ┌─────────────────┐ │control│
│          │  │                       │ │   Message       │ │       │
└──────────┘  └───────────────────────┘ └─────────────────┘ └───────┘
"""
function comms_layout(; resolution=(1920, 1080))
    f = Figure(resolution=resolution)

    main = f[1:5, 1:5] = GridLayout()

    controls = main[1:5, 1] = GridLayout()
    controls.width = 200

    map = main[3:5, 2] = GridLayout()

    # Trajectories and density plot
    trajectory = main[1:2, 2:4] = GridLayout()
    scen_hist = Axis(
        trajectory[1,4]
    )

    temporal = Axis(
        trajectory[1, 2:3],
        title="Scenario Trajectories",
        xlabel="Year",
        ylabel="Mean TAC (m²)"
    )

    # Show y-axis in millions
    temporal.ytickformat = xs -> ["$(x/1e6)M" for x in xs]

    # Time slider for trajectory
    traj_outcome_sld = trajectory[1, 1]
    traj_time_sld = trajectory[2, 2:3]

    # Outcome probabilities
    outcome_view = trajectory[1, 5]
    outcomes = Axis(
        outcome_view,
        title="Outcomes",
        xlabel="Probability"
    )

    # Importance
    # feat_importance = Axis(
    #     main[3:4, 3],
    #     title="Relative Importance (Top 10)"
    # )
    feat_importance = main[3:4, 3]

    # Economics
    econ_disp = Axis(main[3:4, 4])
    econ_ctrl = main[5, 4]

    return (figure=f,
            controls=controls,
            trajectory=(temporal=temporal, outcome_slider=traj_outcome_sld, time_slider=traj_time_sld),
            scen_hist=scen_hist,
            map=map[1,1],
            outcomes=outcomes,
            importance=feat_importance,
            econ_view=econ_disp,
            econ_ctrl=econ_ctrl)
end



"""
┌──┐ ┌─────────────────────────┐ ┌─────┐  ┌──────────────┐
│  │ │                         │ │ H   │  │              │
│  │ │                         │ │ i   │  │              │
│  │ │ Trajectories            │ │ s   │  │ Map          │
│  │ │                         │ │ t   │  │              │
│  │ │                         │ │     │  │              │
└──┘ └─────────────────────────┘ └─────┘  └──────────────┘
     ┌─────────────────────────┐
     │  sliders                │
     └─────────────────────────┘
     ┌───────────────────────────────────────────────────┐
     │                                                   │
     │          Parallel Coordinates                     │
     │                                                   │
     └───────────────────────────────────────────────────┘
     ┌────────────────────────────┐  ┌───────────────────┐
     │                            │  │                   │
     │                            │  │                   │
     │   Pairplot                 │  │   PCP of outcomes │
     │                            │  │                   │
     │                            │  │                   │
     └────────────────────────────┘  └───────────────────┘
"""
function modeler_layout(; resolution=(1920, 1080))
    f = Figure(resolution=resolution)

    # controls = f[1:3, 1] = GridLayout()
    main = f[1:3, 1:6] = GridLayout()

    spatial_temporal = main[1, 1:6] = GridLayout()

    # Time slider for trajectory
    traj_outcome_sld = spatial_temporal[1,1]

    temporal = Axis(
        spatial_temporal[1, 2:4],
        title="Scenario Trajectories",
        xlabel="Year",
        ylabel="Mean TAC (m²)"
    )
    # Show y-axis in millions
    temporal.ytickformat = xs -> ["$(x/1e6)M" for x in xs]

    traj_time_sld = spatial_temporal[2,1:4]

    scen_hist = Axis(
        spatial_temporal[1,5]
    )
    colgap!(spatial_temporal, 5)

    spatial = spatial_temporal[1, 6]

    interv_pcp = Axis(
        main[2, 1:6],
        title="Interventions"
    )

    outcomes = main[3, 1:6] = GridLayout()
    pairplot = outcomes[1:3,1:3] = GridLayout()
    outcome_pcp = Axis(
        outcomes[1:3, 4:6],
        title="Outcomes"
    )

    return (figure=f,
            # controls=controls,
            trajectory=(temporal=temporal, outcome_slider=traj_outcome_sld, time_slider=traj_time_sld),
            scen_hist=scen_hist, map=spatial,
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
