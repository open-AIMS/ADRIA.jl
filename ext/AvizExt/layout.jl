"""


┌─────┐ ┌───────────────────────────────────────────────────────────┐
│     │ │                                                           │
│     │ │                                                           │
│     │ │                                                           │
│  C  │ │                                                           │
│  o  │ │           Trajectories with control sliders               │
│  n  │ │                                                           │
│  t  │ │                                                           │
│  r  │ │                                                           │
│  o  │ └───────────────────────────────────────────────────────────┘
│  l  │ ┌────────────────┐ ┌───────────────────────────┐  ┌─────────┐
│  s  │ │                │ │                           │  │         │
│     │ │                │ │   Relative Importance     │  │ Probas  │
│     │ │                │ │                           │  │         │
│     │ │    Map         │ │                           │  │         │
│     │ │                │ └───────────────────────────┘  └─────────┘
│     │ │                │ ┌────────────────────────────────────────┐
│     │ │                │ │         Message                        │
└─────┘ └────────────────┘ └────────────────────────────────────────┘
"""
function comms_layout(; size=(1920, 1080))
    f = Figure(size=size)

    main = f[1:6, 1:9] = GridLayout()

    controls = main[1:6, 1:2] = GridLayout()
    # controls.width = 400

    # Trajectories and density plot
    trajectory = main[1:2, 3:8] = GridLayout()
    temporal = Axis(
        trajectory[1, 2:7],
        title="Scenario Trajectories",
        xlabel="Year",
        ylabel="Mean TAC (m²)"
    )
    scen_hist = Axis(
        trajectory[1, 8]
    )

    # Show y-axis in millions
    temporal.ytickformat = xs -> ["$(x/1e6)M" for x in xs]

    # Time slider for trajectory
    traj_outcome_sld = trajectory[1, 1]
    traj_time_sld = trajectory[2, 2:7]

    map = main[3:6, 3:5] = GridLayout()

    # Importance
    # feat_importance = Axis(
    #     main[3:4, 3],
    #     title="Relative Importance (Top 10)"
    # )
    feat_importance = main[3:5, 6:7]

    # # Economics
    # econ_disp = Axis(main[3:4, 4])
    # econ_ctrl = main[5, 4]

    # Outcome probabilities
    outcome_view = main[3:5, 8:9]
    outcomes = Axis(
        outcome_view,
        title="Probability Occurrence",
        xlabel="Outcomes",
        xticks=([1, 2, 3, 4, 5],
            ["Very High\n> 80%", "High\n70 - 80%", "Medium\n50 - 70%", "Low\n20 - 50%", "Very Low\n< 20%"])
    )

    messages = Axis(main[6, 5:9])
    hidedecorations!(messages)
    hidespines!(messages)
    text!(messages,
        0.0,
        0.5,
        text="Zoom: Mouse wheel\nPan: Hold right-click\nReset view: Ctrl + Left-click",
        align=(:left, :center),
        justification=:left,
        fontsize=10)

    return (figure=f,
        controls=controls,
        trajectory=(temporal=temporal, outcome_slider=traj_outcome_sld, time_slider=traj_time_sld),
        scen_hist=scen_hist,
        map=map[1, 1],
        importance=feat_importance,
        outcomes=outcomes,
        messages)
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
function modeler_layout(; size=(1920, 1080))
    f = Figure(size=size)

    # controls = f[1:3, 1] = GridLayout()
    main = f[1:3, 1:6] = GridLayout()

    spatial_temporal = main[1, 1:6] = GridLayout()

    # Time slider for trajectory
    traj_outcome_sld = spatial_temporal[1, 1]

    temporal = Axis(
        spatial_temporal[1, 2:4],
        title="Scenario Trajectories",
        xlabel="Year",
        ylabel="Mean TAC (m²)"
    )
    # Show y-axis in millions
    temporal.ytickformat = xs -> ["$(x/1e6)M" for x in xs]

    traj_time_sld = spatial_temporal[2, 1:4]

    scen_hist = Axis(
        spatial_temporal[1, 5]
    )
    colgap!(spatial_temporal, 5)

    spatial = spatial_temporal[1, 6]

    interv_pcp = Axis(
        main[2, 1:6],
        title="Interventions"
    )

    outcomes = main[3, 1:6] = GridLayout()
    pairplot = outcomes[1:3, 1:3] = GridLayout()
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
# function create_layout(; size=(1600, 900))
#     f = Figure(size=size)

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
