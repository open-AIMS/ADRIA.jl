using Test
using ADRIA
using ADRIA.Distributions

@testset "proportional adjustment" begin
    Y = rand(5, 36, 20)
    for i in axes(Y, 1)
        Y[i, :, :] .= Y[i, :, :] / sum(Y[i, :, :]; dims=1)
    end

    tmp = zeros(20)
    for i in axes(Y, 1)
        # No warning should be emitted when values are between 0 and 1
        Test.@test_nowarn ADRIA.proportional_adjustment!(Y[i, :, :], tmp)

        @test all(0.0 .<= Y[i, :, :] .<= 1.0)
    end
end

@testset "Coral Spec" begin
    coral_params = ADRIA.coral_spec().params

    bin_edge_diameters_cm2 = ADRIA.colony_mean_area(ADRIA.bin_edges())
    stored_colony_mean_areas = ADRIA.colony_mean_area(
        coral_params.mean_colony_diameter_m .* 100.0
    )

    # check colony areas in cm^2 are within bounds designated by bin edges
    for k in 1:6
        @test all(
            stored_colony_mean_areas[coral_params.class_id .== k] .>=
            bin_edge_diameters_cm2[k]
        ) ||
            "Some colony areas for size class $k are larger than the size class upper bound."
        @test all(
            stored_colony_mean_areas[coral_params.class_id .== k] .>=
            bin_edge_diameters_cm2[k]
        ) ||
            "Some colony areas for size class $k are smaller than the size class lower bound."
    end
end

@testset "Fecundity" begin
    fec_groups = zeros(6, 216)
    fec_all = zeros(36, 216)
    fec_params = [
        25281.51645394548,
        50989.55542425965,
        78133.52681199001,
        115189.85341730568,
        169820.8550374081,
        250361.6590357049,
        25281.51645394548,
        50989.55542425965,
        78133.52681199001,
        115189.85341730568,
        169820.8550374081,
        250361.6590357049,
        52228.76428701259,
        59199.29777746337,
        63887.49751239493,
        68472.9244216383,
        73387.46329736525,
        78654.73564497223,
        52228.76428701259,
        59199.29777746337,
        63887.49751239493,
        68472.9244216383,
        73387.46329736525,
        78654.73564497223,
        21910.874521191126,
        37082.43894786883,
        51072.30305499843,
        68331.04154366927,
        91421.98332850973,
        122315.9906084096,
        21910.874521191126,
        37082.43894786883,
        51072.30305499843,
        68331.04154366927,
        91421.98332850973,
        122315.9906084096
    ]

    C_cover_t = rand(Uniform(0.0, 0.01), 36, 216)
    total_loc_area = [
        76997.8201778261,
        38180.997513339855,
        334269.6228868989,
        70728.59381575836,
        48824.963081851136,
        87942.62072634231,
        57278.82204914279,
        131481.403591529,
        90463.151137474,
        42261.42923473893,
        312.98931139567867,
        57605.03185068816,
        60083.839003962465,
        54785.65416847123,
        12832.631625673268,
        76044.65694113867,
        100181.29909620434,
        118024.50294493232,
        60109.49596805731,
        242250.00915593235,
        124908.22948251851,
        113635.26297052717,
        91707.8292375924,
        135850.1470950297,
        49141.425121693406,
        53826.22338320641,
        97025.1128987968,
        68525.34328866098,
        148695.41590357665,
        28781.728845587466,
        165585.33163399575,
        23778.652445240412,
        16592.14594898885,
        158322.37248498667,
        118921.10221339483,
        128982.22331511462,
        107034.72890100488,
        86652.49363158084,
        158343.6427825936,
        5318.305293030571,
        9389.681316065602,
        3129.26198370615,
        135152.96035117377,
        23472.247369048186,
        97606.50613648817,
        71946.8830838264,
        35981.50364708854,
        28797.393418124877,
        29107.717398312874,
        53826.99441838125,
        311336.2225115001,
        125505.64010765497,
        99856.55065180548,
        106090.00433640555,
        180018.80202134652,
        326071.049694587,
        190216.44162023207,
        53827.47156010475,
        144629.18991992064,
        148898.01095200004,
        96661.44398395158,
        290148.5026182546,
        114825.04259502981,
        140754.4730709605,
        68829.15950475587,
        95473.48294012994,
        81080.31676690746,
        169308.24664905109,
        114162.37943328498,
        22536.31970276218,
        48824.50477898354,
        64804.19810403744,
        162433.71505506802,
        51000.481316191144,
        150484.32479333598,
        46612.03379469784,
        134619.66478604497,
        54461.06710961368,
        107594.00013558939,
        40370.00313273119,
        62282.677392093,
        111411.61847271444,
        148083.46177229844,
        234284.18705729162,
        96100.27528847847,
        63184.710597992875,
        103282.46208330011,
        126132.27669022558,
        51333.54014409892,
        41937.25823739078,
        70105.24495933158,
        66337.72066151444,
        100498.80730765127,
        22524.106860139407,
        335968.1465102434,
        23157.07392614428,
        64115.71150727989,
        43187.80882960232,
        55396.315229838714,
        322942.3789655925,
        264867.9285628754,
        233662.25014557084,
        134911.29212181736,
        90572.83054631483,
        48411.07756591868,
        87456.35002980288,
        369127.91149331676,
        252347.31258559506,
        231125.33238760522,
        114617.7986012646,
        1561.3605366628617,
        133976.43868495245,
        177710.91558774887,
        261426.5130989817,
        233946.98499754444,
        14987.148259407375,
        68075.72698056,
        69341.32427705498,
        129437.48085331544,
        76901.33963279286,
        111941.78706551343,
        78184.30865436653,
        98454.09477984346,
        52201.226116100326,
        62855.21237831516,
        124458.66966792708,
        24079.841552573256,
        111959.48772720806,
        22512.65185918659,
        74701.63197803684,
        124114.03316707956,
        80338.79890576331,
        41584.86461727973,
        38441.346766835544,
        136971.89531025803,
        167229.85133617045,
        140734.54589663213,
        184158.822707986,
        33770.755155074876,
        17826.207357996143,
        1250.9943127147853,
        101592.19755722815,
        122570.48372718506,
        249215.26396020036,
        183567.23185554985,
        118473.36072853673,
        84114.7206080337,
        252338.2882249197,
        104395.61599875009,
        287106.2325030188,
        248588.7888734066,
        139489.46534616407,
        109694.42252342962,
        226140.37395826913,
        129389.22499938775,
        185781.30174259283,
        106306.2538784002,
        159193.62830397952,
        104134.67320310418,
        86911.49756977474,
        348115.2531043119,
        47815.115320474375,
        190386.1996394787,
        221756.60024294443,
        106927.86914726766,
        89753.67927749828,
        299004.6593301678,
        124114.19568072166,
        120039.92525529955,
        219873.910698622,
        77874.03697757702,
        187571.9804283902,
        58788.913771106396,
        304977.0016628909,
        54074.51778317196,
        75350.34206689568,
        69390.57800343214,
        232402.37505759858,
        126950.81416913401,
        19064.742817895487,
        25021.277749726083,
        14695.997722434346,
        170774.58696733043,
        625.2096516690217,
        130026.698200766,
        205455.53109697672,
        63153.77036182955,
        137544.44021125184,
        107886.94441078696,
        85240.40940979542,
        142395.81966814818,
        60271.87516689906,
        26857.034316257108,
        20922.45744012855,
        226991.3332164348,
        142089.56898094108,
        54014.206533902325,
        144895.9872502829,
        108356.8193304576,
        29666.78814761853,
        27475.359576036688,
        936.6064325589687,
        20608.68716322258,
        47156.65406070184,
        70263.70212964155,
        122069.65583620407,
        9989.258782846853,
        48092.119152385276,
        61209.73700846825,
        189495.98940768326,
        233534.96450603567,
        186725.16725444607,
        140815.23524318123,
        60269.32989888545,
        51815.93369295262,
        49022.921055841725
    ]

    ADRIA.fecundity_scope!(
        fec_groups, fec_all, fec_params, C_cover_t, Matrix(total_loc_area')
    )

    @test any(fec_groups .> 1e8) ||
        "Fecundity is measured in m² and so should be a very large number"
    @test !any(fec_groups .< 0.0) || "Negative fecundity is not allowed"
end

@testset "Larval Production" begin
    tstep = 2
    a_adapt = fill(4.0, 36)
    n_adapt = 0.025
    dhw_scen = fill(4.0, 50)
    LPdhwcoeff = 0.4
    DHWmaxtot = 50.0
    LPDprm2 = 5.0
    n_groups = 6

    LPs = ADRIA.stressed_fecundity(
        tstep,
        a_adapt,
        n_adapt,
        dhw_scen[tstep - 1, :],
        LPdhwcoeff,
        DHWmaxtot,
        LPDprm2,
        n_groups
    )
    @test all(0.0 .<= LPs .< 1.0) || "Larval Production must be between 0 and 1"
end

@testset "Recruitment" begin
    n_locs = 334
    n_groups = 5

    total_loc_area = rand(n_locs) .* 1e6
    max_cover = rand(n_locs)
    avail_area = rand(n_locs)
    larval_pool = rand(n_groups, n_locs) .* 1e8

    recruits_per_m² = ADRIA.recruitment_rate(larval_pool, avail_area)
    abs_recruits = recruits_per_m² .* (avail_area .* max_cover .* total_loc_area)'

    @test any(abs_recruits .> 10^4) || "At least some recruitment values should be > 10,000"

    theoretical_max = ((avail_area .* max_cover .* total_loc_area)' * 51.8)
    for (i, rec) in enumerate(eachrow(abs_recruits))
        @test all(rec' .<= theoretical_max) ||
            "Species group $i exceeded maximum theoretical number of settlers"
    end
end

@testset "Mortality Scaling Tests" begin
    survival_probability::Matrix{Float64} = Matrix{Float64}(undef, 5, 7)
    expected_survival_p::Matrix{Float64} = Matrix{Float64}(undef, 5, 7)

    scaling_params::Vector{Float64} = [1.0, 1.0, 1.0, 1.0, 1.0]
    survival_probability .= 0.7
    expected_survival_p .= 0.85

    @test all(
        expected_survival_p .≈ ADRIA.apply_mortality_scaling(
            survival_probability, scaling_params
        )
    ) || "Given a survival probabilty of 0.7 and a scaling parameter of 1.0, expected 0.85"

    scaling_params .= -1.0
    survival_probability .= 0.7
    expected_survival_p .= 0.55

    @test all(
        expected_survival_p .≈ ADRIA.apply_mortality_scaling(
            survival_probability, scaling_params
        )
    ) || "Given a survival probabilty of 0.7 and a scaling parameter of -1.0, expected 0.65"

    scaling_params .= [1.0, -1.0, 0.5, -0.5, 0.0]
    expected_survival_p[1, :] .= 0.85
    expected_survival_p[2, :] .= 0.55
    expected_survival_p[3, :] .= 0.775
    expected_survival_p[4, :] .= 0.625
    expected_survival_p[5, :] .= 0.7

    @test all(
        expected_survival_p .≈ ADRIA.apply_mortality_scaling(
            survival_probability, scaling_params
        )
    ) || "Unexpected survival probabilities after scaling."
end
