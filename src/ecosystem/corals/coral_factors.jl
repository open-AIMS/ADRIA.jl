# Base values for coral factors

"""
    dist_std(; n_sizes=7)::Vector{Float64}

Natural adaptation / heritability values here informed by [1] and
(unpublished) data from [2].

# References
1. Bairos-Novak, K.R., Hoogenboom, M.O., van Oppen, M.J.H., Connolly, S.R., 2021.
   Coral adaptation to climate change: Meta-analysis reveals high heritability across
     multiple traits.
   Global Change Biology 27, 5694-5710.
   https://doi.org/10.1111/gcb.15829

2. Hughes, T. P., Kerry, J. T., Baird, A. H., Connolly, S. R., Dietzel, A., Eakin, C. M.,
     Heron, S. F., Hoey, A. S., Hoogenboom, M. O., Liu, G., McWilliam, M. J., Pears, R. J.,
     Pratchett, M. S., Skirving, W. J., Stella, J. S., & Torda, G. (2018).
   Global warming transforms coral reef assemblages.
   Nature, 556(7702), 492-496.
   https://doi.org/10.1038/s41586-018-0041-2
"""
function dist_std(; n_sizes=7)::Vector{Float64}
    return repeat(
        Float64[
            # 2.590016677,  # arborescent Acropora
            2.904433676,  # tabular Acropora
            3.159922076,  # corymbose Acropora
            3.474118416,  # Pocillopora + non-Acropora corymbose
            4.773419097,  # Small massives and encrusting
            5.538122776   # Large massives
        ]; inner=n_sizes)
end

"""
    dist_mean(; version=:calib, n_sizes=7)::Vector{Float64}

If `version==:legacy` returns natural adaptation / heritability values informed by
[1] and (unpublished) data from [2].
If `version==:calib` returns values resulting from the model calibration.

# References
1. Bairos-Novak, K.R., Hoogenboom, M.O., van Oppen, M.J.H., Connolly, S.R., 2021.
   Coral adaptation to climate change: Meta-analysis reveals high heritability across
     multiple traits.
   Global Change Biology 27, 5694-5710.
   https://doi.org/10.1111/gcb.15829

2. Hughes, T. P., Kerry, J. T., Baird, A. H., Connolly, S. R., Dietzel, A., Eakin, C. M.,
     Heron, S. F., Hoey, A. S., Hoogenboom, M. O., Liu, G., McWilliam, M. J., Pears, R. J.,
     Pratchett, M. S., Skirving, W. J., Stella, J. S., & Torda, G. (2018).
   Global warming transforms coral reef assemblages.
   Nature, 556(7702), 492-496.
   https://doi.org/10.1038/s41586-018-0041-2
"""
function dist_mean(; version=:calib, n_sizes=7)::Vector{Float64}
    if version == :legacy
        return repeat(
            Float64[
                # 3.345484656,  # arborescent Acropora
                3.751612251,  # tabular Acropora
                4.081622683,  # corymbose Acropora
                4.487465256,  # Pocillopora + non-Acropora corymbose
                6.165751937,  # Small massives and encrusting
                7.153507902   # Large massives
            ]; inner=n_sizes)
    elseif version == :calib
        return [
            3.70071, 3.44978, 3.68859, 3.23064, 3.54517, 3.61326, 3.20368,   # tabular Acropora
            4.06432, 3.93149, 4.00685, 3.75739, 3.93149, 3.87897, 3.88147,   # corymbose Acropora
            4.11692, 4.34996, 4.19347, 4.07487, 3.826, 4.10226, 4.29725,   # Pocillopora + non-Acropora corymbose
            5.8865, 5.96419, 5.48963, 5.7387, 5.55563, 5.77022, 5.5164,   # Small massives and encrusting
            6.70074, 6.87894, 6.55428, 6.39785, 6.47497, 6.54627, 6.07672    # Large massives
        ]
    end
    return error("Invalid param value `version`.")
end

"""
    mortality_base_rate(; version=:calib)
If `version==:legacy` returns values informed by EcoRRAP (unpublished) data.
If `version==:calib` returns values resulting from the model calibration.
"""
function mortality_base_rate(; version=:calib)
    return if version == :legacy
        return 1.0 .- [
            0.60 0.76 0.81 0.76 0.85 0.86 0.86;     # Tabular Acropora
            0.60 0.76 0.77 0.87 0.83 0.90 0.90;     # Corymbose Acropora
            0.52 0.77 0.77 0.87 0.89 0.98 0.98;     # Corymbose non-Acropora
            0.72 0.87 0.77 0.98 0.99 0.99 0.99;     # Small massives and encrusting
            0.58 0.87 0.78 0.98 0.98 0.98 0.98      # Large massives
        ]
    elseif version == :calib
        return [
            0.282038  0.175418   0.190753   0.184469   0.16062    0.144075   0.141913;  # Tabular Acropora
            0.203026  0.120358   0.0826026  0.113109   0.10184    0.102131   0.110558;  # Corymbose Acropora
            0.200871  0.11627    0.0709936  0.0872101  0.0826969  0.0775613  0.0770466;  # Corymbose non-Acropora
            0.214515  0.0729543  0.043159   0.02497    0.012377   0.0145106  0.0252974;  # Small massives and encrusting
            0.289532  0.0733323  0.0415098  0.0253114  0.0127778  0.0154392  0.0284772  # Large massives
        ]
    end
    return error("Invalid param value `version`.")
end
# survival_rate::Matrix{Float64} = [
#     0.859017851 0.858528906 0.857044217 0.856477498 0.856104353 0.855852241 0.855852241;    # Tabular Acropora
#     0.865006527 0.87915437 0.892044073 0.905304164 0.915373252 0.925707536 0.925707536;     # Corymbose Acropora
#     0.953069031 0.959152694 0.964460394 0.968306361 0.972598906 0.97621179 0.97621179;     # Corymbose non-Acropora
#     0.869976692 0.938029324 0.977889252 0.987199004 0.99207702 0.996931548 0.996931548;     # Small massives and encrusting
#     0.9782479 0.979496637 0.980850254 0.982178103 0.983568572 0.984667677 0.984667677       # Large massives
# ]

"""
    linear_extensions()

Linear extensions. The values are converted from `cm` to the desired unit.
The default unit is `m`.
If `version==:legacy` returns values informed by EcoRRAP (unpublished) data.
If `version==:calib` returns values resulting from the model calibration.
"""
function linear_extensions(; unit=:m, version=:calib)::Matrix{Float64}
    if version == :legacy
        return [
            0.609456 1.071840 2.551490 5.079880 9.450910 16.8505 0.0;       # Tabular Acropora
            0.768556 1.220850 1.864470 2.822970 3.529380 3.00422 0.0;       # Corymbose Acropora
            0.190455 0.343747 0.615467 0.974770 1.700790 2.91729 0.0;       # Corymbose non-Acropora
            0.318034 0.473850 0.683729 0.710587 0.581085 0.581085 0.0;      # Small massives and encrusting
            0.122478 0.217702 0.382098 0.718781 1.241720 2.08546 0.0        # Large massives
        ] .* linear_scale(:cm, unit)
    elseif version == :calib
        return [
            0.0268946  0.0456585   0.0539942   0.0692523  0.072782   0.0758184  0.0;    # Tabular Acropora
            0.0266253  0.0314461   0.030562    0.0309425  0.0342009  0.0384155  0.0;    # Corymbose Acropora
            0.0174539  0.0184296   0.014432    0.0173405  0.017545   0.0182012  0.0;    # Corymbose non-Acropora
            0.0127291  0.00804491  0.0075679   0.0102268  0.0145612  0.0150638  0.0;    # Small massives and encrusting
            0.0127959  0.00810474  0.00818616  0.0100728  0.0143641  0.0151507  0.0    # Large massives
        ] .* linear_scale(:m, unit)
    end
    return error("Invalid value $version for param `version`.")
end

"""
    bin_edges()

Helper function defining coral colony diameter bin edges. The values are converted from `cm`
to the desired unit. The default unit is `m`.
"""
function bin_edges(; unit=:m)
    return Matrix(
        [
            2.5 7.5 12.5 25.0 50.0 80.0 120.0 160.0;
            2.5 7.5 12.5 20.0 30.0 60.0 100.0 150.0;
            2.5 7.5 12.5 20.0 30.0 40.0 50.0 60.0;
            2.5 5.0 7.5 10.0 20.0 40.0 50.0 100.0;
            2.5 5.0 7.5 10.0 20.0 40.0 50.0 100.0
        ]
    ) .* linear_scale(:cm, unit)
end
#   Matrix([
#         0.0 1.0 2.0 6.0 15.0 36.0 89.0 90.0;
#         0.0 1.0 2.0 4.0  9.0 18.0 38.0 39.0;
#         0.0 1.0 2.0 4.0  7.0 14.0 27.0 28.0;
#         0.0 1.0 2.0 5.0  8.0 12.0 26.0 27.0;
#         0.0 1.0 2.0 4.0  9.0 19.0 40.0 41.0
#   ]

"""
    planar_area_params()

Colony planar area parameters (see Fig 2B in [1])
First column is `b`, second column is `a`
log(S) = b + a * log(x)

# References
1. Aston Eoghan A., Duce Stephanie, Hoey Andrew S., Ferrari Renata (2022).
    A Protocol for Extracting Structural Metrics From 3D Reconstructions of Corals.
    Frontiers in Marine Science, 9.
    https://doi.org/10.3389/fmars.2022.854395
"""
function planar_area_params()
    return Array{Float64,2}([
        # -8.97 3.14    # Abhorescent Acropora (using branching porites parameters as similar method of growing ever expanding colonies).
        -8.95 2.80      # Tabular Acropora
        -9.13 2.94      # Corymbose Acropora
        -8.90 2.94      # Corymbose non-Acropora (using branching pocillopora values from fig2B)
        -8.87 2.30      # Small massives
        -8.87 2.30      # Large massives
    ])
end

"""
    linear_extension_group_scale_factors()::Matrix{Float64}

Matrix with dimensions (functional_groups ⋅ cb_calib_groups) where each element represents
the scale factor to be applied to the `linear_extensions`` of all size classes of that
`functional_group` and `cb_calib_group`.
"""
function linear_extension_group_scale_factors()::Matrix{Float64}
    return [
        0.806835  0.876033  1.36114   0.93817   0.846152  1.06008   1.31028   0.818468  1.154     0.782852  1.09145   1.48357;
        0.769908  1.31033   1.39064   0.873504  0.704694  0.721531  0.722839  1.47735   1.05525   0.756351  0.859923  0.726549;
        1.21256   0.92972   0.815334  1.24875   1.00857   0.710315  1.07801   1.4278    1.40837   1.49243   0.982459  0.791796;
        0.863153  0.77229   1.19286   1.25985   1.0634    1.38277   0.978835  0.929562  0.734949  1.07017   0.785061  0.878392;
        0.818373  0.870761  0.900768  1.10552   1.17317   0.933766  0.9214    1.02332   1.08906   1.41145   1.27386   1.29516
    ]
end

"""
    mb_rate_group_scale_factors()::Matrix{Float64}

Matrix with dimensions (functional_groups ⋅ cb_calib_groups) where each element represents
the scale factor to be applied to the `mb_rates` of all size classes of that
`functional_group` and `cb_calib_group`.
"""
function mb_rate_group_scale_factors()::Matrix{Float64}
    return [
        -0.250654   0.157892   0.511453   -0.520191   -0.91249     0.959655   0.572635   0.704358  0.19608   -0.510319   0.996532    0.99122;
        -0.722423   0.648037   0.181218   -0.869958    0.021476   -0.326979  -0.558814  -0.472231  0.25488   -0.473727  -0.0745634  -0.831083;
        -0.453511  -0.91636   -0.763698   -0.0440058   0.424088   -0.60318    0.441798   0.999925  0.651033   0.973715   0.621306   -0.92931;
        -0.178026  -0.700117   0.772009    0.961315   -0.0707886  -0.159285  -0.277358  -0.266131  0.863607   0.660176  -0.552144   -0.705929;
        -0.406945  -0.466123   0.0755071  -0.212037    0.214754   -0.250614   0.359056   0.165012  0.105881   0.982599   0.387922    0.925376
    ]
end
