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
            3.6288, 3.39071, 3.51495, 3.05432, 3.46573, 3.54919, 2.99916,   # tabular Acropora
            4.04999, 3.90801, 3.96954, 3.55288, 3.87375, 3.79475, 3.77829,  # corymbose Acropora
            3.947, 4.03984, 4.06396, 4.01402, 3.76394, 3.81475, 4.00364,    # Pocillopora + non-Acropora corymbose
            5.70734, 5.64845, 5.39778, 5.72893, 5.36046, 5.51469, 5.01105,  # Small massives and encrusting
            6.63257, 6.38167, 6.18104, 5.90905, 6.01206, 6.11049, 6.0298    # Large massives
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
            0.283821  0.175423   0.190262   0.174482   0.166568   0.133093   0.137608;     # Tabular Acropora
            0.201831  0.118862   0.0851508  0.121045   0.100549   0.108431   0.118181;     # Corymbose Acropora
            0.197234  0.120599   0.0784736  0.0774731  0.0789747  0.0775685  0.077262;     # Corymbose non-Acropora
            0.219038  0.0732543  0.0403903  0.0241569  0.0132131  0.0172465  0.0252122;    # Small massives and encrusting
            0.259963  0.0753171  0.0439998  0.0252023  0.0143514  0.0161733  0.0289463     # Large massives
        ]
    end
    return error("Invalid param value `version`.")
end

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
            0.0266831  0.0456663   0.055042    0.0699477  0.0631659  0.0726078  0.0;   # Tabular Acropora
            0.0267012  0.0314565   0.0305842   0.0312019  0.0343453  0.0392649  0.0;   # Corymbose Acropora
            0.0192487  0.0183389   0.016142    0.0181973  0.0178266  0.0160711  0.0;   # Corymbose non-Acropora
            0.0127533  0.00815496  0.00798983  0.0100552  0.0143346  0.0153482  0.0;   # Small massives and encrusting
            0.0127479  0.00800844  0.00814327  0.0100825  0.0146049  0.0154561  0.0    # Large massives
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
        0.882941  1.08778  1.08675  0.922802  1.00579   0.778656  1.07741   0.818478  1.06056   1.08393   0.887496  0.962425;
        0.748283  1.25239  1.46479  0.737451  0.863211  0.715931  1.40074   0.884426  1.15045   1.13255   1.14274   1.06435;
        1.38788   1.4012   1.40572  0.702632  0.764635  0.82063   0.951614  0.728251  0.949598  1.03868   1.05603   1.12285;
        0.757405  1.16121  1.19449  0.739176  0.871387  1.23361   1.49762   1.0755    1.02776   0.73994   1.03886   1.21785;
        0.714787  1.26304  0.98837  1.44196   1.1571    0.714038  0.789416  1.48924   0.854136  0.709678  1.18619   1.20667
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
        -0.998468   0.951675    0.991806   -0.937867  -0.256171  -0.497685   0.5383    -0.956278   -0.095815   0.842615    0.0471974  0.920477;
        -0.058027   0.990021    0.986846   -0.46469   -0.97378   -0.536303   0.966195  -0.63045     0.405114  -0.241912   -0.198515   0.233474;
        -0.264616   0.973657   -0.0846454  -0.622472  -0.280372   0.659781   0.440735  -0.809034   -0.800206  -0.0415746   0.489644   0.876568;
        -0.814235  -0.395024    0.915013   -0.119188  -0.817617   0.244052   0.934592   0.978878   -0.616667  -0.59078     0.215869   0.336434;
        -0.965108  -0.0261854   0.474846    0.688      0.239802  -0.402938  -0.521073   0.0838699   0.367821   0.180813   -0.465029   0.935718
    ]
end
