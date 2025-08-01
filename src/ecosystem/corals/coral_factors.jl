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
            3.709682, 3.500189, 3.750201, 3.392543, 3.695615, 3.683136, 3.531146,   # tabular Acropora
            4.067916, 4.080612, 4.023438, 3.786189, 3.978531, 4.059314, 3.909487,   # corymbose Acropora
            4.247176, 4.383617, 4.374969, 4.171898, 4.250910, 4.407923, 4.458594,   # Pocillopora + non-Acropora corymbose
            5.903581, 6.085730, 5.790553, 5.809598, 5.969619, 6.012185, 5.889642,   # Small massives and encrusting
            7.149330, 6.967213, 6.859941, 7.017632, 6.892476, 7.07118, 6.739256     # Large massives
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
            0.284 0.177 0.193 0.183 0.157 0.149 0.151;  # Tabular Acropora
            0.204 0.119 0.096 0.112 0.100 0.107 0.111;  # Corymbose Acropora
            0.199 0.117 0.073 0.091 0.082 0.075 0.081;  # Corymbose non-Acropora
            0.220 0.073 0.046 0.024 0.013 0.015 0.026;  # Small massives and encrusting
            0.257 0.074 0.043 0.024 0.013 0.015 0.026   # Large massives
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
            0.026947 0.045658 0.055006 0.068457 0.060761 0.081604 0.0;    # Tabular Acropora
            0.026485 0.031391 0.030247 0.031083 0.032695 0.035496 0.0;    # Corymbose Acropora
            0.019048 0.018135 0.016260 0.018028 0.015115 0.018014 0.0;    # Corymbose non-Acropora
            0.012374 0.008204 0.008210 0.010310 0.013633 0.014312 0.0;    # Small massives and encrusting
            0.012564 0.008186 0.008138 0.010208 0.014300 0.015012 0.0     # Large massives
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
        0.82164   1.08925   1.09167   0.999908  1.00538   1.01946   0.937651  0.784938  0.956751  1.08086   0.991705  0.974071 0.974071;
        0.923883  1.49121   1.49719   1.16706   0.955033  1.02107   0.70478   0.793097  1.24596   0.806298  0.964745  0.91959  0.91959;
        1.17701   0.84639   1.19024   0.786443  1.35553   1.42635   0.834344  0.799714  0.760481  0.767771  1.07771   1.46848  1.46848;
        0.865704  1.12302   1.38518   0.711301  0.911399  1.42821   0.780267  1.49374   1.0102    0.70714   0.830893  1.31417  1.31417;
        0.819216  0.983586  0.884345  1.46174   0.818091  0.904022  0.715186  1.49572   1.32292   0.806857  1.16899   1.1661   1.1661
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
        -0.617747   0.979436   0.99745   -0.789161    0.391622   -0.149784   0.00874778  -0.914765   0.00325658  -0.195597  -0.237917   0.813194 0.813194;
        -0.378105   0.999229   0.999292  -0.967839   -0.795322    0.795927  -0.98076      0.161166  -0.0933056    0.365728   0.300143   0.087848 0.087848;
         0.508711  -0.535707   0.894669   0.0662139  -0.268865    0.622204  -0.622528    -0.555281   0.00126368  -0.396586  -0.169007   0.961894 0.961894;
        -0.537044  -0.139579   0.972233  -0.408186   -0.0325474   0.961359  -0.933316     0.535899   0.279357    -0.810235   0.0632827  0.471081 0.471081;
        -0.960488   0.0119783  0.24806    0.951627   -0.0412475  -0.788822  -0.507666     0.961015  -0.848863     0.705065  -0.207037   0.86103  0.86103
    ]
end
