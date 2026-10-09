using Plots
using ADRIA
using Statistics  # Required for the mean() function

#Double check I am using the dev version
ADRIA.hello_dev_version()

# --- Load Parameters directly from ADRIA package ---

# Define the names and number of functional groups in the order they appear in the functions
functional_group_names = [
    "Tabular Acropora",
    "Corymbose Acropora",
    "Pocillopora etc.",
    "Small Massives",
    "Large Massives",
]
n_sizes = 7
n_groups = length(functional_group_names)

# Get the standard deviations. This function returns one value per group,
# repeated for each size class. We can take every n_sizes-th value to get the unique std.
all_stds = ADRIA.dist_std(; n_sizes=n_sizes)
unique_stds = all_stds[1:n_sizes:end]

# Get the calibrated means. This returns a unique value for each size class.
all_means = ADRIA.dist_mean(; version=:calib, n_sizes=n_sizes)

# To plot a single representative curve per group, we calculate the average mean
# tolerance across all size classes for that group.
avg_means = zeros(n_groups)
for i in 1:n_groups
    start_idx = (i - 1) * n_sizes + 1
    end_idx = i * n_sizes
    avg_means[i] = mean(all_means[start_idx:end_idx])
end


# --- Constants from the ADRIA model ---
const DHW_range = 1.0:0.1:20.0
const LOWER_TRUNC_BOUND = 4.0  # Mortality is not calculated for DHW < 4.0
const HEAT_UB = 20.0           # An upper bound constant from the source code
const depth = 7.0              # Surface depth (m) for maximum bleaching effect

# Calculate the depth coefficient once
depth_coeff = ADRIA.depth_coefficient(depth)


# --- Generate and Display Plot ---
# Initialize the plot with titles and labels
plt = plot(
    title="Bleaching Mortality Curves from ADRIA Parameters",
    xlabel="Degree Heating Weeks (DHW)",
    ylabel="Proportional Mortality",
    xlims=(1, 20),
    ylims=(0, 1.05),
    legend=:topleft,
    size=(1000, 700) # Increased size for better legend readability
)

# Loop through each functional group and plot its mortality curve
for i in 1:n_groups
    name = functional_group_names[i]
    μ = avg_means[i]
    σ = unique_stds[i]

    # The upper bound for the truncated distribution is the mean + a constant
    upper_bound = μ + HEAT_UB

    # Calculate the mortality for the entire DHW range using broadcasting (.)
    affected_pop = ADRIA.truncated_normal_cdf.(
        DHW_range,
        μ,
        σ,
        LOWER_TRUNC_BOUND,
        upper_bound
    )
    mortality_curve = affected_pop .* depth_coeff

    # Add the curve to the plot with formatted labels
    plot!(
        plt,
        DHW_range,
        mortality_curve,
        label="$(name) (avg μ=$(round(μ, digits=2)), σ=$(round(σ, digits=2)))",
        linewidth=2.5
    )
end

# Add a vertical line to show the minimum DHW threshold
vline!(plt, [LOWER_TRUNC_BOUND], linestyle=:dash, color=:red, label="Mortality Threshold")

# Display the final plot
display(plt)