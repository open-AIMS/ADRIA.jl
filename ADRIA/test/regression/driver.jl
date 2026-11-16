# Pinned baseline driver: runs a small fixed-seed ensemble on an ADRIA checkout.
#
#   julia --project=<checkout>/ADRIA -t4 driver.jl <checkout>/ADRIA <outdir>
#
# Writes one result store per case under <outdir>/<case>/ (ADRIA_OUTPUT_DIR is set per
# case). Nothing is written inside the checkout. Cases:
#   mixed    guided + unguided scenarios with interventions, COTS off and on (plus juvenile
#            predation on for some rows): exercises outcomes, ranks, intervention logs and
#            the COTS logs.
#   cots_off the same guided scenarios, all with COTS off (no COTS logs in the store).
#
# Per-stage variants: if a stage renames accessors, copy this file to e.g. driver_stage1.jl
# and point ADRIA_REGRESSION_DRIVER at it (see README.md). Keep the ensemble identical.
using Random
using ADRIA

const N_GUIDED = 8
const N_UNGUIDED = 4
const SEED = 20261116

function build_cases(dom)
    Random.seed!(SEED)
    guided = ADRIA.sample_guided(dom, N_GUIDED)
    unguided = ADRIA.sample_unguided(dom, N_UNGUIDED)

    off = copy(guided)
    off.cots_on .= 0

    on = copy(guided)
    on.cots_on .= 1
    on.cots_juvenile_predation[1:2:end] .= 1

    unguided_on = copy(unguided)
    unguided_on.cots_on .= 1

    mixed = vcat(off, on, unguided, unguided_on)
    return Dict("mixed" => mixed, "cots_off" => off)
end

function main(checkout::String, outdir::String)
    mkpath(outdir)
    # Apply the default configuration now: otherwise the first run_scenarios call does it
    # and overrides ADRIA_OUTPUT_DIR (the cwd must not hold a config.toml; the runner uses
    # an empty scratch dir)
    ADRIA.setup()
    dom = ADRIA.load_domain(joinpath(checkout, "test", "data", "Test_domain"), "45")
    cases = build_cases(dom)
    for name in sort!(collect(keys(cases)))
        case_dir = joinpath(outdir, name)
        mkpath(case_dir)
        withenv("ADRIA_OUTPUT_DIR" => case_dir) do
            t = @elapsed ADRIA.run_scenarios(
                dom, cases[name], "45"; show_progress=false
            )
            println(
                "case $(name): $(size(cases[name], 1)) scenarios, $(round(t; digits=1)) s"
            )
        end
    end
    return nothing
end

main(abspath(ARGS[1]), abspath(ARGS[2]))
