using Test
using ADRIA
using ADRIA.decision:
    build_state,
    strategy_status,
    revisit_cadence_mask,
    filter_candidate_locations,
    PeriodicStrategy,
    ReactiveStrategy
using ADRIA: At

using DataStructures: CircularBuffer

if !@isdefined(ADRIA_DIR)
    const ADRIA_DIR = pkgdir(ADRIA)
    const TEST_DOMAIN_PATH = joinpath(ADRIA_DIR, "test", "data", "Test_domain")
end

if !@isdefined(ADRIA_DOM_45)
    const ADRIA_DOM_45 = ADRIA.load_domain(TEST_DOMAIN_PATH, 45)
end

# ─────────────────────────────────────────────────────────────────────────────
# revisit_cadence_mask
# ─────────────────────────────────────────────────────────────────────────────
@testset "revisit_cadence_mask" begin
    last_deployment = [10, 8, 0, 5]
    timestep = 10

    @testset "cadence <= 0 disables restriction" begin
        @test all(revisit_cadence_mask(last_deployment, timestep, 0)) ||
            "cadence == 0 should mark all locations eligible"
        @test all(revisit_cadence_mask(last_deployment, timestep, -1)) ||
            "negative cadence should mark all locations eligible"
    end

    @testset "eligibility threshold" begin
        # timestep - last_deployment: [0, 2, 10, 5]; last_deployment==0 -> always eligible
        mask = revisit_cadence_mask(last_deployment, timestep, 3)
        @test mask == BitVector([false, false, true, true]) ||
            "Expected only the never-deployed and long-idle locations to be eligible"
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# build_state ordering fix (Verification item 1)
# ─────────────────────────────────────────────────────────────────────────────
@testset "build_state ordering" begin
    dom = deepcopy(ADRIA_DOM_45)
    loc_ids = dom.loc_ids
    n = length(loc_ids)

    # Deliberately reorder a subset of loc_ids so the relative order differs from
    # domain.loc_ids's own order (interleaved, not a contiguous first/second half).
    reorder_idx = [5, 1, 4, 2, 3]
    target_locations = loc_ids[reorder_idx]

    # States keyed by domain.loc_ids's own order: value at position i identifies
    # location i unambiguously (its 1-based index in loc_ids).
    current_cover = Float64.(1:n)
    recent_cover_losses = CircularBuffer{Vector{Float64}}(1)
    push!(recent_cover_losses, Float64.(10 .* (1:n)))
    last_deployment = collect(100 .* (1:n))

    states = (
        current_cover=current_cover,
        recent_cover_losses=recent_cover_losses,
        last_deployment=last_deployment
    )

    # Expected: for each location in unique(target_locations) (in its own order),
    # the state values must correspond to THAT location's index in domain.loc_ids,
    # not the position within target_locations/domain.loc_ids traversal order.
    expected_domain_idx = [findfirst(==(loc), loc_ids) for loc in unique(target_locations)]

    @testset "ReactiveStrategy" begin
        strategy = ReactiveStrategy(target_locations, 0, 75, 0.2, 0.3, 0.0, 0, 0)
        state = build_state(dom, strategy, states)

        @test collect(state.current_cover) == Float64.(expected_domain_idx) ||
            "current_cover not aligned to target_locations' own order for ReactiveStrategy"
        @test collect(state.recent_cover_losses) == Float64.(10 .* expected_domain_idx) ||
            "recent_cover_losses not aligned to target_locations' own order for ReactiveStrategy"
        @test collect(state.last_deployment) == 100 .* expected_domain_idx ||
            "last_deployment not aligned to target_locations' own order for ReactiveStrategy"

        # Cross-check actual IDs: the location with the highest current_cover in the
        # returned state must be the target with the highest index in domain.loc_ids.
        max_pos = argmax(state.current_cover)
        expected_loc = unique(target_locations)[argmax(expected_domain_idx)]
        @test unique(target_locations)[max_pos] == expected_loc ||
            "Wrong location ID returned at the max-cover position"
    end

    @testset "PeriodicStrategy" begin
        strategy = PeriodicStrategy(target_locations, 0, 75, 1, 75, 0, 0)
        state = build_state(dom, strategy, states)

        @test collect(state.current_cover) == Float64.(expected_domain_idx) ||
            "current_cover not aligned to target_locations' own order for PeriodicStrategy"
        @test collect(state.recent_cover_losses) == Float64.(10 .* expected_domain_idx) ||
            "recent_cover_losses not aligned to target_locations' own order for PeriodicStrategy"
        @test collect(state.last_deployment) == 100 .* expected_domain_idx ||
            "last_deployment not aligned to target_locations' own order for PeriodicStrategy"
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# filter_candidate_locations — PeriodicStrategy (Verification item 2)
# ─────────────────────────────────────────────────────────────────────────────
@testset "filter_candidate_locations — PeriodicStrategy" begin
    targets = ["A", "B", "C", "D", "E"]

    @testset "cadence exclusion + backfill with stable tie-break" begin
        timestep = 10
        cadence = 3
        min_locations = 4

        # A, B, C deployed at t=10 (0 idle); D deployed at t=8 (idle=2); E never (idle=inf)
        last_deployment = [10, 10, 10, 8, 0]

        strategy = PeriodicStrategy(targets, 0, 75, 1, 75, cadence, min_locations)
        state = (last_deployment=last_deployment,)

        result = filter_candidate_locations(strategy, timestep, state)

        # Only E is eligible outright (count=1 < min_locations=4). Backfill needed: 3 more.
        # Excluded = [A,B,C,D], idle times = [0,0,0,2]. Sorted by idle desc, stable
        # tie-break preserves original order among ties: D (idle=2) first, then A,B,C
        # (tied at idle=0, in original order). Backfill picks the first 3: D, A, B.
        @test result == ["E", "D", "A", "B"] ||
            "Backfill did not select longest-idle locations with stable tie-break " *
              "(got $result)"
    end

    @testset "no backfill needed when cadence-eligible count already meets minimum" begin
        timestep = 10
        cadence = 3
        min_locations = 1

        last_deployment = [10, 10, 10, 8, 0]
        strategy = PeriodicStrategy(targets, 0, 75, 1, 75, cadence, min_locations)
        state = (last_deployment=last_deployment,)

        result = filter_candidate_locations(strategy, timestep, state)
        @test result == ["E"] ||
            "Expected only cadence-eligible locations when count meets min_locations " *
              "(got $result)"
    end

    @testset "revisit_cadence == 0 is a no-op (no filtering)" begin
        timestep = 10
        min_locations = 100  # deliberately impossible, to prove cadence short-circuits first
        last_deployment = [10, 10, 10, 10, 10]  # would all be excluded under any cadence > 0

        strategy = PeriodicStrategy(targets, 0, 75, 1, 75, 0, min_locations)
        state = (last_deployment=last_deployment,)

        result = filter_candidate_locations(strategy, timestep, state)
        @test result == targets ||
            "revisit_cadence == 0 should return all target locations unfiltered " *
              "(got $result)"
    end

    @testset "nothing state is also a no-op" begin
        timestep = 10
        strategy = PeriodicStrategy(targets, 0, 75, 1, 75, 5, 3)
        result = filter_candidate_locations(strategy, timestep, nothing)
        @test result == targets || "nothing state should return all target locations"
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# filter_candidate_locations — ReactiveStrategy (Verification item 2)
# ─────────────────────────────────────────────────────────────────────────────
@testset "filter_candidate_locations — ReactiveStrategy" begin
    targets = ["A", "B", "C", "D", "E"]

    @testset "cadence excludes locations that would otherwise pass reactive thresholds, no backfill" begin
        timestep = 10
        cadence = 3
        absolute_threshold = 0.2
        loss_threshold = 0.3
        min_cover_remaining = 0.0

        # All locations pass the absolute-threshold reactive trigger.
        current_cover = fill(0.1, 5)
        recent_cover_losses = zeros(5)

        # Same cadence setup as the periodic test: only E (never deployed) is eligible.
        last_deployment = [10, 10, 10, 8, 0]

        strategy = ReactiveStrategy(
            targets, 0, 75, absolute_threshold, loss_threshold,
            min_cover_remaining, 0, cadence
        )
        state = (
            current_cover=current_cover,
            recent_cover_losses=recent_cover_losses,
            last_deployment=last_deployment
        )

        result = filter_candidate_locations(strategy, timestep, state)

        # Reactive strategy does NOT backfill: even though only 1 of 5 locations is
        # eligible, the result must stay at that reduced count.
        @test result == ["E"] ||
            "Cadence should exclude on-cooldown locations without backfill " *
              "(got $result); NO_BACKFILL_BUG if a topped-up count is observed"
        @test length(result) == 1 ||
            "ReactiveStrategy must not backfill when cadence drops candidates below " *
              "any threshold (got $(length(result)) candidates)"
    end

    @testset "revisit_cadence == 0 behaves as pure reactive (no cadence restriction)" begin
        timestep = 10
        absolute_threshold = 0.2
        loss_threshold = 0.3
        min_cover_remaining = 0.0

        current_cover = fill(0.1, 5)
        recent_cover_losses = zeros(5)
        last_deployment = [10, 10, 10, 10, 10]  # would exclude everyone under cadence > 0

        strategy = ReactiveStrategy(
            targets, 0, 75, absolute_threshold, loss_threshold,
            min_cover_remaining, 0, 0
        )
        state = (
            current_cover=current_cover,
            recent_cover_losses=recent_cover_losses,
            last_deployment=last_deployment
        )

        result = filter_candidate_locations(strategy, timestep, state)
        @test result == targets ||
            "revisit_cadence == 0 should not restrict candidates by deployment recency " *
              "(got $result)"
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Multi-share end-to-end scenario test (Verification item 3)
# ─────────────────────────────────────────────────────────────────────────────
@testset "Multi-share revisit cadence — end-to-end" begin
    dom = deepcopy(ADRIA_DOM_45)
    loc_ids = dom.loc_ids
    n_locs = length(loc_ids)  # 10 for Test_domain

    half = n_locs ÷ 2
    # Share 1: reversed order relative to domain.loc_ids (out-of-canonical-order share)
    share_1_locs = reverse(loc_ids[1:half])
    # Share 2: natural order
    share_2_locs = loc_ids[(half + 1):end]

    weight_1 = 0.5
    weight_2 = 0.5

    # cadence=2 with min_iv_locations=3 over a 10-location aggregate pool (below)
    # guarantees at least 10-3=7 locations remain eligible every decision year
    # (>= min_iv_locations), so Config A never needs to invoke PeriodicStrategy's
    # backfill path — keeping this subtest a clean test of cadence exclusion alone.
    cadence = 2
    seed_year_start = 1
    seed_years = 10

    function _configure!(dom, min_iv_locations)
        ADRIA.set_seed_target_locations!(
            dom,
            [
                (weight=weight_1, target_locs=share_1_locs),
                (weight=weight_2, target_locs=share_2_locs)
            ]
        )
        ADRIA.fix_factor!(
            dom;
            N_seed_TA=500_000.0,
            N_seed_CA=500_000.0,
            N_seed_CNA=500_000.0,
            N_seed_SM=500_000.0,
            N_seed_LM=500_000.0,
            seed_year_start=Float64(seed_year_start),
            seed_years=Float64(seed_years),
            seed_deployment_freq=1.0,  # decision every year
            seed_revisit_cadence=Float64(cadence),
            seed_strategy=Float64(ADRIA.DECISION_STRATEGY[:periodic]),
            min_iv_locations=Float64(min_iv_locations)
        )
        return dom
    end

    num_samples = 2

    # ── Config A: min_iv_locations set low relative to the 10-location aggregate
    # pool, so that a clean 10/3-location round-robin can (in principle) satisfy
    # both the per-share minimum and the cadence window without needing backfill
    # on every decision year. Used to check cadence + correct-ID/ordering behaviour.
    dom_a = _configure!(deepcopy(dom), 3)
    scens_a = ADRIA.sample_guided(dom_a, num_samples)
    rs_a = ADRIA.run_scenarios(dom_a, scens_a, "45")
    # rs.ranks structure: (timesteps, locations, interventions, scenarios); rank > 0
    # (stored as UInt16) indicates the location was selected for deployment at that
    # timestep for that intervention type.
    seed_ranks_a = rs_a.ranks[intervention = At(:seed)]

    @testset "cadence is respected (no re-selection within the window)" begin
        # Config A is deliberately sized so PeriodicStrategy's min_locations backfill
        # (which intentionally re-offers on-cooldown locations when the eligible pool
        # would otherwise drop below min_iv_locations — see PeriodicStrategy's
        # docstring) never triggers here; this isolates cadence exclusion itself from
        # the backfill escape hatch, which is separately exercised in Config B below.
        for s = 1:num_samples
            for loc_idx = 1:n_locs
                deployed_ts = findall(
                    >(0), collect(seed_ranks_a[locations = loc_idx, scenarios = At(s)])
                )
                deployed_ts = deployed_ts[
                    seed_year_start .<= deployed_ts .<= (seed_year_start + seed_years - 1)
                ]
                isempty(deployed_ts) && continue
                gaps = diff(sort(deployed_ts))
                @test all(gaps .>= cadence) ||
                    "Location $(loc_ids[loc_idx]) (scenario $s) was re-selected within " *
                      "the $cadence-year cadence window (deployment timesteps: " *
                      "$deployed_ts)"
            end
        end
    end

    @testset "correct location IDs selected despite non-canonical share ordering" begin
        # Every location that was ever selected for seeding must belong to one of the
        # two target-location sets (cross-checked by actual ID, not just index/count).
        all_target_locs = Set(vcat(share_1_locs, share_2_locs))
        for s = 1:num_samples
            for loc_idx = 1:n_locs
                any_deployment = any(
                    >(0), collect(seed_ranks_a[locations = loc_idx, scenarios = At(s)])
                )
                if any_deployment
                    @test loc_ids[loc_idx] in all_target_locs ||
                        "Location $(loc_ids[loc_idx]) was seeded but is not part of any " *
                          "configured target-location share"
                end
            end
        end

        # Reversed share_1 specifically: confirm at least one of its members (which sit
        # at non-canonical positions in the aggregate target list) was actually seeded,
        # proving the reversed ordering did not silently drop/misalign them.
        share_1_idx = findall(in(share_1_locs), loc_ids)
        share_1_deployed = any(
            any(>(0), collect(seed_ranks_a[locations = li, scenarios = At(s)]))
            for li in share_1_idx, s = 1:num_samples
        )
        @test share_1_deployed ||
            "No location in the reversed-order share (share_1_locs) was ever seeded — " *
              "possible ordering/alignment regression"
    end

    # ── Config B: min_iv_locations set high relative to the 10-location aggregate
    # pool (unsustainable under a strict 3-year cadence round-robin), to force
    # PeriodicStrategy's backfill path on essentially every decision year.
    min_iv_locations_b = 8
    dom_b = _configure!(deepcopy(dom), min_iv_locations_b)
    scens_b = ADRIA.sample_guided(dom_b, num_samples)
    rs_b = ADRIA.run_scenarios(dom_b, scens_b, "45")
    seed_ranks_b = rs_b.ranks[intervention = At(:seed)]

    @testset "Periodic backfill fills the aggregate pool when needed" begin
        found_backfill_evidence = false
        for s = 1:num_samples
            for tstep = (seed_year_start + 1):(seed_year_start + seed_years - 1)
                deployed_this_year = findall(
                    >(0), collect(seed_ranks_b[timesteps = tstep, scenarios = At(s)])
                )
                if length(deployed_this_year) >= min_iv_locations_b
                    found_backfill_evidence = true
                end
            end
        end
        @test found_backfill_evidence ||
            "Expected at least one decision year where the aggregate selected-location " *
              "count reached min_iv_locations ($min_iv_locations_b), evidencing " *
              "PeriodicStrategy backfill"
    end
end
