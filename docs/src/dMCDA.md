# Dynamic Multi-Criteria Decision Analysis

## Site Selection in ADRIA
Site selection in ADRIA is performed by iteratively updating how each site/reef in the domain performs against a range of criteria. Currently, ADRIA allows the following criteria to be used in site selection:
1.	Connectivity, calculated as centrality*(site coral cover)*(site area),
2. 	Wave damage probability,
3.	Heat stress probability,
4.	Priority predecessors, or sites with the strongest connectivity to ‘priority sites’ (defined by the user in the prioritysites parameter),
5.	Priority zones, such as GBRMPA management zones (defined by the user in the priorityzones parameter),
6.	Coral cover relative to the site carrying capacity.

In the case of site selection for seeding, the coral cover criterium is defined as the area not occupied by coral relative to the site carrying capacity, thus giving a measure of the real estate available to baby corals. For shading, the coral cover criterium is defined as the area already covered with coral relative to the carrying capacity, giving a measure of the coral population available to benefit from shading. All other criteria are the same for seeding and shading, such that sites selected for these two interventions tend to be similar, unless coral cover is weighted very highly relative to the other criteria. The priority predecessors and priority zones criteria allow the user to incorporate spatial objectives, such as through GBRMPA management zone layers or data largers from other spatial management software such as MARXAN.

These criteria are used within methods based on multi-criteria decision analysis (MCDA) to rank sites from most to least suitable for intervention at each intervention time step when ADRIA is run. The criteria weights can be set in ADRIA through the input scenario parameters csv, or inline, as in the following example.

```julia
@info "Loading data package"
here = @__DIR__
ex_domain = ADRIA.load_domain(joinpath(here, "Example_domain"))

@info "Creating 5 scenarios based on parameter bounds using the Sobol' method"
scens = ADRIA.sample(ex_domain, 5)

@info "Adjusting heat_stress and wave_stress weightings"
scens.wave_stress .= [0.7,0.6,0.1,0.4]
scens.heat_stress .= [0.0,0.2,0.1,0.7]

```
The criteria weights and parameters which the user sets in ADRIA are described in the table below.

|Parameter              |Description                                                                              |
|-----------------------|-----------------------------------------------------------------------------------------|
|heat_stress            |Weight for heat stress when selecting seeding or shading sites.                          |
|wave_stress            |Weight for wave stress when selecting seeding or shading sites.                          |
|shade_connectivity     |Weight for larval connectivity when selecting shading sites.                             | 
|in_seed_connectivity   |Weight for in-degree larval connectivity when selecting seeding sites.                   |
|out_seed_connectivity  |Weight for out-degree larval connectivity when selecting seeding sites.                  |
|coral_cover_high       |Weight of coral cover area covered when selecting shading sites.                         |
|coral_cover_low        |Weight of coral cover area available for coral gorwth when selecting seeding sites.      |
|seed_priority          |Weight of being a priority predecessor site to a priority site when seeding.             |
|shade_priority         |Weight of being a priority predecessor site to a priority site when shading.             |
|coral_cover_tol        |% of area of corals to be seeded allowed on a site when selecting sites for seeding.     |
|deployed_coral_risk_tol|risk tolerance for heat and wave risk when seeding or shading corals.                    |
|depth_min              |Minimum depth of sites when selecting for seeding or shading corals.                     |
|depth_offset           |depth_min+depth_offset is maximum depth of sites when seeding or shading corals.         |

## MCDA methods
Multi-criteria decision analysis encompasses a series of algorithms for conducting decision making in a way that is formalised, structured and transparent. It evaluates the decision objective according to numerical measures of decision criteria assembled by the decision maker. Constructing explicit criteria over which to evaluate alternatives allows a clear evaluation of benefits, negatives and trade-offs in coming to a decision solution.

Commonly, comparison of criteria in MCDA is acheived through the construction of a “values matrix” or “decision matrix”, which contains measures of value for each alternative against an array of criteria. The measure of value, for example, could be a relative ranking of worth, probability of success, or risk of failure. Formally, these values form an $N$ by $M$ matrix where entry $x_{i,j}$ is the value of the ith alternative as evaluated against the jth criterium. The decision matrix is generally normalised, so that

$$\sum{i=1}^N\sqrt{x_{i,j}^2}=1 \forall j.$$

The user also weights each criterium with a weighting representing it's relative importance in the decison problem (which could be decided through stakeholder engagement or expert opinion for example). These weightings $w_j, j =1,..,M$, are also normalised so that
$$\sum{j=1}^M w_j=1.$$

The final decision matrix used in decision algorithms in ADRIA has the form,
$$X_{i,j} = w_j x_{i,j}$$

The particular MCDA technique chosen prescribes a strategy for aggregating the values measures in the decision matrix and forming rankings of alternatives. Perhaps the simplest way to aggregate criteria values, referred to here as ‘order ranking’, is to simply add up the criteria values for each alternative. This additive score can then be ordered from highest to lowest, giving an order of preference for the alternatives considered in the decision. Simply adding criteria values, however, can mask trade-offs between different criteria. Many ‘compensatory aggregation’ algorithms exist to combat this issue, including TOPSIS and VIKOR. 


TOPSIS, or Technique for Order of Preference by Similarity to Ideal Solution, ranks alternatives by comparing them to the ‘Positive Ideal Solution’ (PIS) and ‘Negative Ideal Solution’ (NIS). The PIS is the value of the highest valued alternative for a criterium we want to maximise (or lowest if we want to minimise). The NIS is the value of the lowest valued alternative for a criterium we want to maximise (or highest if we want to minimise). TOPSIS ranks alternatives according to a ratio calculated from the geometric distance of the alternative to the PIS and NIS in each of the criteria. This allows alternatives which perform poorly under some criteria to balance their ranking with a high performance in other criteria, providing a measure of trade-offs between criteria. If it is desired to avoid a very bad performance in all criteria, however, this method can mask alternatives which perform very badly in only a few criteria but very well in others.

VIKOR is another very popular MCDA algorithm which attempts to combat the hidden extremes issue of TOPSIS in its formulation. It does this using two distance formulations; the Manhattan distance (or L1 distance) and Chebyshev distance. The Manhattan distance, S, for each alternative is the sum over all criteria of the distance between each criteria’s value and its PIS, giving a measure of the overall utility of each alternative when considering all criteria. I.e. this is the general performance of an alternative when considering all criteria, with possible hidden trade-offs. The Chebyshev distance, R, is the maximum distance between the PIS and the value of a criteria for each alternative, giving a measure of the maximum deviance from the best solution for that alternative. I.e. this gives a measure of the worst hidden trade-off in the Manhattan Distance (which is not considered in TOPSIS). The final ranking is calculated from a linear combination of R and S where the weightings of R and S are determined by the balance of ‘group utility’ (measured by S) and ‘individual regret’ (measured by R) a decision maker wants to incorporate in their decision solution. Additional steps to the algorithm give guidelines for choosing between a single best alternative and n (n>1) best alternatives.

The decision strategies for each of these algorithms can be summarised as follows:
* Order ranking: when overall performance matters and trade-offs between criteria need not be considered (also least computationally expensive).
* TOPSIS: when trade-offs between criteria should be considered, but it is not necessary to avoid hidden value extremes.
* VIKOR: When trade-offs between criteria should be considered and the decision-maker wants to weight against hidden very poorly performing criteria.

The algorithm used in an ADRIA run is set using the 'guided' parameter, as in the example below.

```julia
@info "Loading data package"
here = @__DIR__
ex_domain = ADRIA.load_domain(joinpath(here, "Example_domain"))

@info "Creating 5 scenarios based on parameter bounds using the Sobol' method"
scens = ADRIA.sample(ex_domain, 5)

@info "Setting the guided parameter for each run"
scens.guided[1] = -1.0 # counterfactual run
scens.guided[2] = 0.0 # unguided selection (randomised)
scens.guided[3] = 1.0 # order ranking
scens.guided[4] = 2.0 # TOPSIS
scens.guided[5] = 3.0 # VIKOR
```

## Using site selection separately
ADRIA's site selection algorithms can also be used outside of the ecological model, using the site_selection function. The user provides data layers representing the selection criteria and weights to produce a set of site ranks for seeding and/or shading. An example of the function's usage is shown below.

```julia
@info "Loading data package"
here = @__DIR__
ex_domain = ADRIA.load_domain(joinpath(here, "Example_domain"))

alg_ind = 2.0 # use TOPSIS

# define weights for criteria
criteria = DataFrame(wave_stress = 0.2,
                    heat_stress = 0.4,
                    shade_connectivity = 0.5,
                    out_seed_connectivity = 0.4,
                    in_seed_connectivity = 0.8,
                    coral_cover_high = 0.8,
                    coral_cover_low = 0.6,
                    seed_priority = 0.5,
                    shade_priority = 0.6,
                    coral_cover_tol = 0.2,
                    deployed_coral_risk_tol = 1.0,
                    depth_min = 5,
                    depth_offset = 10)

area_to_seed = 5/10^6 # area to seed in km^2
ts = 5 # year/time step you want to make a decision at
n_reps = 50 # number of dhw/wave replicates

ranks = site_selection(ex_domain, criteria, area_to_seed, ts, n_reps, alg_ind)
```