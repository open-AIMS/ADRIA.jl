# Dynamic Multi-Criteria Decision Analysis


## MCDA methods
Multi-criteria decision analysis encompasses a series of methods for conducting decision making in a way that 
is formalised, structured and transparent. It evaluates the decision objective according to numerical measures 
of decision criteria assembled by the decision maker. Constructing explicit criteria over which to evaluate 
alternatives allows a transperant evaluation of benefits, negatives and trade-offs in coming to a decision 
solution.

Commonly, comparison of criteria in MCDA is acheived through the construction of a ''values matrix'' or ''decision 
matrix'', which contains measures of value for each alternative against an array of criteria. 
The measure of value, for example, could be a relative ranking of worth, probability of success, or risk 
of failure. Formally, these values form an $N$ by $M$ matrix where entry $x_{i,j}$ is the value of the 
$i$th alternative as evaluated against the $j$th criterium. The decision matrix is generally normalised, 
so that

$$\sum_{i=1}^N\sqrt{x_{i,j}^2}=1 \forall j.$$

The user also weights each criterium with a weighting representing it's relative importance in the decison 
problem (which could, for example, be decided through stakeholder engagement or expert opinion). These 
weightings $w_j, j =1,...,M$, are also normalised so that
$$\sum_{j=1}^M w_j=1.$$

The final decision matrix used in the criteria ranking algorithms in ADRIA has the form,
$$X_{i,j} = w_j x_{i,j}$$

The particular MCDA technique chosen prescribes a strategy for aggregating criteria in the decision matrix and 
forming rankings of alternatives. ADRIA allows the user to select from 3 main aggregation methods, which are set 
through the ```guided``` parameter, which can be set in the criteria DataFrame. This is a regular DataFrame, with 
values which can be set as in the example below.

```julia
# Loading data package.
dom = ADRIA.load_domain("path to a domain data package")

# Creating 5 scenarios based on parameter bounds using the Sobol' method
scens = ADRIA.sample(dom, 5)

# Setting the guided parameter for each run
scens.guided[1] = -1.0 # counterfactual run
scens.guided[2] = 0.0 # unguided selection (randomised)
scens.guided[3] = 1.0 # order ranking
scens.guided[4] = 2.0 # TOPSIS
scens.guided[5] = 3.0 # VIKOR
```

## Aggregation methods
# Order Ranking
The simplest method available in ADRIA is referred to as ‘order ranking’, which simply sums the criteria values 
for each alternative. This additive score can then be ordered from highest to lowest, giving an order of preference 
for the alternatives considered in the decision. In this case, the rank $r_i$ of site $i$ is,

$$r_i = \sum_{j=1}^M X_{i,j}.$$

Simply adding criteria values, however, can mask trade-offs between different criteria. Many ‘compensatory 
aggregation' algorithms exist to combat this issue, including TOPSIS and VIKOR. 

# TOPSIS
TOPSIS, or Technique for Order of Preference by Similarity to Ideal Solution, ranks alternatives by comparing 
them to the 'Positive Ideal Solution' ($PIS$) and 'Negative Ideal Solution' ($NIS$). The $PIS$ is the value of 
the highest valued alternative for a criterium we want to maximise (or lowest if we want to minimise). The $NIS$ 
is the value of the lowest valued alternative for a criterium we want to maximise (or highest if we want to 
minimise). TOPSIS ranks alternatives according to a ratio calculated from the geometric distance of the alternative 
to the $PIS$ and $NIS$ in each of the criteria. 

The $PIS$ and $NIS$ for each criteria $j$ can be defined as,

$$PIS_j = \textbf{max}_{i} X_{i,j},$$
$$NIS_j = \textbf{min}_{i} X_{i,j}.$$

The final aggregate score in TOPSIS uses the $L^2$ distance between alternative $i$ and the NIS and PIS:

$$d_{NIS} = \sqrt{\sum_{j=1}^M(X_{i,j}-NIS_j)^2},$$
$$d_{PIS} = \sqrt{\sum_{j=1}^M(X_{i,j}-PIS_j)^2}.$$

The TOPSIS aggregate score used to rank alternatives is then,
$$C = d_{NIS}/(d_{NIS}+d_{PIS}).$$

The quantity $C$ measures the closeness of a particular alternative to the worst and best alternatives, such 
that the higher $C$, the closer to the best alternative and further from the worst. This method allows 
alternatives which perform poorly under some criteria to balance their ranking with a high performance in other 
criteria, providing a measure of trade-offs between criteria.

# VIKOR
VIKOR is another very popular MCDA algorithm, which incorporates two distance formulations; the Manhattan distance 
(or $L^1$ distance) and Chebyshev distance. The Manhattan distance, $S$, for each alternative is the sum over all 
criteria of the distance between each criteria’s value and its $PIS$, giving a measure of the overall utility of 
each alternative when considering all criteria. 

$$S_i = \sum_{j=1}^M \frac{PIS_j-X_{i,j}}{PIS_j-NIS_j}$$

I.e. this is the general performance of an alternative when considering all criteria, with possible trade-offs. 
The Chebyshev distance, $R$, is the maximum distance between the PIS and the value of a criteria for each 
alternative, giving a measure of the maximum deviance from the best solution for that alternative. 

$$R_i = \textbf{max}_{j} \frac{PIS_j-X_{i,j}}{PIS_j-NIS_j}$$

I.e. this gives a measure of the worst hidden trade-off in the Manhattan Distance (which is not considered in 
TOPSIS). The final ranking $Q$ is calculated from a linear combination of $R$ and $S$ where the weighting of $R$ 
and $S$ are determined by the balance of ‘group utility’ (measured by $S$) and ‘individual regret’ (measured by 
$R$) a decision maker wants to incorporate in their decision solution.

$$S^- = \textbf{min}_{i}(S_i),S^+ = \textbf{max}_{i}(S_i)$$
$$R^- = \textbf{min}_{i}(R_i),R^+ = \textbf{max}_{i}(R_i)$$
$$Q_i = \delta \frac{S_i-S^-}{S^+-S^-} + (1-\delta)\frac{R_i-R^-}{R^+-R^-}$$

Here, $\delta$ is a weighting which balances strategy between maximising all criteria and regret due to the worst 
performing criteria. Additional steps to the algorithm give guidelines for choosing between a single best 
alternative and n (n>1) best alternatives.

## Site Selection in ADRIA
At each time step where an intervention is performed at a site, the sites for intervention can be 
selected either randomly or by using a site selection algorithm derived from Multi-criteria 
Decision Analysis methods. When an algorithm is used, site selection is performed by iteratively 
updating how each site/reef in the domain performs against a range of criteria. Currently, 
ADRIA allows the following criteria to be used in site selection:

* Connectivity, calculated as centrality*(site coral cover)*(site area).
* Wave damage probability.
* Heat stress probability.
* Priority predecessors, or sites with the strongest connectivity to "priority sites" (defined by the user in the `priority_sites` parameter).
* Priority zones, such as GBRMPA management zones (defined by the `priority_zones` parameter in order of preference).
* Coral cover relative to the site carrying capacity.
* Site depth, for which median site depth is used.

For seeding, the coral cover criterium is defined as the area not occupied by coral relative to the site carrying 
capacity. For shading, the inverse is used for the coral cover criterium. The priority predecessors and 
priority zones criteria allow the user to incorporate spatial objectives, such as through GBRMPA management 
zone layers or data layers from other spatial management software such as MARXAN.

These criteria are used within to rank sites from most to least suitable for intervention at each intervention 
time step when ADRIA is run. The criteria are weighted according to the user's preference and the weights can 
be defined in ADRIA through the input scenario parameters, or inline, as in the following example:

```julia
# Loading data package.
dom = ADRIA.load_domain("path to a domain data package")

# Creating 5 scenarios based on parameter bounds using the Sobol' method.
scens = ADRIA.sample(dom, 5)

# Adjusting heat_stress and wave_stress weightings.
scens.wave_stress .= [0.7,0.6,0.1,0.4]
scens.heat_stress .= [0.0,0.2,0.1,0.7]

```
The criteria weights and parameters which the user defines in ADRIA are described in the table below.

|Parameter                    |Description                                                                                                                                               |
|-----------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------|
|```heat_stress```            |Weight for heat stress when selecting seeding or shading sites. Higher values increase priority of cooler locations when seeding or shading.              |
|```wave_stress```            |Weight for wave stress when selecting seeding or shading sites. Higher values increase priority of low wave damage when seeding or shading.               |
|```shade_connectivity```     |Weight for larval connectivity when selecting shading sites. Higher values increase priority of connectivity when shading.                                | 
|```in_seed_connectivity```   |Weight for in-degree larval connectivity when selecting seeding sites. Higher values increase priority of in-going connectivity when seeding.             |                                                          
|```out_seed_connectivity```  |Weight for out-degree larval connectivity when selecting seeding sites. Higher values increase priority of out-going connectivity when shading.           |
|```coral_cover_high```       |Weight of coral cover area covered when selecting shading sites. Higher values increase priority of coral cover when shading.                             |
|```coral_cover_low```        |Weight of coral cover area available for coral gorwth when selecting seeding sites. Higher values increase priority of space for coral cover when seeding.|
|```seed_priority```          |Weight of being a priority predecessor site to a priority site when seeding. Higher values increase importance of priority predecessors when seeding.     |
|```shade_priority```         |Weight of being a priority predecessor site to a priority site when shading. Higher values increase importance of priority predecessors when shading.     |
|```zone_seed```              |Weight of being in/connected to a site in a zone when seeding. Higher values increase importance of zoning when seeding.                                  |
|```zone_shade```             |Weight of being in/connected to a site in a zone when seeding. Higher values increase importance of zoning when seeding.                                  |
|```coral_cover_tol```        |Proportion of area of corals to be seeded allowed on a site when selecting sites for seeding. If a site has less space available it will be filtered out. |                                                                    |
|```deployed_coral_risk_tol```|Risk tolerance for heat and wave risk when seeding or shading corals. Sites with greater heat/wave risk will be filtered out                              |
|```depth_min```              |Minimum depth of sites when selecting for seeding or shading corals. Sites shallower will be filtered out.                                                |
|```depth_offset```           |depth_min+depth_offset is maximum depth of sites when seeding or shading corals. Sites deeper will be filtered out.                                       |
|```dist_thresh```            |If distance sorting is used, sites closer than (median distance between sites) - dist_thresh% will be replaced.                                           |
|```top_n```                  |Number of top ranked alternatives to select from when replacing sites due to the current selection being too close together.                              |

## Post-ranking distance sorting 
If highly ranked sites happen to be very close together in a particular domain, this can represent an issue for 
site selection. This is because the site selection algorithm will seed or shade in a very small region, 
representing a significant loss of resources if a high damage, localised event, such as a ship grounding, 
were to occur. To represent an insurance policy against this possibility, the user can also specify for ADRIA to 
replace sites which are too close together with similarly ranked but further apart sites. This distance sorting 
of the selected sites is performed if the parameter `use_dist = 1`.
The parameter `dist_thresh` represents a threshold for the acceptible distance between sites and the algorithm 
will seek to replace sites further apart than this distance until the set of sites satisfies the threshold. 
The distance sorting algorithm will only search within the first `top_n` sites for the highest ranking, suitably 
far apart set of sites. If an insufficient number of sites are far enough apart, the remaining sites are simply set 
as the highest ranking sites not already selected.

## Using site selection separately
ADRIA's site selection algorithms can also be applied separately from the ecological model, using the 
`run_site_selection()` function. The user provides data layers representing the selection criteria and 
weights to produce a set of site ranks for seeding and/or shading. An example of the function's usage is 
shown below. First an ADRIA domain is loaded, containing environmental data such as wave and heat stress 
layers and site data. Then a criteria DataFrame is created, which specifies criteria weightings and wave 
and dhw scenarios. Finally, coral cover layers for each scenario are loaded. These are used within 
`run_site_selection()` to generate reccommended intervention sites for seeding at time `ts`, for each parameter 
scenario.

```julia
using ADRIA
# Loading data package.
dom = ADRIA.load_domain("path to a domain data package")

criteria_df = ADRIA.sample(dom, 5) # get scenario dataframe

area_to_seed = 1.5 * 10^-6  # area of seeded corals in km^2
ts = 5  # time step to perform site selection at

# initial coral cover matching number of criteria samples (size = (no. criteria scens, no. of sites))
sum_cover = fill(0.1, nrow(criteria_df), nrow(dom.site_data))
ranks = ADRIA.run_site_selection(dom, criteria_df[criteria_df.guided.>0, :], sum_cover[criteria_df.guided.>0, :], area_to_seed, ts)

```