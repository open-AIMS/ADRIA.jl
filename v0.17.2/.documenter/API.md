
# ADRIA API {#ADRIA-API}

## Metrics
<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics._absolute_juveniles-Union{Tuple{T}, Tuple{YAXArrays.Cubes.YAXArray{T, N, A} where {N, A<:AbstractArray{T, N}}, AbstractVector{T}}} where T<:Real' href='#ADRIA.metrics._absolute_juveniles-Union{Tuple{T}, Tuple{YAXArrays.Cubes.YAXArray{T, N, A} where {N, A<:AbstractArray{T, N}}, AbstractVector{T}}} where T<:Real'><span class="jlbinding">ADRIA.metrics._absolute_juveniles</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
absolute_juveniles(X::AbstractArray{T,3}, coral_spec::DataFrame, k_area::AbstractVector{T})::AbstractArray{T,2} where {T<:Real}
absolute_juveniles(rs::ResultSet)::AbstractArray{<:Real,2}
```


Juvenile coral cover in m².

**Arguments**
- `X` : Raw model results for a single scenario. Dimensions (n_timesteps, n_group, n_sizes,
  

n_locations)
- `coral_spec` : Coral spec DataFrame
  
- `k_area` : The coral habitable area.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/metrics.jl#L345-L356" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics._absolute_shelter_volume-Union{Tuple{T}, Tuple{YAXArrays.Cubes.YAXArray{T, 4, A} where A<:AbstractArray{T, 4}, Vector{T}, DataFrames.DataFrameRow}} where T<:Real' href='#ADRIA.metrics._absolute_shelter_volume-Union{Tuple{T}, Tuple{YAXArrays.Cubes.YAXArray{T, 4, A} where A<:AbstractArray{T, 4}, Vector{T}, DataFrames.DataFrameRow}} where T<:Real'><span class="jlbinding">ADRIA.metrics._absolute_shelter_volume</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
absolute_shelter_volume(X::YAXArray{T,4}, k_area::Vector{T}, inputs::DataFrameRow)::AbstractArray{T} where {T<:Real}
absolute_shelter_volume(X::YAXArray{T,4}, k_area::Vector{T}, inputs::YAXArray)::AbstractArray{T} where {T<:Real}
absolute_shelter_volume(X::YAXArray{T,5}, k_area::Vector{T}, inputs::DataFrame)::AbstractArray{T} where {T<:Real}
absolute_shelter_volume(X::YAXArray{T,5}, k_area::Vector{T}, inputs::YAXArray)::AbstractArray{T} where {T<:Real}
absolute_shelter_volume(rs::ResultSet)
```


Provide indication of shelter volume in volume of cubic meters.

The metric applies log-log linear models developed by Urbina-Barreto et al., [1] which uses colony diameter and planar area (2D metrics) to estimate shelter volume (a 3D metric).

**Arguments**
- `X` : raw results
  
- `k_area` : area in m^2 for each site
  
- `max_cover` : maximum possible coral cover for each site (in percentage of loc_area)
  
- `inputs` : DataFrame of scenario inputs
  

**References**
1. Urbina-Barreto, I., Chiroleu, F., Pinel, R., Fréchon, L., Mahamadaly, V.,   Elise, S., Kulbicki, M., Quod, J.-P., Dutrieux, E., Garnier, R.,   Henrich Bruggemann, J., Penin, L., &amp; Adjeroud, M. (2021). Quantifying the shelter capacity of coral reefs using photogrammetric   3D modeling: From colonies to reefscapes. Ecological Indicators, 121, 107151. https://doi.org/10.1016/j.ecolind.2020.107151
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/metrics.jl#L531-L559" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics._collate_ranked_locs-Tuple{YAXArrays.Cubes.YAXArray}' href='#ADRIA.metrics._collate_ranked_locs-Tuple{YAXArrays.Cubes.YAXArray}'><span class="jlbinding">ADRIA.metrics._collate_ranked_locs</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_collate_ranked_locs(data::YAXArray)::Matrix{Int64}
```


Collates number of ranked locations.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/ranks.jl#L78-L82" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics._collate_ranks-Tuple{ADRIA.ResultSet, Any}' href='#ADRIA.metrics._collate_ranks-Tuple{ADRIA.ResultSet, Any}'><span class="jlbinding">ADRIA.metrics._collate_ranks</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_collate_ranks(rs, selected)
```


Collates ranks into seed/fog ranking results into a common structure.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/ranks.jl#L10-L14" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics._coral_diversity-Union{Tuple{YAXArrays.Cubes.YAXArray{T, N, A} where {N, A<:AbstractArray{T, N}}}, Tuple{T}} where T<:Real' href='#ADRIA.metrics._coral_diversity-Union{Tuple{YAXArrays.Cubes.YAXArray{T, N, A} where {N, A<:AbstractArray{T, N}}}, Tuple{T}} where T<:Real'><span class="jlbinding">ADRIA.metrics._coral_diversity</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
coral_diversity(ce::AbstractArray{T})::AbstractArray{T} where {T}
coral_diversity(rs::ResultSet)::AbstractArray{T} where {T}
```


Calculates coral diversity metric as the Gini-Simpson index. This is calculated from coral evenness (which is the inverse Simpson&#39;s index, `1/D`) as `1 - 1/evenness`, which is equivalent to `1 - D`.

**Arguments**
- `ce` : Coral evenness (inverse Simpson&#39;s index).
  
- `rs` : A ResultSet object.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/metrics.jl#L498-L509" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics._coral_evenness-Union{Tuple{AbstractArray{T, 3}}, Tuple{T}} where T<:Real' href='#ADRIA.metrics._coral_evenness-Union{Tuple{AbstractArray{T, 3}}, Tuple{T}} where T<:Real'><span class="jlbinding">ADRIA.metrics._coral_evenness</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
coral_evenness(r_taxa_cover::AbstractArray{T})::AbstractArray{T} where {T<:Real}
coral_evenness(rs::ResultSet)::AbstractArray{T} where {T}
```


Calculates evenness across functional coral groups in ADRIA as a diversity metric. Inverse Simpsons diversity indicator.

**References**
1. Hill, M. O. (1973).
  

Diversity and Evenness: A Unifying Notation and Its Consequences. Ecology, 54(2), 427-432. https://doi.org/10.2307/1934352


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/metrics.jl#L459-L471" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics._extract_axes_values-Tuple{YAXArrays.Cubes.YAXArray}' href='#ADRIA.metrics._extract_axes_values-Tuple{YAXArrays.Cubes.YAXArray}'><span class="jlbinding">ADRIA.metrics._extract_axes_values</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Helper method to extract pairs of YAXArray axes (names and their values).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/metrics.jl#L76-L78" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics._get_ranks-Tuple{ADRIA.ResultSet, Symbol}' href='#ADRIA.metrics._get_ranks-Tuple{ADRIA.ResultSet, Symbol}'><span class="jlbinding">ADRIA.metrics._get_ranks</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_get_ranks(rs::ResultSet, intervention::Int64; kwargs...)
```


Extracts results for a specific intervention (:seed or :fog)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/ranks.jl#L1-L5" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics._juvenile_indicator-Union{Tuple{T}, Tuple{AbstractArray{T}, DataFrames.DataFrame, Vector{Float64}}} where T<:Real' href='#ADRIA.metrics._juvenile_indicator-Union{Tuple{T}, Tuple{AbstractArray{T}, DataFrames.DataFrame, Vector{Float64}}} where T<:Real'><span class="jlbinding">ADRIA.metrics._juvenile_indicator</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
juvenile_indicator(X::AbstractArray{T}, coral_spec::DataFrame, k_area::Vector{Float64})::AbstractArray{T,2} where {T<:Real}
juvenile_indicator(rs::ResultSet)::AbstractArray{<:Real,2}
```


Indicator for juvenile density (0 - 1), where 1 indicates the maximum theoretical density for juveniles have been achieved.

**Arguments**
- `X` : Raw model results for a single scenario. Dimensions (n_timesteps, n_group, n_sizes,
  

n_locations).
- `coral_spec` : Coral spec DataFrame.
  
- `k_area` : The coral habitable area.
  
- `max_juvenile_density` : Maximum density of juveniles defaulting to 51.8 juveniles / m²
  

**Notes**

Maximum density is 51.8 juveniles / m², where juveniles are defined as &lt; 5cm diameter. See email correspondence from: Dr. A Thompson; to: Dr. K. Anthony Subject: RE: Max density of juvenile corals on the GBR Sent: Friday, 14 October 2022 2:58 PM


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/metrics.jl#L394-L414" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics._max_juvenile_area' href='#ADRIA.metrics._max_juvenile_area'><span class="jlbinding">ADRIA.metrics._max_juvenile_area</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
_max_juvenile_area(coral_params::DataFrame, max_juv_density::Float64=51.8)
```


Calculate the maximum possible area that can be covered by juveniles for a given m².


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/metrics.jl#L382-L386" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics._reef_biodiversity_condition_index-Tuple{AbstractArray{<:Real, 3}, AbstractArray{<:Real, 3}, AbstractArray{<:Real, 3}}' href='#ADRIA.metrics._reef_biodiversity_condition_index-Tuple{AbstractArray{<:Real, 3}, AbstractArray{<:Real, 3}, AbstractArray{<:Real, 3}}'><span class="jlbinding">ADRIA.metrics._reef_biodiversity_condition_index</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_reef_biodiviersity_condition_index(rs::ResultSet)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/reef_indices.jl#L329-L331" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics._reef_condition_index-Tuple{AbstractArray{<:Real, 3}, AbstractArray{<:Real, 3}, AbstractArray{<:Real, 3}}' href='#ADRIA.metrics._reef_condition_index-Tuple{AbstractArray{<:Real, 3}, AbstractArray{<:Real, 3}, AbstractArray{<:Real, 3}}'><span class="jlbinding">ADRIA.metrics._reef_condition_index</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
reef_condition_index(ltmp_cover::AbstractArray, sv::AbstractArray, juves::AbstractArray,)::AbstractArray
reef_condition_index(ltmp_cover::AbstractArray, juves::AbstractArray, sv::AbstractArray, rubble::AbstractArray)::AbstractArray
reef_condition_index(rs::ResultSet)::AbstractArray{<:Real}
reef_condition_index(rs::ResultSet, rubble::AbstractArray)::AbstractArray{<:Real}
```


Estimates a Reef Condition Index (RCI) using either the 3-metric version using relative cover, juveniles, shelter volume or the 4-metric versions with rubble added.

The RCI is a single value that indicates the condition of a reef.

**Notes**

Juveniles are made relative to maximum observed juvenile density (15.0/m²) See table 1 in reference 1.

**Arguments**
- `ltmp_cover` : LTMP coral cover across all groups
  
- `juves` : Abundance of coral juveniles &lt; 5 cm diameter
  
- `sv` : Shelter volume based on coral sizes and abundances
  
- `rubble` : Cover of rubble (optional)
  
- `rs` : A ResultSet object.
  

**Returns**

YAXArray[timesteps ⋅ locations ⋅ scenarios]

**References**
1. Ryan F. Heneghan, Gabriela Scheufele, Yves-Marie Bozec et al. A framework to inform economic valuation of non-use benefits from coral-reef intervention efforts, 02 October 2025, PREPRINT (Version 1) available at Research Square [https://doi.org/10.21203/rs.3.rs-7644150/v1]
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/reef_indices.jl#L3-L33" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics._reef_fish_index-Union{Tuple{AbstractArray{T}}, Tuple{T}} where T<:Real' href='#ADRIA.metrics._reef_fish_index-Union{Tuple{AbstractArray{T}}, Tuple{T}} where T<:Real'><span class="jlbinding">ADRIA.metrics._reef_fish_index</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
reef_fish_index(rc::AbstractArray)
reef_fish_index(rs::ResultSet)
```


The Reef Fish Index (RFI) estimates fish biomass from relative coral cover.

A linear regression (developed by Dr. R. Heneghan, Queensland University of Technology) is used to indicate the relationship between coral cover and fish biomass. The regression was developed with digitized data from Figures 4a and 6b in Graham &amp; Nash (2013; see [1]).

Values are provided ∈ [0, 1], where 1 indicates maximum fish biomass.

Note: Coral cover here is relative to coral habitable area ($k$ area).

**Arguments**
- `rc` : Relative cover
  

**Returns**

YAXArray[timesteps ⋅ locations ⋅ scenarios], values in kg/km²

**References**
1. Graham, N.A.J., Nash, K.L., 2013.
  

The importance of structural complexity in coral reef ecosystems. Coral Reefs 32, 315–326. https://doi.org/10.1007/s00338-012-0984-y


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/reef_indices.jl#L249-L275" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics._reef_tourism_index-Tuple{AbstractArray{<:Real, 3}, AbstractArray{<:Real, 3}, AbstractArray{<:Real, 3}, AbstractArray{<:Real, 3}}' href='#ADRIA.metrics._reef_tourism_index-Tuple{AbstractArray{<:Real, 3}, AbstractArray{<:Real, 3}, AbstractArray{<:Real, 3}, AbstractArray{<:Real, 3}}'><span class="jlbinding">ADRIA.metrics._reef_tourism_index</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
reef_tourism_index(rc::AbstractArray{<:Real,3}, sv::AbstractArray{<:Real,3}, juves::AbstractArray{<:Real,3}, cots::AbstractArray{<:Real,3}, rubble::AbstractArray{<:Real,3})::AbstractArray
reef_tourism_index(rc::AbstractArray{<:Real,3}, ce::AbstractArray{<:Real,3}, sv::AbstractArray{<:Real,3}, juves::AbstractArray{<:Real,3})::AbstractArray
reef_tourism_index(rs::ResultSet, cots::YAXArray, rubble::YAXArray)::AbstractArray
reef_tourism_index(rs::ResultSet)::AbstractArray
```


Estimate tourism index.

This metric is a variation of the Reef Condition Index, but weighted by metrics known to be of importance to tourists. This version uses 5 metrics: relative cover, shelter volume, juvenile abundance, CoTS, and rubble.

**Arguments**
- `rs` : ResultSet
  
- `cots` : Outbreak status of Crown-of-Thorns Starfish
  
- `rubble` : Cover of rubble
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/reef_indices.jl#L143-L160" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics._relative_cover-Tuple{YAXArrays.Cubes.YAXArray{var"#s504", 4, A} where {var"#s504"<:Real, A<:AbstractArray{var"#s504", 4}}}' href='#ADRIA.metrics._relative_cover-Tuple{YAXArrays.Cubes.YAXArray{var"#s504", 4, A} where {var"#s504"<:Real, A<:AbstractArray{var"#s504", 4}}}'><span class="jlbinding">ADRIA.metrics._relative_cover</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
relative_cover(X::AbstractArray{<:Real}, loc_area::AbstractVector{<:Real})::AbstractArray{<:Real}
relative_cover(rs::ResultSet)::AbstractArray{<:Real}
```


Indicate coral cover relative to available hard substrate ($k$ area).

**Arguments**
- `X` : Matrix with dimensions (n_timesteps, n_functional_groups * n_size_classes,
  

n_locations) of raw model results (coral cover relative to available space)

**Returns**

Coral cover [0 - 1], relative to available $k$ area for the entire study area.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/metrics.jl#L87-L99" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics._relative_juveniles-Union{Tuple{YAXArrays.Cubes.YAXArray{T, 4, A} where A<:AbstractArray{T, 4}}, Tuple{T}} where T<:Real' href='#ADRIA.metrics._relative_juveniles-Union{Tuple{YAXArrays.Cubes.YAXArray{T, 4, A} where A<:AbstractArray{T, 4}}, Tuple{T}} where T<:Real'><span class="jlbinding">ADRIA.metrics._relative_juveniles</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
relative_juveniles(X::AbstractArray{T,3}, coral_spec::DataFrame)::AbstractArray{T,2} where {T<:Real}
relative_juveniles(rs::ResultSet)::AbstractArray{<:Real,2}
```


Juvenile coral cover relative to the location&#39;s area.

**Arguments**
- `X` : Raw model results for a single scenario. Dimensions (n_timesteps, n_group, n_sizes,
  

n_locations)
- `coral_spec` : Coral spec DataFrame
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/metrics.jl#L299-L309" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics._relative_loc_taxa_cover-Union{Tuple{AbstractArray{T, 4}}, Tuple{T}} where T<:Real' href='#ADRIA.metrics._relative_loc_taxa_cover-Union{Tuple{AbstractArray{T, 4}}, Tuple{T}} where T<:Real'><span class="jlbinding">ADRIA.metrics._relative_loc_taxa_cover</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
relative_loc_taxa_cover(X::AbstractArray{T}, k_area::Vector{T}, n_groups::Int64)::AbstractArray{T,3} where {T<:Real}
```


**Arguments**
- `X` : Raw model results for a single scenario. Dimensions (n_timesteps, n_group, n_sizes,
  

n_locations)
- `k_area` : The coral habitable area.
  
- `n_groups` : Number of function coral groups.
  

**Returns**

Coral cover, grouped by taxa for the given scenario, for each timestep and location, relative to location k area.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/metrics.jl#L266-L278" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics._relative_shelter_volume-Union{Tuple{T}, Tuple{YAXArrays.Cubes.YAXArray{T, 4, A} where A<:AbstractArray{T, 4}, Vector{T}, DataFrames.DataFrameRow}} where T<:Real' href='#ADRIA.metrics._relative_shelter_volume-Union{Tuple{T}, Tuple{YAXArrays.Cubes.YAXArray{T, 4, A} where A<:AbstractArray{T, 4}, Vector{T}, DataFrames.DataFrameRow}} where T<:Real'><span class="jlbinding">ADRIA.metrics._relative_shelter_volume</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
relative_shelter_volume(X::AbstractArray{T,4}, k_area::Vector{T}, inputs::DataFrameRow)::AbstractArray{T} where {T<:Real}
relative_shelter_volume(X::AbstractArray{T,4}, k_area::Vector{T}, inputs::YAXArray)::AbstractArray{T} where {T<:Real}
relative_shelter_volume(X::AbstractArray{T,5}, k_area::Vector{T}, inputs::DataFrame)::AbstractArray{T} where {T<:Real}
relative_shelter_volume(X::AbstractArray{T,5}, k_area::Vector{T}, inputs::YAXArray)::AbstractArray{T} where {T<:Real}
relative_shelter_volume(rs::ResultSet)
```


Provide indication of shelter volume relative to theoretical maximum volume for the area covered by coral.

The metric applies log-log linear models developed by Urbina-Barreto et al., [1] which uses colony diameter and planar area (2D metrics) to estimate shelter volume (a 3D metric).

$$RSV = \begin{cases}
TASV / MSV & MSV > 0, \\
0 & \text{otherwise}
\end{cases}$$

where $TASV$ represents Total Absolute Shelter Volume and $MSV$ represents the maximum shelter volume possible.

**Arguments**
- `X` : raw results
  
- `k_area` : area in m^2 for each site
  
- `scens` : DataFrame of scenario inputs
  

**Returns**

Shelter volume relative to a theoretical maximum volume for the available $k$ area. The maximum volume is defined as the volume occupied by corals where there is 1 95cm diameter Tabular Acropora per m².

**References**
1. Urbina-Barreto, I., Chiroleu, F., Pinel, R., Fréchon, L., Mahamadaly, V., Elise, S.,
  

Kulbicki, M., Quod, J.-P., Dutrieux, E., Garnier, R., Henrich Bruggemann, J., Penin, L., &amp; Adjeroud, M. (2021). Quantifying the shelter capacity of coral reefs using photogrammetric 3D modeling: From colonies to reefscapes. Ecological Indicators, 121, 107151. https://doi.org/10.1016/j.ecolind.2020.107151


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/metrics.jl#L666-L705" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics._relative_taxa_cover-Tuple{AbstractArray{<:Real, 4}, Vector{<:Real}}' href='#ADRIA.metrics._relative_taxa_cover-Tuple{AbstractArray{<:Real, 4}, Vector{<:Real}}'><span class="jlbinding">ADRIA.metrics._relative_taxa_cover</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
relative_taxa_cover(X::AbstractArray{<:Real}, k_area::Vector{<:Real}, n_groups::Int64)::AbstractArray{<:Real,2}
relative_taxa_cover(rs::ResultSet)::AbstractArray{<:Real,2}
```


Relative coral cover grouped by groups summed up across all locations.

**Arguments**
- `X` : Raw model results for a single scenario. Dimensions (n_timesteps, n_group, n_sizes,
  

n_locations).
- `k_area` : The coral habitable area.
  
- `n_groups` : Number of function coral groups.
  

**Returns**

Coral cover, grouped by taxa for the given scenario, summed up across all locations, relative to total k area.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/metrics.jl#L191-L206" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics._scenario_absolute_juveniles-Tuple{YAXArrays.Cubes.YAXArray}' href='#ADRIA.metrics._scenario_absolute_juveniles-Tuple{YAXArrays.Cubes.YAXArray}'><span class="jlbinding">ADRIA.metrics._scenario_absolute_juveniles</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
scenario_absolute_juveniles(data::YAXArray, coral_spec::DataFrame, k_area::AbstractVector{<:Real}; kwargs...)::AbstractArray{<:Real}
scenario_absolute_juveniles(rs::ResultSet; kwargs...)::AbstractArray{<:Real}
```


Calculate the mean absolute juvenile population for each scenario for the entire domain.

**Arguments**
- `aj` : Raw data for a single scenario.
  
- `k_area` : K_area.
  
- `rs` : Resultset.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/scenario.jl#L185-L195" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics._scenario_asv-Tuple{YAXArrays.Cubes.YAXArray}' href='#ADRIA.metrics._scenario_asv-Tuple{YAXArrays.Cubes.YAXArray}'><span class="jlbinding">ADRIA.metrics._scenario_asv</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
scenario_asv(sv::YAXArray; kwargs...)::AbstractArray{<:Real}
scenario_asv(rs::ResultSet; kwargs...)::AbstractArray{<:Real}
```


Calculate the mean absolute shelter volumes for each scenario for the entire domain.

**Arguments**
- `asv` : Absolute shelter volume.
  
- `rs` : Resultset.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/scenario.jl#L243-L252" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics._scenario_evenness-Tuple{YAXArrays.Cubes.YAXArray}' href='#ADRIA.metrics._scenario_evenness-Tuple{YAXArrays.Cubes.YAXArray}'><span class="jlbinding">ADRIA.metrics._scenario_evenness</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
scenario_evenness(ev::YAXArray; kwargs...)::AbstractArray{<:Real}
scenario_evenness(rs::ResultSet; kwargs...)::AbstractArray{<:Real}
```


Calculate the mean coral evenness for each scenario for the entire domain.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/scenario.jl#L294-L299" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics._scenario_juvenile_indicator-Tuple{YAXArrays.Cubes.YAXArray{var"#s504", 3, A} where {var"#s504"<:Real, A<:AbstractArray{var"#s504", 3}}}' href='#ADRIA.metrics._scenario_juvenile_indicator-Tuple{YAXArrays.Cubes.YAXArray{var"#s504", 3, A} where {var"#s504"<:Real, A<:AbstractArray{var"#s504", 3}}}'><span class="jlbinding">ADRIA.metrics._scenario_juvenile_indicator</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
scenario_juvenile_indicator(data::YAXArray, coral_spec::DataFrame, k_area::AbstractVector{<:Real}; kwargs...)::AbstractArray{<:Real}
scenario_juvenile_indicator(rs::ResultSet; kwargs...)::AbstractArray{<:Real}
```


Determine juvenile indicator ∈ [0, 1], where 1 indicates maximum mean juvenile density (51.8) has been achieved.

**Arguments**
- `ji` : Juvenile Indicator for each location.
  
- `rs` : Resultset.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/scenario.jl#L216-L225" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics._scenario_rci-Tuple{YAXArrays.Cubes.YAXArray, YAXArrays.Cubes.YAXArray}' href='#ADRIA.metrics._scenario_rci-Tuple{YAXArrays.Cubes.YAXArray, YAXArrays.Cubes.YAXArray}'><span class="jlbinding">ADRIA.metrics._scenario_rci</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
scenario_rci(rci::YAXArray, tac::YAXArray; kwargs...)
scenario_rci(rci::YAXArray, rubble::YAXArray; kwargs...)
scenario_rci(rs::ResultSet; kwargs...)
```


Extract the total populated area of locations with Reef Condition Index of &quot;Good&quot; or higher for each scenario for the entire domain.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/reef_indices.jl#L101-L108" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics._scenario_relative_cover-Tuple{ADRIA.ResultSet}' href='#ADRIA.metrics._scenario_relative_cover-Tuple{ADRIA.ResultSet}'><span class="jlbinding">ADRIA.metrics._scenario_relative_cover</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
scenario_relative_cover(rs::ResultSet; kwargs...)::AbstractArray{<:Real}
```


Calculate the mean relative coral cover for each scenario for the entire domain.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/scenario.jl#L61-L65" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics._scenario_relative_juveniles-Tuple{YAXArrays.Cubes.YAXArray{var"#s503", 3, A} where {var"#s503"<:Real, A<:AbstractArray{var"#s503", 3}}, AbstractVector{<:Real}}' href='#ADRIA.metrics._scenario_relative_juveniles-Tuple{YAXArrays.Cubes.YAXArray{var"#s503", 3, A} where {var"#s503"<:Real, A<:AbstractArray{var"#s503", 3}}, AbstractVector{<:Real}}'><span class="jlbinding">ADRIA.metrics._scenario_relative_juveniles</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
scenario_relative_juveniles(X::YAXArray{<:Real,3}, coral_spec::DataFrame, k_area::AbstractVector{<:Real}; kwargs...)::AbstractArray{<:Real}
scenario_relative_juveniles(rs::ResultSet; kwargs...)::YAXArray
```


Calculate the mean relative juvenile population for each scenario for the entire domain.

**Arguments**
- `X` : Raw data for a single scenario.
  
- `rs` : Resultset.
  
- `coral_spec` : Coral spec DataFrame.
  
- `k_area` : K_area.
  

**Examples**

```julia
num_scens = 2^5
scens = ADRIA.sample(dom, num_scens)

_coral_spec = ADRIA.to_coral_spec(scens[1,:])
_k_area = loc_k_area(dom)

# X contains raw coral cover results for a single scenario
ADRIA.metrics.scenario_relative_juveniles(X, _coral_spec, _k_area)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/scenario.jl#L140-L163" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics._scenario_rfi-Tuple{YAXArrays.Cubes.YAXArray}' href='#ADRIA.metrics._scenario_rfi-Tuple{YAXArrays.Cubes.YAXArray}'><span class="jlbinding">ADRIA.metrics._scenario_rfi</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



scenario_rfi(rfi::YAXArray; kwargs...) scenario_rfi(rs::ResultSet; kwargs...)

Calculate the mean Reef Fish Index (RFI) for each scenario for the entire domain.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/reef_indices.jl#L308-L313" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics._scenario_rsv-Tuple{YAXArrays.Cubes.YAXArray}' href='#ADRIA.metrics._scenario_rsv-Tuple{YAXArrays.Cubes.YAXArray}'><span class="jlbinding">ADRIA.metrics._scenario_rsv</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
scenario_rsv(sv::YAXArray; kwargs...)::AbstractArray{<:Real}
scenario_rsv(rs::ResultSet; kwargs...)::AbstractArray{<:Real}
```


Calculate the mean relative shelter volumes for each scenario for the entire domain.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/scenario.jl#L271-L276" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics._scenario_rti-Tuple{YAXArrays.Cubes.YAXArray}' href='#ADRIA.metrics._scenario_rti-Tuple{YAXArrays.Cubes.YAXArray}'><span class="jlbinding">ADRIA.metrics._scenario_rti</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
scenario_rti(rs::ResultSet, cots::YAXArray, rubble::YAXArray; kwargs...)
scenario_rti(rs::ResultSet; kwargs)
```


Calculate the mean Reef Tourism Index (RTI) for each scenario for the entire domain.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/reef_indices.jl#L224-L229" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics._scenario_total_cover-Tuple{AbstractArray}' href='#ADRIA.metrics._scenario_total_cover-Tuple{AbstractArray}'><span class="jlbinding">ADRIA.metrics._scenario_total_cover</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
scenario_total_cover(rs::ResultSet; kwargs...)::AbstractArray{<:Real}
```


Calculate the mean absolute coral for each scenario for the entire domain.

**Arguments**
- `tac` : Total absolute cover
  
- `rs` : ResultSet
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/scenario.jl#L33-L41" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics._summarize_inner-Union{Tuple{NK}, Tuple{NA}, Tuple{F}, Tuple{T}, Tuple{AbstractArray{T}, F, NTuple{NK, Base.OneTo{Int64}}, Val{NA}, Val{NK}}} where {T, F<:Function, NA, NK}' href='#ADRIA.metrics._summarize_inner-Union{Tuple{NK}, Tuple{NA}, Tuple{F}, Tuple{T}, Tuple{AbstractArray{T}, F, NTuple{NK, Base.OneTo{Int64}}, Val{NA}, Val{NK}}} where {T, F<:Function, NA, NK}'><span class="jlbinding">ADRIA.metrics._summarize_inner</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Barrier function: Val{NA}/Val{NK} allow Julia to specialize on ndims at compile time, eliminating per-iteration dynamic dispatch on view/metric calls.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/spatial.jl#L74-L77" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics._total_absolute_cover-Tuple{YAXArrays.Cubes.YAXArray{var"#s504", 3, A} where {var"#s504"<:Real, A<:AbstractArray{var"#s504", 3}}, Vector{<:Real}}' href='#ADRIA.metrics._total_absolute_cover-Tuple{YAXArrays.Cubes.YAXArray{var"#s504", 3, A} where {var"#s504"<:Real, A<:AbstractArray{var"#s504", 3}}, Vector{<:Real}}'><span class="jlbinding">ADRIA.metrics._total_absolute_cover</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
total_absolute_cover(relative_cover::AbstractArray{<:Real}, k_area::Vector{<:Real})::AbstractArray{<:Real}
total_absolute_cover(rs::ResultSet)::AbstractArray{<:Real}
```


The Total Absolute Coral Cover. Sum of proportional area taken up by all corals, multiplied by the location area.

**Arguments**
- `relative_cover` : Array with relative_cover
  
- `k_area` : Proportional area, with locations following the same order as given indicated in `relative_cover`.
  

**Returns**

Absolute coral cover for a given location in m².


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/metrics.jl#L157-L170" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics._years_above_threshold-Tuple{YAXArrays.Cubes.YAXArray{var"#s500", 2, A} where {var"#s500"<:Real, A<:AbstractMatrix{var"#s500"}}}' href='#ADRIA.metrics._years_above_threshold-Tuple{YAXArrays.Cubes.YAXArray{var"#s500", 2, A} where {var"#s500"<:Real, A<:AbstractMatrix{var"#s500"}}}'><span class="jlbinding">ADRIA.metrics._years_above_threshold</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
years_above_threshold(data::YAXArray{<:Real,2}; threshold::Real=0.1, kwargs...)::AbstractArray{<:Real}
years_above_threshold(data::YAXArray{<:Real,3}; threshold::Real=0.1, kwargs...)::AbstractArray{<:Real}
```


Count the number of years a metric exceeded a threshold.

For 2D input (timesteps x scenarios): returns Vector[scenarios]. For 3D input (timesteps x locations x scenarios): returns Matrix[locations x scenarios].

**Arguments**
- `data` : Metric values over time.
  
- `threshold` : Threshold value (default: 0.1).
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/scenario.jl#L80-L92" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics.axes_units-Tuple{Union{Tuple, Vector{Symbol}}}' href='#ADRIA.metrics.axes_units-Tuple{Union{Tuple, Vector{Symbol}}}'><span class="jlbinding">ADRIA.metrics.axes_units</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
axes_units(axes_names::Union{Vector{Symbol},Tuple})::Tuple
```


Get units for each metric axis.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/metadata.jl#L76-L80" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics.call_metric-Tuple{Union{Function, ADRIA.metrics.Metric}, YAXArrays.Cubes.YAXArray, Vararg{Any}}' href='#ADRIA.metrics.call_metric-Tuple{Union{Function, ADRIA.metrics.Metric}, YAXArrays.Cubes.YAXArray, Vararg{Any}}'><span class="jlbinding">ADRIA.metrics.call_metric</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
call_metric(metric::Union{Function,Metric}, data::YAXArray, args...; kwargs...)
```


Convenience method that slices the data in the specified manner.

**Arguments**
- `metric` : Function, the metric function to apply to &quot;raw&quot; data.
  
- `data` : YAXArray, data to pass into `metric`
  
- `args` : Additional positional arguments to pass into `metric`
  
- `kwargs` : Additional keyword arguments to pass into `slice_results`
  - `dims` : dummy keyword argument, not used but defined to allow use with other methods
    
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/utils.jl#L65-L76" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics.deployed_locations-Tuple{ADRIA.ResultSet}' href='#ADRIA.metrics.deployed_locations-Tuple{ADRIA.ResultSet}'><span class="jlbinding">ADRIA.metrics.deployed_locations</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
deployed_locations(rs::ResultSet; intervention::Symbol=:seed)::Vector{Int}
```


Return the integer indices of locations that received at least one deployment of the given intervention type, across all timesteps and scenarios.

**Arguments**
- `rs`           : ResultSet
  
- `intervention` : Intervention type — one of `:seed`, `:fog`, `:mc` (default: `:seed`)
  

**Returns**

Vector of 1-based integer location indices. Pass directly to metrics that accept a `locations` keyword argument, e.g.:

```julia
locs = ADRIA.metrics.deployed_locations(rs; intervention=:seed)
rc   = ADRIA.metrics.scenario_relative_cover(rs; locations=locs)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/ranks.jl#L229-L247" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics.dims-Tuple{ADRIA.metrics.Metric}' href='#ADRIA.metrics.dims-Tuple{ADRIA.metrics.Metric}'><span class="jlbinding">ADRIA.metrics.dims</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
dims(m::Metric)::Tuple
```


Get dimension names for a given outcome/metric.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/utils.jl#L47-L51" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics.dominates-Tuple{AbstractVector{<:Real}, AbstractVector{<:Real}}' href='#ADRIA.metrics.dominates-Tuple{AbstractVector{<:Real}, AbstractVector{<:Real}}'><span class="jlbinding">ADRIA.metrics.dominates</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
dominates(x::Vector{<:Real}, y::Vector{<:Real})::Vector
```


Adapted from: https://discourse.julialang.org/t/fast-optimized-non-dominated-sorting-algorithms/86793/7

Original function name is `dominates2()`


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/pareto.jl#L3-L10" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics.ensemble_loc_difference-Union{Tuple{T}, Tuple{YAXArrays.Cubes.YAXArray{T, 3, A} where A<:AbstractArray{T, 3}, DataFrames.DataFrame}} where T' href='#ADRIA.metrics.ensemble_loc_difference-Union{Tuple{T}, Tuple{YAXArrays.Cubes.YAXArray{T, 3, A} where A<:AbstractArray{T, 3}, DataFrames.DataFrame}} where T'><span class="jlbinding">ADRIA.metrics.ensemble_loc_difference</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
ensemble_loc_difference(outcome::YAXArray{T,3}, scens::DataFrame; agg_metric::Union{Function,AbstractFloat}=median, diff_target=:guided, conf::Float64=0.95, rng::AbstractRNG=Random.GLOBAL_RNG)::YAXArray where {T}
```


Mean bootstrapped difference (counterfactual - target) between some outcome aggregated for each location.

**Arguments**
- `outcome` : Metric outcome with dimensions (:timesteps, :locations, :scenarios).
  
- `scens` : Scenarios DataFrame.
  
- `agg_metric` : Metric used to aggregate scenarios when comparing between counterfactual and
  

target. If it is an `AbstractFloat` between 0 and 1, it uses the `agg_metric`-th quantile. Defaults to `median`.
- `diff_target` : Target group of scenarios to compare with. Valid options are `:guided` and
  

`:unguided`. Defaults to `:guided`
- `conf` : Percentile used for the confidence interval. Defaults to 0.95.
  
- `rng` : Pseudorandom number generator.
  

**Example**

```julia
# Load domain
dom = ADRIA.load_domain(path_to_domain, "<RCP>")

# Create scenarios
num_scens = 2^6
scens = ADRIA.sample(dom, num_scens)

# Run model
rs = ADRIA.run_scenarios(dom, scens, "45")

# Calculate difference to the counterfactual for given metric
_relative_cover = metrics.relative_cover(rs)

# Compute difference between guided and counterfactual using the 0.6-th quantile
gd_res = metrics.ensemble_loc_difference(r_cover, scens; agg_metric=0.6)

# Compute difference between unguided and counterfactual using the median
ug_res = metrics.ensemble_loc_difference(r_cover, scens; diff_target=:unguided)

# Plot maps of difference to the counterfactual
ADRIA.viz.map(rs, gd_res[summary=At(:agg_value)]; diverging=true)
ADRIA.viz.map(rs, ug_res[summary=At(:agg_value)]; diverging=true)
```


**Returns**

Vector with bootstrapped difference (counterfactual - guided) for each location.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/spatial.jl#L157-L202" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics.fill_axes_metadata!-Tuple{YAXArrays.Cubes.YAXArray}' href='#ADRIA.metrics.fill_axes_metadata!-Tuple{YAXArrays.Cubes.YAXArray}'><span class="jlbinding">ADRIA.metrics.fill_axes_metadata!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
fill_axes_metadata!(outcomes::YAXArray)::Nothing
```


Fill outcomes axes metadata.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/metadata.jl#L60-L64" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics.fill_metadata!-Union{Tuple{A}, Tuple{N}, Tuple{T}, Tuple{YAXArrays.Cubes.YAXArray{T, N, A}, ADRIA.metrics.Metric}} where {T, N, A}' href='#ADRIA.metrics.fill_metadata!-Union{Tuple{A}, Tuple{N}, Tuple{T}, Tuple{YAXArrays.Cubes.YAXArray{T, N, A}, ADRIA.metrics.Metric}} where {T, N, A}'><span class="jlbinding">ADRIA.metrics.fill_metadata!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
fill_metadata!(outcomes::YAXArray{T,N,A}, metric::Metric)::YAXArray{T,N,A} where {T,N,A}
fill_metadata!(outcomes::YAXArray{T,N,A}, metadata::Dict{Symbol,Any})::YAXArray{T,N,A} where {T,N,A}
```


Fill outcomes YAXArray metadata (`properties` attribute).

**Arguments**
- `outcomes` : YAXArray datacube of metric outcomes.
  
- `metric` : ADRIA.metrics.Metric object.
  
- `metadata` : Dict to be used to fill outcomes metrics metadata.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/metadata.jl#L13-L23" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics.fog_ranks-Tuple{ADRIA.ResultSet}' href='#ADRIA.metrics.fog_ranks-Tuple{ADRIA.ResultSet}'><span class="jlbinding">ADRIA.metrics.fog_ranks</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
fog_ranks(rs::ResultSet; kwargs...)
```


**Arguments**
- rs : ResultSet
  
- kwargs : named dimensions to slice across
  

**Returns**

YAXArray[timesteps, sites, scenarios]

**Example**

```julia
ADRIA.metrics.fog_ranks(rs; timesteps=1:10, scenarios=3:5)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/ranks.jl#L58-L72" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics.loc_trajectory-Union{Tuple{A}, Tuple{N}, Tuple{T}, Tuple{D}, Tuple{Any, YAXArrays.Cubes.YAXArray{D, T, N, A}}} where {D, T, N, A}' href='#ADRIA.metrics.loc_trajectory-Union{Tuple{A}, Tuple{N}, Tuple{T}, Tuple{D}, Tuple{Any, YAXArrays.Cubes.YAXArray{D, T, N, A}}} where {D, T, N, A}'><span class="jlbinding">ADRIA.metrics.loc_trajectory</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
loc_trajectory(metric, data::YAXArray{D,T,N,A})::YAXArray where {D,T,N,A}
```


Alias for summarize(data, [:scenarios], metric). Collate trajectory for each location, applying `metric` across values for all scenarios.

**Examples**

```julia
using Statistics

rs = ADRIA.load_results("some results")
tac = ADRIA.metrics.total_absolute_cover(rs)

# Get median trajectory for each site
ADRIA.metrics.loc_trajectory(median, tac)
#75×216 YAXArray{Float64,2} with dimensions:
#  Dim{:timesteps} Categorical{Any} Any[1, 2, …, 74, 75] Unordered,
#  Dim{:locations} Categorical{Any} Any[1, 2, …, 215, 216] Unordered
#Total size: 126.56 KB

# Get upper 95% CI for each site
ADRIA.metrics.loc_trajectory(x -> quantile(x, 0.975), tac)
#75×216 YAXArray{Float64,2} with dimensions:
#  Dim{:timesteps} Categorical{Any} Any[1, 2, …, 74, 75] Unordered,
#  Dim{:locations} Categorical{Any} Any[1, 2, …, 215, 216] Unordered
#Total size: 126.56 KB
```


**Arguments**
- metric : Any function (nominally from the Statistics package) to be applied to `data`
  
- data : Data set to apply metric to
  

**Returns**

2D array of $T ⋅ S$, where $T$ is total number of time steps and $S$ is number of locations


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/spatial.jl#L29-L62" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics.metadata-Tuple{YAXArrays.Cubes.YAXArray}' href='#ADRIA.metrics.metadata-Tuple{YAXArrays.Cubes.YAXArray}'><span class="jlbinding">ADRIA.metrics.metadata</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
metadata(outcomes::YAXArray)::Dict{Symbol,Any}
```


Helper function to extract metadata from YAXArrays.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/metadata.jl#L4-L8" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics.metric_label-Tuple{ADRIA.metrics.Metric}' href='#ADRIA.metrics.metric_label-Tuple{ADRIA.metrics.Metric}'><span class="jlbinding">ADRIA.metrics.metric_label</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
metric_label(m::Metric)::String
metric_label(f::Function, unit::String)
```


Return name of metric in the format: &quot;Title Case [Unit]&quot;, suitable for use as a label.

**Example**

```julia
m_label = metric_label(scenario_total_cover)
# "Scenario Total Cover [m²]"
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/utils.jl#L23-L34" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics.n_fog_locations-Tuple{ADRIA.ResultSet}' href='#ADRIA.metrics.n_fog_locations-Tuple{ADRIA.ResultSet}'><span class="jlbinding">ADRIA.metrics.n_fog_locations</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
n_fog_locations(rs::ResultSet; kwargs...)::Matrix{Int64}
```


Determine the number of locations fogged at each time step, for each scenario.

**Returns**

YAXArray[timesteps ⋅ scenarios] indicating the number of locations fogged at each time step.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/ranks.jl#L119-L126" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics.n_seed_locations-Tuple{ADRIA.ResultSet}' href='#ADRIA.metrics.n_seed_locations-Tuple{ADRIA.ResultSet}'><span class="jlbinding">ADRIA.metrics.n_seed_locations</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
n_seed_locations(rs::ResultSet; kwargs...)::Matrix{Int64}
```


Determine the number of locations seeded at each time step, for each scenario.

**Returns**

YAXArray[timesteps ⋅ scenarios] indicating the number of locations seeded at each time step.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/ranks.jl#L101-L108" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics.nds' href='#ADRIA.metrics.nds'><span class="jlbinding">ADRIA.metrics.nds</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
nds(X::AbstractArray{<:Real}, dist::Int64=0)::Vector{Vector{<:Int}}
```


Naive n-dimensional non-dominated sorting.

Adapted from: https://discourse.julialang.org/t/fast-optimized-non-dominated-sorting-algorithms/86793/7

Original function name is `nds4()`

**Arguments**

X : outcomes, where rows are scenarios and columns are metric results. dist : distance from front, where 0 is on the frontier.

**Returns**

Vector of Vectors with row indices for each `dist` from frontier, where 0 is on the frontier.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/pareto.jl#L21-L37" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics.per_loc-Union{Tuple{A}, Tuple{N}, Tuple{T}, Tuple{D}, Tuple{Any, YAXArrays.Cubes.YAXArray{D, T, N, A}}} where {D, T, N, A}' href='#ADRIA.metrics.per_loc-Union{Tuple{A}, Tuple{N}, Tuple{T}, Tuple{D}, Tuple{Any, YAXArrays.Cubes.YAXArray{D, T, N, A}}} where {D, T, N, A}'><span class="jlbinding">ADRIA.metrics.per_loc</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
per_loc(metric, data::YAXArray{D,T,N,A})::YAXArray where {D,T,N,A}
```


Alias for summarize(data, [:scenarios, :timesteps], metric). Get metric results applied to the location-level at indicated time (or across timesteps).

**Arguments**
- metric : Any function (nominally from the Statistics package) to be applied to `data`
  
- data : Data set to apply metric to
  
- timesteps : timesteps to apply `metric` across
  

**Returns**

Named Vector of $N$ elements, where $N$ is the number of locations.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/spatial.jl#L6-L19" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics.scenario_outcomes-Tuple{ADRIA.ResultSet, Vector{<:ADRIA.metrics.Metric}}' href='#ADRIA.metrics.scenario_outcomes-Tuple{ADRIA.ResultSet, Vector{<:ADRIA.metrics.Metric}}'><span class="jlbinding">ADRIA.metrics.scenario_outcomes</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
scenario_outcomes(rs::ResultSet, metrics::Vector{Metric})::YAXArray
```


Get outcomes for a given list of metrics and a result set.

**Arguments**
- `rs` : ResultSet
  
- `metrics` : Vector of scenario Metrics (the ones that start with `scenario_`)
  

**Returns**

YAXArray with (:timesteps, :scenarios, :outcomes)

**Examples**

```julia
metrics::Vector{ADRIA.metrics.Metric} = [
    ADRIA.metrics.scenario_total_cover,
    ADRIA.metrics.scenario_asv,
    ADRIA.metrics.scenario_absolute_juveniles,
]

# 3-dimensional Array of outcomes
outcomes = ADRIA.metrics.scenario_outcomes(rs, metrics)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/scenario.jl#L316-L339" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics.scenario_trajectory-Tuple{AbstractArray}' href='#ADRIA.metrics.scenario_trajectory-Tuple{AbstractArray}'><span class="jlbinding">ADRIA.metrics.scenario_trajectory</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
scenario_trajectory(data::AbstractArray; metric=mean)::YAXArray{<:Real}
```


Produce scenario trajectories using the provided metric/aggregation function.

**Arguments**
- `data` : Results to aggregate
  
- `metric` : Function or Callable used to summarize data
  

**Returns**

Matrix[timesteps ⋅ scenarios]


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/scenario.jl#L12-L23" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics.seed_ranks-Tuple{ADRIA.ResultSet}' href='#ADRIA.metrics.seed_ranks-Tuple{ADRIA.ResultSet}'><span class="jlbinding">ADRIA.metrics.seed_ranks</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
seed_ranks(rs::ResultSet; kwargs...)
```


**Arguments**
- rs : ResultSet
  
- kwargs : named dimensions to slice across
  

**Returns**

YAXArray[timesteps, sites, scenarios]

**Example**

```julia
ADRIA.metrics.seed_ranks(rs; timesteps=1:10, scenarios=3:5)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/ranks.jl#L38-L52" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics.slice_results-Tuple{YAXArrays.Cubes.YAXArray}' href='#ADRIA.metrics.slice_results-Tuple{YAXArrays.Cubes.YAXArray}'><span class="jlbinding">ADRIA.metrics.slice_results</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
slice_results(data::YAXArray; timesteps=(:), species=(:), locations=(:), scenarios=(:))
```


Slice data as indicated. Dimensions not found in target data are ignored.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/utils.jl#L86-L90" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics.summarize-Union{Tuple{A}, Tuple{N}, Tuple{T}, Tuple{D}, Tuple{YAXArrays.Cubes.YAXArray{D, T, N, A}, Vector{Symbol}, Function}} where {D, T, N, A}' href='#ADRIA.metrics.summarize-Union{Tuple{A}, Tuple{N}, Tuple{T}, Tuple{D}, Tuple{YAXArrays.Cubes.YAXArray{D, T, N, A}, Vector{Symbol}, Function}} where {D, T, N, A}'><span class="jlbinding">ADRIA.metrics.summarize</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
summarize(data::YAXArray{<:Real}, alongs_axis::Vector{Symbol}, metric::Function)::YAXArray{<:Real}
summarize(data::YAXArray{<:Real}, alongs_axis::Vector{Symbol}, metric::Function, timesteps::Union{UnitRange,Vector{Int64},BitVector})::YAXArray{<:Real}
```


Apply summary metric along some axis of a data set across some or all timesteps.

**Arguments**
- `data` : Data set to apply metric to.
  
- `alongs_axis` : which axis will be replaced with (:) when slicing.
  
- `metric` : Any function (nominally from the Statistics package) to be applied to `data`.
  
- `timesteps` : timesteps to apply `metric` across.
  

**Returns**

YAXArray with summary metric for the remaining axis.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/spatial.jl#L93-L107" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics.summarize_absolute_shelter_volume-Tuple{YAXArrays.Cubes.YAXArray}' href='#ADRIA.metrics.summarize_absolute_shelter_volume-Tuple{YAXArrays.Cubes.YAXArray}'><span class="jlbinding">ADRIA.metrics.summarize_absolute_shelter_volume</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
summarize_absolute_shelter_volume(sv::YAXArray; kwargs...)::Dict{Symbol,AbstractArray{<:Real}}
summarize_absolute_shelter_volume(rs::ResultSet, kwargs...)::Dict{Symbol,AbstractArray{<:Real}}
```


Calculate summarized coral evenness.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/temporal.jl#L132-L137" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics.summarize_coral_evenness-Tuple{YAXArrays.Cubes.YAXArray}' href='#ADRIA.metrics.summarize_coral_evenness-Tuple{YAXArrays.Cubes.YAXArray}'><span class="jlbinding">ADRIA.metrics.summarize_coral_evenness</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
summarize_coral_evenness(raw::YAXArray; kwargs...)::Dict{Symbol,AbstractArray{<:Real}}
summarize_coral_evenness(rs::ResultSet, kwargs...)::Dict{Symbol,AbstractArray{<:Real}}
```


Calculate summarized coral evenness.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/temporal.jl#L114-L119" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics.summarize_raw-Tuple{YAXArrays.Cubes.YAXArray}' href='#ADRIA.metrics.summarize_raw-Tuple{YAXArrays.Cubes.YAXArray}'><span class="jlbinding">ADRIA.metrics.summarize_raw</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
summarize_raw(data::YAXArray; kwargs...)::Dict{Symbol,AbstractArray{<:Real}}
```


Summarize raw data, aggregating the specified dimensions (e.g., `timesteps`, `scenarios`, etc.) and collapsing given `dims`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/temporal.jl#L43-L48" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics.summarize_relative_cover-Tuple{YAXArrays.Cubes.YAXArray}' href='#ADRIA.metrics.summarize_relative_cover-Tuple{YAXArrays.Cubes.YAXArray}'><span class="jlbinding">ADRIA.metrics.summarize_relative_cover</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
summarize_relative_cover(rc::YAXArray; kwargs...)::Dict{Symbol,AbstractArray{<:Real}}
summarize_relative_cover(rs::ResultSet, kwargs...)::Dict{Symbol,AbstractArray{<:Real}}
```


Calculate summarized relative cover.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/temporal.jl#L96-L101" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics.summarize_relative_shelter_volume-Tuple{YAXArrays.Cubes.YAXArray}' href='#ADRIA.metrics.summarize_relative_shelter_volume-Tuple{YAXArrays.Cubes.YAXArray}'><span class="jlbinding">ADRIA.metrics.summarize_relative_shelter_volume</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
summarize_relative_shelter_volume(sv::YAXArray; kwargs...)::Dict{Symbol,AbstractArray{<:Real}}
summarize_relative_shelter_volume(rs::ResultSet, kwargs...)::Dict{Symbol,AbstractArray{<:Real}}
```


Calculate summarized coral evenness.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/temporal.jl#L150-L155" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics.summarize_total_cover-Tuple{YAXArrays.Cubes.YAXArray, AbstractArray{<:Real}}' href='#ADRIA.metrics.summarize_total_cover-Tuple{YAXArrays.Cubes.YAXArray, AbstractArray{<:Real}}'><span class="jlbinding">ADRIA.metrics.summarize_total_cover</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
summarize_total_cover(raw::YAXArray, areas::AbstractArray{<:Real}; kwargs...)::Dict{Symbol,AbstractArray{<:Real}}
summarize_total_cover(rs::ResultSet; kwargs...)::Dict{Symbol,AbstractArray{<:Real}}
```


Calculate summarized total absolute cover.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/temporal.jl#L74-L79" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics.to_string-Tuple{ADRIA.metrics.Metric}' href='#ADRIA.metrics.to_string-Tuple{ADRIA.metrics.Metric}'><span class="jlbinding">ADRIA.metrics.to_string</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
to_string(m::Metric)::String
```


Get name of metric as a string.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/utils.jl#L1-L5" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics.to_symbol-Tuple{ADRIA.metrics.Metric}' href='#ADRIA.metrics.to_symbol-Tuple{ADRIA.metrics.Metric}'><span class="jlbinding">ADRIA.metrics.to_symbol</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
to_symbol(m::Metric)::String
```


Get name of metric as a symbol.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/utils.jl#L14-L18" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics.top_N_sites-Tuple{ADRIA.ResultSet, Int64}' href='#ADRIA.metrics.top_N_sites-Tuple{ADRIA.ResultSet, Int64}'><span class="jlbinding">ADRIA.metrics.top_N_sites</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
top_N_sites(rs::ResultSet; N::Int64; metric::relative_cover)
top_N_sites(data::AbstractArray{Real}, N::Int64; stat=mean)
```


Return the top `N` sites according to the provided metric (defaulting to `mean` of `relative_cover`).

**Arguments**
- rs : ResultSet
  
- N : Number of best performing sites to be selected
  
- metric : Metric to use to order sites from best to worst,          must take ResultSet as input
  
- stat : Summary statistic to use for comparison (default: mean)
  

**Returns**

YAXArray[:scenarios, :locations], where `locations` indicates order of location ranking.

**Example**

```julia
ADRIA.metrics.top_N_sites(rs, 5)
ADRIA.metrics.top_N_sites(rs, 5; metric=ADRIA.metric.relative_cover)
ADRIA.metrics.top_N_sites(rs, 5; metric=ADRIA.metric.relative_cover, stat=median)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/ranks.jl#L186-L208" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics.top_n_seeded_sites-Tuple{ADRIA.ResultSet, Int64}' href='#ADRIA.metrics.top_n_seeded_sites-Tuple{ADRIA.ResultSet, Int64}'><span class="jlbinding">ADRIA.metrics.top_n_seeded_sites</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
top_n_seeded_sites(rs::ResultSet, n::Int64; kwargs...)
```


Get the top n seeded sites over time by their unique location id. Lower rank values are better (e.g., 1 = first choice)

**Arguments**
- rs : ResultSet
  
- n : `n` locations to retrieve
  
- kwargs : dimensions to slice across
  

**Returns**

YAXArray[locations, [loc_id, loc_name, rank], scenarios]


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/ranks.jl#L137-L150" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics.trajectory_heatmap-Tuple{YAXArrays.Cubes.YAXArray}' href='#ADRIA.metrics.trajectory_heatmap-Tuple{YAXArrays.Cubes.YAXArray}'><span class="jlbinding">ADRIA.metrics.trajectory_heatmap</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
trajectory_heatmap(data::YAXArray)::HeatMap
```


Estimate heatmap of trajectories from a 2D dataset.

**Arguments**
- data : An N*D matrix where N is time steps and D is the scenario outcome for the given timestep in N
  

**Returns**

OnlineStats.HeatMap


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/temporal.jl#L173-L183" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics.trajectory_heatmap_data-Tuple{YAXArrays.Cubes.YAXArray}' href='#ADRIA.metrics.trajectory_heatmap_data-Tuple{YAXArrays.Cubes.YAXArray}'><span class="jlbinding">ADRIA.metrics.trajectory_heatmap_data</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
trajectory_heatmap_data(data::YAXArray)::Tuple{Vector{Float64},Vector{Float64},Matrix{Int64}}
```


Estimate heatmap of trajectories from a 2D dataset.

**Arguments**
- data : An N*D matrix where N is time steps and D is the scenario outcome for the given timestep in N
  

**Returns**

Tuple of xedges, yedges, and bi-dimensional histogram matrix


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/temporal.jl#L191-L201" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='Base.ndims-Tuple{ADRIA.metrics.Metric}' href='#Base.ndims-Tuple{ADRIA.metrics.Metric}'><span class="jlbinding">Base.ndims</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
ndims(m::Metric)::Int64
```


Infer the number of dimensions for a given outcome/metric.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/utils.jl#L56-L60" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.metrics.Metric-Tuple{Any, Vararg{Any}}' href='#ADRIA.metrics.Metric-Tuple{Any, Vararg{Any}}'><span class="jlbinding">ADRIA.metrics.Metric</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
(f::Metric)(raw, args...; kwargs...)
(f::Metric)(rs::ResultSet, args...; kwargs...)
```


Makes Metric types callable with arbitrary arguments that are passed to associated function.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/metrics.jl#L48-L53" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Performance indicators {#Performance-indicators}
<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.performance.RMSE-Tuple{Any, Any}' href='#ADRIA.performance.RMSE-Tuple{Any, Any}'><span class="jlbinding">ADRIA.performance.RMSE</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Root Mean Square Error


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/performance.jl#L20" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.performance.environmental_diversity-Tuple{Any, Any}' href='#ADRIA.performance.environmental_diversity-Tuple{Any, Any}'><span class="jlbinding">ADRIA.performance.environmental_diversity</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
environmental_diversity(ms, inputs_i)
```


Obtain an indication of environmental factor diversity for a scenario set. Higher values indicate a greater of mix of environmental conditions were experienced between scenarios.

This is referred to as $E$.

**Arguments**
- ms : model spec
  
- inputs_i : inputs used for scenarios of interest
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/performance.jl#L174-L186" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.performance.gmd-Tuple{AbstractVector{<:Real}}' href='#ADRIA.performance.gmd-Tuple{AbstractVector{<:Real}}'><span class="jlbinding">ADRIA.performance.gmd</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
gmd(vals::AbstractVector{<:Real})::Float64
gmd(vals::AbstractMatrix{<:Real})
```


Gini&#39;s Mean Difference.

The absolute mean of all pairwise distances between elements in a given set.

**References**
1. La Haye, R., &amp; Zizler, P. (2019). The Gini mean difference and variance. METRON, 77(1), 43-52. https://doi.org/10.1007/s40300-019-00149-2
  
2. Yitzhaki, S. (2003). Gini&#39;s Mean difference: A superior measure of variability for non-normal   distributions. Metron - International Journal of Statistics, LXI(2), 285-316. https://ideas.repec.org/a/mtn/ancoec/030208.html
  
3. Kashif, M., Aslam, M., Al-Marshadi, A. H., &amp; Jun, C.-H. (2016). Capability Indices for Non-Normal Distribution Using Gini&#39;s Mean Difference as Measure of Variability. IEEE Access, 4, 7322-7330. https://doi.org/10.1109/ACCESS.2016.2620241
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/performance.jl#L32-L56" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.performance.intervention_diversity-Tuple{Any, Any}' href='#ADRIA.performance.intervention_diversity-Tuple{Any, Any}'><span class="jlbinding">ADRIA.performance.intervention_diversity</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
intervention_diversity(ms, inputs_i)
```


Obtain an indication of intervention diversity for a scenario. Higher values indicate a greater of mix of interventions options were applied.

This is referred to as $D$.

**Arguments**
- ms : model spec
  
- inputs_i : inputs used for scenarios of interest
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/performance.jl#L158-L169" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.performance.intervention_effort-Tuple{Any, Any, Any}' href='#ADRIA.performance.intervention_effort-Tuple{Any, Any, Any}'><span class="jlbinding">ADRIA.performance.intervention_effort</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
intervention_effort(ms, inputs_i)
```


Obtain an indication of intervention effort for each scenario and intervention type. This is referred to as $F$.

**Arguments**
- ms : model spec
  
- inputs_i : inputs used for scenarios of interest
  

**Returns**

Matrix of `s * 8`, where `s` is the number of scenarios and columns are: `N_seed_TA`, `N_seed_CA`, `N_seed_CNA`, `N_seed_SM`, `N_seed_LM`, `fogging`, `SRM`, `seed_years`, `shade_years`, `fog_years`


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/performance.jl#L109-L123" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.performance.normalize-Tuple{AbstractArray{<:Real}}' href='#ADRIA.performance.normalize-Tuple{AbstractArray{<:Real}}'><span class="jlbinding">ADRIA.performance.normalize</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
normalize(vals::AbstractArray{<:Real})
```


Normalize values using feature scaling such that values are bound between 0 and 1, where 1 is equivalent to the maximum value found.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/performance.jl#L6-L11" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.performance.probability-Tuple{AbstractArray{<:Real}}' href='#ADRIA.performance.probability-Tuple{AbstractArray{<:Real}}'><span class="jlbinding">ADRIA.performance.probability</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
probability(vals::AbstractArray{<:Real})
```


Calculate probability of individual trajectories, given a scenario ensemble $S$.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/performance.jl#L23-L27" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.performance.temporal_variability-Tuple{AbstractVector{<:Real}}' href='#ADRIA.performance.temporal_variability-Tuple{AbstractVector{<:Real}}'><span class="jlbinding">ADRIA.performance.temporal_variability</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
temporal_variability(x::AbstractVector{<:Real})
temporal_variability(x::AbstractArray{<:Real, 2})
temporal_variability(x::AbstractArray{<:Real}, func_or_data...)
```


The V meta-metric.

As a meta-metric, it can be applied to any combination of metrics (including itself), assuming $x$ is bound between 0 and 1. If this is not the case, consider normalizing values first.

By default (`detrend=true`), variability is assessed on the first differences of $x$ (i.e. step-to-step change), so a smoothly trending series is treated as stable while an erratic one is not. Set `detrend=false` to instead assess variability on the raw values of $x$ (order-independent spread of levels).

**Examples**

```julia
# Apply V to a time series
julia> temporal_variability(rand(50))

# Apply V to an ensemble of time series
julia> x = rand(50, 200)
julia> temporal_variability(x)

# Create and apply a modified V metric to an ensemble of time series.
# Where the argument is an array and not a function, the data is used directly
# and so it is assumed all matrices are of the same size and shape.
julia> temporal_variability(x, temporal_variabilty, temporal_variability(P(x)))
julia> temporal_variability(x, temporal_variabilty, P(x), D(x), E(x))
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/metrics/performance.jl#L66-L97" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Sensitivity
<details class='jldocstring custom-block' open>
<summary><a id='ADRIAanalysis.sensitivity._cramers_v-Tuple{AbstractVector, AbstractVector{Bool}}' href='#ADRIAanalysis.sensitivity._cramers_v-Tuple{AbstractVector, AbstractVector{Bool}}'><span class="jlbinding">ADRIAanalysis.sensitivity._cramers_v</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_cramers_v(col_vals, selection_mask) -> (statistic::Float64, effect_size::Float64)
```


Cramer&#39;s V association strength between a nominal/unordered-categorical feature column `col_vals` and a binary `selection_mask`, via `HypothesisTests.ChisqTest` on the `(n_levels x 2)` contingency table of factor level vs. group membership.

`V = sqrt(chi2 / (n * (min(n_rows, n_cols) - 1)))`, unsigned, in `[0, 1]`. `n_cols` is always 2 (the two `selection_mask` groups); `n_rows` is the number of unique values in `col_vals`. When `col_vals` has fewer than 2 unique values, `min(n_rows, n_cols) - 1` would be zero (division by zero) – mirroring the zero-variance handling in `rsa`, this case emits a `@warn` and returns the sentinel `(0.0, 0.0)` instead.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIAanalysis/src/sensitivity/rsa.jl#L232-L244" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIAanalysis.sensitivity._rank_aligned_delta-Tuple{Vector{Float64}, Vector{Float64}}' href='#ADRIAanalysis.sensitivity._rank_aligned_delta-Tuple{Vector{Float64}, Vector{Float64}}'><span class="jlbinding">ADRIAanalysis.sensitivity._rank_aligned_delta</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_rank_aligned_delta(y_iv, y_cf_raw) -> (Vector{Float64}, Vector{Int})
```


Sort both outcome vectors ascending and compute element-wise delta at matched rank positions.  When lengths differ, the CF empirical quantile function is evaluated at the midpoint probability for each of the `n_iv` rank positions (`p_r = (r - 0.5) / n_iv`) to avoid boundary extrapolation.

Returns `(y_delta, perm)` where `perm` is the sort permutation applied to `y_iv`; apply `fs_iv[perm, :]` to keep the feature matrix aligned with `y_delta`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIAanalysis/src/sensitivity/rsa.jl#L546-L557" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIAanalysis.sensitivity._stratified_rsa-Tuple{Any, DataFrames.DataFrame, AbstractVector, Symbol, Int64}' href='#ADRIAanalysis.sensitivity._stratified_rsa-Tuple{Any, DataFrames.DataFrame, AbstractVector, Symbol, Int64}'><span class="jlbinding">ADRIAanalysis.sensitivity._stratified_rsa</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
_stratified_rsa(compute_si, fs, labels, strat_col, n_strata) -> DataFrame
```


Shared stratification/aggregation core for `stratified_rsa`. Bins `fs` into `n_strata` equal-frequency quantile groups on `strat_col`, drops DHW stat columns and zero-variance columns within each stratum, and calls `compute_si(X_s, labels_s)` per stratum – `labels` is either the scalar outcome vector `y` or a `selection_mask`, sliced to the stratum&#39;s rows. `compute_si` returns a per-stratum `rsa` result DataFrame, or `nothing` to skip the stratum (after emitting its own `@warn`).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIAanalysis/src/sensitivity/rsa.jl#L462-L472" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIAanalysis.sensitivity.convergence-Tuple{DataFrames.DataFrame, YAXArrays.Cubes.YAXArray, Vector{Symbol}}' href='#ADRIAanalysis.sensitivity.convergence-Tuple{DataFrames.DataFrame, YAXArrays.Cubes.YAXArray, Vector{Symbol}}'><span class="jlbinding">ADRIAanalysis.sensitivity.convergence</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
convergence(X::DataFrame, y::YAXArray, target_factors::Vector{Symbol}; n_steps::Int64=10)::YAXArray
convergence(rs::ResultSet, X::DataFrame, y::YAXArray, components::Vector{String}; n_steps::Int64=10)::YAXArray
```


Calculates the PAWN sensitivity index for an increasing number of scenarios where the maximum is the total number of scenarios in scens. Number of scenario subsets determined by N_steps. Can be calculated for individual factors or aggregated over factors for specified model components.

**Arguments**
- `rs` : Result set (only needed if aggregating over model components).
  
- `X` : Model inputs
  
- `y` : Model outputs
  
- `target_factors` : Names of target factors represented by columns in `X`.
  
- `components` : Names of model components to aggregate over (e.g. [:Intervention, :Criteria]).
  
- `n_steps` : Number of steps to cut the total number of scenarios into.
  

**Returns**

YAXArray, of min, lower bound, mean, median, upper bound, max, std, and cv summary statistics for an increasing number of scenarios.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIAanalysis/src/sensitivity/sensitivity.jl#L212-L232" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIAanalysis.sensitivity.counterfactual_delta-Tuple{ADRIA.ResultSet, ADRIA.ResultSet, Any}' href='#ADRIAanalysis.sensitivity.counterfactual_delta-Tuple{ADRIA.ResultSet, ADRIA.ResultSet, Any}'><span class="jlbinding">ADRIAanalysis.sensitivity.counterfactual_delta</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
counterfactual_delta(
    rs_intervention::ResultSet,
    rs_counterfactual::ResultSet,
    metric_fn;
    bootstrap_n::Int=0
) -> NamedTuple
```


Compute per-scenario intervention lift (delta) and return the intervention feature matrix with DHW columns dropped, ready for `rsa`.

Both `rs_intervention` and `rs_counterfactual` are assumed to be independently Sobol&#39;-sampled over the same parameter bounds.

Two estimators are available via `delta_method`:

**`:rank_aligned`** (default): Both outcome vectors are sorted ascending; the CF vector is interpolated to `n_iv` quantile positions if lengths differ. The lift for rank position `r` is:

```julia
y_delta[r] = y_iv^(r) - y_cf^(r)
```


`fs_intervention` rows are permuted to the same sort order so that `y_delta[r]` and `fs_intervention[r, :]` always correspond to the same intervention scenario.  RSA on this delta asks: &quot;which factors characterise scenarios that beat the same-rank counterfactual outcome?&quot;

**`:mean_difference`**: The lift for intervention scenario `i` is:

```julia
y_delta[i] = y_iv[i] - mean(y_cf)
```


Note: because `mean(y_cf)` is a constant offset, RSA on this delta is mathematically equivalent to RSA on `y_iv` directly. `fs_intervention` row order is unchanged.

**Arguments**
- `rs_intervention`   : ResultSet from an intervention run.
  
- `rs_counterfactual` : ResultSet from a no-intervention (counterfactual) run.
  
- `metric_fn`         : Function `rs -> AbstractVector{<:Real}` mapping a                       ResultSet to a scalar outcome per scenario.
  
- `delta_method`      : `:rank_aligned` (default) or `:mean_difference`; see above.
  
- `bootstrap_n`       : Number of bootstrap resamples for a 95% CI on                       `mean(y_delta)`.  0 (default) skips bootstrapping.
  
- `fs`                : Precomputed `feature_set(rs_intervention)` to reuse, for                       callers that invoke `counterfactual_delta` repeatedly                       against the same `rs_intervention` with different                       `metric_fn`/location restrictions (e.g. per-period or                       per-location-set repeats). `feature_set` is a pure                       function of `rs_intervention` alone – independent of                       `metric_fn` – so recomputing it on every call                       (deployment-log summaries, DHW stats) is pure overhead                       once a caller needs more than one delta from the same                       `rs_intervention`. Must have exactly `length(metric_fn(                       rs_intervention))` rows, in `rs_intervention`&#39;s original                       scenario order (i.e. `feature_set(rs_intervention)`                       itself, unpermuted). `nothing` (default) computes it                       internally, as before.
  

**Returns**

`NamedTuple` with fields:
- `y_delta`        : `Vector{Float64}` – per-scenario lift (length = n intervention                    scenarios).  With `:rank_aligned`, sorted ascending by outcome.
  
- `fs_intervention`: `DataFrame` – `feature_set(rs_intervention)` with DHW stat columns                    (`:dhw_mean`, `:dhw_stdev`, `:dhw_complexity`) removed.                    With `:rank_aligned`, rows are permuted to match `y_delta` order.
  
- `iv_perm`        : `Vector{Int}` – permutation applied to the intervention rows,                    i.e. `fs_intervention`/`y_delta` row `i` corresponds to the `i`-th                    original intervention scenario (in `metric_fn(rs_intervention)`                    order) at index `iv_perm[i]`. Identity for `:mean_difference`.                    Callers holding other per-original-intervention-row arrays (e.g.                    a `guided`/`unguided` mask) must index them with `iv_perm` before                    lining them up against `fs_intervention`/`y_delta`.
  
- `bootstrap_ci`   : `Union{Nothing, Tuple{Float64,Float64}}` – bootstrapped 95% CI on                    `mean(y_delta)`, or `nothing` if `bootstrap_n == 0`
  

**Edge cases**
- If `all(y_delta .== 0)`, a `@warn` is emitted (ATE is zero; RSA will return sentinels).
  
- If either ResultSet has zero scenarios, throws `ArgumentError`.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIAanalysis/src/sensitivity/rsa.jl#L576-L654" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIAanalysis.sensitivity.counterfactual_delta-Tuple{ADRIA.ResultSet, AbstractVector{Bool}, Any}' href='#ADRIAanalysis.sensitivity.counterfactual_delta-Tuple{ADRIA.ResultSet, AbstractVector{Bool}, Any}'><span class="jlbinding">ADRIAanalysis.sensitivity.counterfactual_delta</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
counterfactual_delta(
    rs::ResultSet,
    cf_mask::AbstractVector{Bool},
    metric_fn;
    bootstrap_n::Int=0
) -> NamedTuple
```


Variant of `counterfactual_delta` for result sets where intervention and counterfactual scenarios are stored together (e.g. `rs.inputs.guided .== -1` marks the counterfactual group).

`metric_fn` is called once with the full `rs`; its output (one value per scenario) is then split by `cf_mask`.  Intervention scenarios are all rows where `cf_mask` is `false`.

All other behaviour (`delta_method`, DHW column removal, optional bootstrap CI, the `iv_perm` return field) is identical to the two-ResultSet overload. Here `iv_perm` indexes into `findall(.!cf_mask)` order, i.e. the intervention rows of `rs` in their original order.

`fs` (optional): precomputed `feature_set(rs)` to reuse across repeated calls against the same `rs` with different `metric_fn`/location restrictions – see the two-ResultSet overload&#39;s docstring for the rationale. Must have `length(cf_mask)` rows in `rs`&#39;s original scenario order (i.e. `feature_set(rs)` itself, unpermuted and NOT pre-filtered to `iv_mask`). `nothing` (default) computes it internally, as before.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIAanalysis/src/sensitivity/rsa.jl#L719-L746" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIAanalysis.sensitivity.ks_statistic-Tuple{HypothesisTests.ApproximateKSTest}' href='#ADRIAanalysis.sensitivity.ks_statistic-Tuple{HypothesisTests.ApproximateKSTest}'><span class="jlbinding">ADRIAanalysis.sensitivity.ks_statistic</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
ks_statistic(ks)
```


Calculate the Kolmogorov-Smirnov test statistic.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIAanalysis/src/sensitivity/sensitivity.jl#L23-L27" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIAanalysis.sensitivity.pawn-Tuple{AbstractMatrix{<:Real}, AbstractVector{<:Real}, Vector{String}}' href='#ADRIAanalysis.sensitivity.pawn-Tuple{AbstractMatrix{<:Real}, AbstractVector{<:Real}, Vector{String}}'><span class="jlbinding">ADRIAanalysis.sensitivity.pawn</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
pawn(rs::ResultSet, y::Union{NamedDimsArray,AbstractVector{<:Real}}; S::Int64=10)::NamedDimsArray
pawn(X::AbstractMatrix{<:Real}, y::AbstractVector{<:Real}, factor_names::Vector{String}; S::Int64=10)::NamedDimsArray
pawn(X::DataFrame, y::AbstractVector{<:Real}; S::Int64=10)::NamedDimsArray
pawn(X::NamedDimsArray, y::Union{NamedDimsArray,AbstractVector{<:Real}}; S::Int64=10)::NamedDimsArray
pawn(X::Union{DataFrame,AbstractMatrix{<:Real}}, y::AbstractMatrix{<:Real}; S::Int64=10)::NamedDimsArray
```


Calculates the PAWN sensitivity index.

The PAWN method (by Pianosi and Wagener) is a moment-independent approach to Global Sensitivity Analysis. Outputs are characterized by their Cumulative Distribution Function (CDF), quantifying the variation in the output distribution after conditioning an input over &quot;slices&quot; ($S$) - the conditioning intervals. If both distributions coincide at all slices (i.e., the distributions are similar or identical), then the factor is deemed non-influential.

This implementation applies the Kolmogorov-Smirnov test as the distance measure and returns summary statistics (min, lower bound, mean, median, upper bound, max, std, and cv) over the slices.

**Arguments**
- `rs` : ResultSet
  
- `X` : Model inputs
  
- `y` : Model outputs
  
- `factor_names` : Names of each factor represented by columns in `X`
  
- `S` : Number of slides (default: 10)
  

**Returns**

YAXArray, of min, mean, lower bound, median, upper bound, max, std, and cv summary statistics.

**Examples**

```julia
dom = ADRIA.load_domain("example_domain", "<RCP>")
scens = ADRIA.sample(dom, 128)
rs = ADRIA.run_scenarios(dom, scens, "45")

# Get mean coral cover over time and locations
μ_tac = mean(ADRIA.metrics.scenario_total_cover(rs), dims=:timesteps)

ADRIAanalysis.sensitivity.pawn(rs, μ_tac)
```


**References**
1. Pianosi, F., Wagener, T., 2018. Distribution-based sensitivity analysis from a generic input-output sample. Environmental Modelling &amp; Software 108, 197-207. https://doi.org/10.1016/j.envsoft.2018.07.019
  
2. Baroni, G., Francke, T., 2020. GSA-cvd Combining variance- and distribution-based global sensitivity analysis https://github.com/baronig/GSA-cvd
  
3. Puy, A., Lo Piano, S., &amp; Saltelli, A. 2020. A sensitivity analysis of the PAWN sensitivity index. Environmental Modelling &amp; Software, 127, 104679. https://doi.org/10.1016/j.envsoft.2020.104679
  
4. https://github.com/SAFEtoolbox/Miscellaneous/blob/main/Review_of_Puy_2020.pdf
  

**Extended help**

Pianosi and Wagener have made public their review responding to a critique of their method by Puy et al., (2020). A key criticism by Puy et al. was that the PAWN method is sensitive to its tuning parameters and thus may produce biased results. The tuning parameters referred to are the number of samples ($N$) and the number of conditioning points - $n$ in Puy et al., but denoted as $S$ here.

Puy et al., found that the ratio of $N$ (number of samples) to $S$ has to be sufficiently high ($N/S > 80$) to avoid biased results. Pianosi and Wagener point out this requirement is not particularly difficult to meet. Using the recommended value ($S := 10$), a sample of 1024 runs (small for purposes of Global Sensitivity Analysis) meets this requirement ($1024/10 = 102.4$). Additionally, lower values of $N/S$ is more an indication of faulty experimental design moreso than any deficiency of the PAWN method.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIAanalysis/src/sensitivity/sensitivity.jl#L34-L106" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIAanalysis.sensitivity.quantile_strata_edges-Tuple{AbstractVector{<:Real}, Int64}' href='#ADRIAanalysis.sensitivity.quantile_strata_edges-Tuple{AbstractVector{<:Real}, Int64}'><span class="jlbinding">ADRIAanalysis.sensitivity.quantile_strata_edges</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
quantile_strata_edges(vals::AbstractVector{<:Real}, n_strata::Int) -> Vector{Float64}
```


Equal-frequency quantile bin edges for `vals` (length `n_strata + 1`), with the final edge nudged up by `nextfloat` so the maximum value falls inside the last bin: stratum `s` is `edges[s] <= v < edges[s+1]`.

Used internally by `stratified_rsa` to bin `strat_col`; exposed so callers can derive a matching bin assignment for a _different_ vector against the same edges (e.g. binning counterfactual scenarios&#39; DHW values into the same strata as the intervention scenarios used to derive `edges`, as `stratified_cf_mask` does).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIAanalysis/src/sensitivity/rsa.jl#L275-L287" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIAanalysis.sensitivity.rsa-Tuple{DataFrames.DataFrame, AbstractVector{<:Real}}' href='#ADRIAanalysis.sensitivity.rsa-Tuple{DataFrames.DataFrame, AbstractVector{<:Real}}'><span class="jlbinding">ADRIAanalysis.sensitivity.rsa</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
rsa(X::DataFrame, y::AbstractVector{<:Real};
    top_proportion::Float64=0.9, method::Symbol=:auto) -> DataFrame
rsa(X::DataFrame, selection_mask::Union{BitVector,AbstractVector{Bool}};
    method::Symbol=:auto) -> DataFrame
```


Rank-based Regional Sensitivity Analysis: score each input factor by how well it discriminates between high- and low-outcome scenario groups.

rsa is a rank-based, scenario-conditioned sensitivity method – it answers &quot;which factors&#39; distributions differ most between high- and low-outcome scenarios?&quot; This complements PAWN/KS-based methods (which measure how the full output distribution shifts across the input space) rather than replacing them. Both are forms of sensitivity analysis; they answer related but distinct questions.

**Test dispatch per column**

Two tests are available:
- **Mann-Whitney U** (`HypothesisTests.MannWhitneyUTest`): appropriate for continuous, ordered-categorical, and ordered-discrete factors, where the two groups&#39; _ranks_ are meaningfully comparable.
  
- **Cramer&#39;s V** (`HypothesisTests.ChisqTest` on a factor-level x group contingency table): appropriate for unordered/nominal categorical factors, where there is no meaningful notion of &quot;rank&quot;.
  

With `method=:auto` (the default), each column of `X` is dispatched individually based on its `DataFrames.colmetadata(X, col, "ptype", "continuous")` value: columns tagged `"unordered categorical"` use Cramer&#39;s V; everything else (`"continuous"`, `"ordered categorical"`, `"ordered discrete"`, `"discrete"`, or no metadata at all) uses Mann-Whitney. Passing `method=:mann_whitney` or `method=:cramers_v` overrides this and forces every column through that single test regardless of its `ptype` tag.

**Backwards-compatible fallback**: if `X` has no `colmetadata` attached to _any_ column (e.g. a plain, untagged `DataFrame`) and at least one column&#39;s `eltype` looks non-numeric (a heuristic proxy for &quot;might actually be nominal but untagged&quot;), a single `@warn` is emitted once before the per-feature loop, noting that no ptype metadata was found and Mann-Whitney is being applied indiscriminately – which is statistically invalid for nominal categorical factors. This is informational only: every untagged column still falls back to `"continuous"`/Mann-Whitney, matching pre-existing behaviour.

**Primary dispatch: scalar outcomes**
- `X`              : Feature matrix (DataFrame, columns = factors)
  
- `y`              : Scalar outcome per scenario; scenarios above the `top_proportion`                    quantile form the selected group
  
- `top_proportion` : Quantile threshold for the high-outcome group (default: 0.9)
  
- `method`         : `:auto` (default, per-column ptype dispatch), `:mann_whitney`, or                    `:cramers_v` (the latter two force uniform application to all columns)
  

**Escape-hatch dispatch: pre-computed mask**
- `X`              : Feature matrix (DataFrame, columns = factors)
  
- `selection_mask` : Boolean mask (true = selected/&quot;high-outcome&quot; group)
  
- `method`         : Ranking method, as above
  

**Returns**

DataFrame sorted descending by `prob_superiority`:
- feature          (Symbol)  : factor name
  
- test             (Symbol)  : `:mann_whitney` or `:cramers_v` – which test produced                              this row
  
- statistic        (Float64) : raw Mann-Whitney U, or the chi-squared statistic for                              Cramer&#39;s V rows
  
- prob_superiority (Float64) : U / (n1 * n2), in [0, 1], for Mann-Whitney rows. Cramer&#39;s                              V has no equivalent &quot;probability of superiority&quot; concept                              (it is an unsigned association strength, not a rank                              comparison), so this is `NaN` for `:cramers_v` rows –                              following the same &quot;not meaningful for this row&quot; convention                              used for the zero-variance sentinel below, except NaN                              rather than 0.5 since there is no neutral value to report.
  
- effect_size      (Float64) : 1 - 2*U / (n1 * n2), in `[-1, 1]` (signed), for Mann-Whitney                              rows; Cramer&#39;s V itself, in `[0, 1]` (unsigned), for                              `:cramers_v`rows. These two are NOT on a comparable scale                              -- distinguish them via the`test`column, do not compare`effect_size` values across test types directly.
  

Because Mann-Whitney and Cramer&#39;s V effect sizes are not comparable, a result table containing BOTH test types does not produce a single meaningfully-ranked-together ordering: `sort!(result, :prob_superiority; rev=true)` will place all `:cramers_v` rows (prob_superiority = NaN) according to Julia&#39;s NaN sort ordering, not by association strength. Callers wanting a proper within-statistic ranking should filter by `test` first (e.g. `filter(:test =&gt; ==(:mann_whitney), result)`).

HypothesisTests.MannWhitneyUTest applies a normal approximation with tie correction. A @warn is emitted when more than 20% of values in a feature column are tied, as the effect_size formula becomes less reliable in that case.

**Known confound**: for factors like `mcda_method` (sentinel-zeroed/inapplicable when `guided <= 0`), Cramer&#39;s V will partly re-detect &quot;is guided active at all&quot;, a signal `guided`&#39;s own Mann-Whitney effect size already captures separately. This is a known, expected confound – it is documented here rather than engineered around.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIAanalysis/src/sensitivity/rsa.jl#L1-L88" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIAanalysis.sensitivity.stratified_cf_mask-Tuple{AbstractVector{<:Real}, AbstractVector{<:Real}, AbstractVector{<:Real}, AbstractVector{<:Real}}' href='#ADRIAanalysis.sensitivity.stratified_cf_mask-Tuple{AbstractVector{<:Real}, AbstractVector{<:Real}, AbstractVector{<:Real}, AbstractVector{<:Real}}'><span class="jlbinding">ADRIAanalysis.sensitivity.stratified_cf_mask</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
stratified_cf_mask(y_iv::AbstractVector{<:Real}, dhw_iv::AbstractVector{<:Real},
                    y_cf::AbstractVector{<:Real}, dhw_cf::AbstractVector{<:Real};
                    n_strata::Int=4, quantile_level::Float64=0.8,
                    min_stratum_n::Int=10) -> BitVector
```


Behavioural mask for intervention scenarios, with the &quot;behavioural&quot; bar set _locally per DHW stratum_: scenario `i` is behavioural iff `y_iv[i]` beats the `quantile_level` quantile of counterfactual outcomes (`y_cf`) drawn from the SAME DHW stratum as scenario `i`.

Bin edges are the equal-frequency quantile edges of `dhw_iv` (see `quantile_strata_edges`); `dhw_cf` is binned into those same edges so each stratum&#39;s threshold is computed from counterfactual scenarios under comparable climate conditions.

This is deliberately local rather than global: a single global `quantile(y_cf, quantile_level)` threshold, especially at a high `quantile_level`, is disproportionately clearable only by scenarios drawn under mild DHW – passed to `stratified_rsa`, that collapses the &quot;which factors matter within this DHW band&quot; question back into &quot;which DHW band is mild&quot;, which stratification exists to avoid. Since DHW is fixed within a stratum here, &quot;beat the local top `1 - quantile_level` of counterfactual outcomes&quot; isolates intervention-parameter effects from climate severity instead.

Strata with fewer than `min_stratum_n` counterfactual scenarios fall back to the GLOBAL `quantile_level` quantile of `y_cf`, with a `@warn`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIAanalysis/src/sensitivity/rsa.jl#L302-L329" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIAanalysis.sensitivity.stratified_rsa-Tuple{DataFrames.DataFrame, AbstractVector{<:Real}}' href='#ADRIAanalysis.sensitivity.stratified_rsa-Tuple{DataFrames.DataFrame, AbstractVector{<:Real}}'><span class="jlbinding">ADRIAanalysis.sensitivity.stratified_rsa</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
stratified_rsa(fs::DataFrame, y::AbstractVector{<:Real};
               strat_col::Symbol=:dhw_mean,
               n_strata::Int=4,
               top_proportion::Float64=0.9) -> DataFrame
stratified_rsa(fs::DataFrame, selection_mask::Union{BitVector,AbstractVector{Bool}};
               strat_col::Symbol=:dhw_mean,
               n_strata::Int=4) -> DataFrame
```


Run `rsa` independently within each DHW quantile stratum and return a long-format DataFrame summarising factor importance across strata.

Scenarios are binned into `n_strata` equal-frequency quantile groups on `strat_col` (default `:dhw_mean`).  DHW stat columns (`dhw_mean`, `dhw_stdev`, `dhw_complexity`) are dropped from the feature matrix within each stratum before calling `rsa`, so they do not confound the within-stratum ranking.

**Primary dispatch: scalar outcomes**
- `fs`             : Feature matrix as returned by `feature_set(rs)`, which                    must contain `strat_col` as a column.
  
- `y`              : Scalar outcome per scenario (same row order as `fs`).
  
- `strat_col`      : Column of `fs` used to form strata (default `:dhw_mean`).
  
- `n_strata`       : Number of equal-frequency quantile bins (default 4).
  
- `top_proportion` : Passed through to `rsa`; quantile threshold for the                    high-outcome group within each stratum (default 0.9).
  

**Escape-hatch dispatch: pre-computed mask**
- `fs`             : Feature matrix, as above.
  
- `selection_mask` : Boolean mask (true = selected/&quot;behavioural&quot; group), same                    row order as `fs`. Sliced per stratum and passed straight                    to `rsa`&#39;s mask dispatch – the caller decides how the                    behavioural split is defined (e.g. against a counterfactual                    threshold) rather than `stratified_rsa` re-deriving a                    per-stratum quantile from `y` itself.
  
- `strat_col`, `n_strata` : As above.
  

**Returns**

Long-format `DataFrame` with columns:
- `feature`          (Symbol)  : factor name
  
- `stratum`          (Int)     : stratum index 1..n_strata (low-&gt;high DHW)
  
- `prob_superiority` (Float64) : Mann-Whitney P(X1 &gt; X2) for this factor in this stratum
  
- `effect_size`      (Float64) : 1 - 2U/(n1*n2)
  
- `mean_importance`  (Float64) : mean of abs(prob_superiority - 0.5) across all strata                                 (same value repeated per row). Measures average deviation                                 from neutrality; a factor with reversed importance across                                 strata (e.g. [0.9, 0.1]) scores the same as a consistently                                 important one ([0.9, 0.9]), both higher than a neutral                                 factor ([0.5, 0.5]).
  

Strata that contain fewer than 10 scenarios are skipped with a `@warn`. For the scalar-outcome dispatch, strata where all outcome values are identical are also skipped; for the mask dispatch, strata where the sliced mask is all-true or all-false are skipped instead. Skipped strata are absent from the returned DataFrame.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIAanalysis/src/sensitivity/rsa.jl#L374-L429" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIAanalysis.sensitivity.tsa-Tuple{DataFrames.DataFrame, AbstractMatrix{<:Real}}' href='#ADRIAanalysis.sensitivity.tsa-Tuple{DataFrames.DataFrame, AbstractMatrix{<:Real}}'><span class="jlbinding">ADRIAanalysis.sensitivity.tsa</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
tsa(X::DataFrame, y::AbstractMatrix)::YAXArray
```


Perform Temporal (or time-varying) Sensitivity Analysis using the PAWN sensitivity index.

The sensitivity index value for time $t$ is inclusive of all time steps prior to $t$. Alternate approaches use a moving window, or only data for time $t$.

**Examples**

```julia
rs = ADRIA.load_results("a ResultSet of interest")

# Get scenario outcomes over time (shape: `time × scenarios`)
y_tac = ADRIA.metrics.scenario_total_cover(rs)

# Calculate sensitivity of outcome to factors for each time step
ADRIAanalysis.sensitivity.tsa(rs.inputs, y_tac)
```


**Arguments**
- `X` : Scenario specification
  
- `y` : scenario outcomes over time
  

**Returns**

YAXArray, of shape $D$ × 6 × $T$, where
- 
  $$D$$
  is the number of dimensions/factors
  
- 6 corresponds to the min, mean, median, max, std, and cv of the PAWN indices
  
- 
  $$T$$
  is the number of time steps
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIAanalysis/src/sensitivity/sensitivity.jl#L297-L325" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## General API {#General-API}
<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.bin_edges-Tuple{}' href='#ADRIA.bin_edges-Tuple{}'><span class="jlbinding">ADRIA.bin_edges</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
bin_edges()
```


Helper function defining coral colony diameter bin edges. The values are converted from `cm` to the desired unit. The default target unit is `cm`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/ecosystem/corals/coral_factors.jl#L139-L144" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.coral_spec-Tuple{}' href='#ADRIA.coral_spec-Tuple{}'><span class="jlbinding">ADRIA.coral_spec</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
coral_spec()::NamedTuple
```


Return the coral parameter specification for ADRIA. Includes &quot;vital&quot; bio/ecological parameters for each coral taxa, functional group, and size class.

Results are cached: the value is computed once per session and persisted to disk at `~/.julia/cache/ADRIA/coral_spec-julia<X.Y>-adria<V>.cache`. The cache filename encodes the Julia minor version and ADRIA version, so stale files from other versions are never loaded. Within a file, the cache is also keyed to a hash of `Corals.jl`, so it is invalidated automatically whenever the source file changes.

**Returns**

A `NamedTuple` with fields:
- `taxa_names`  : names of functional groups
  
- `param_names` : names of perturbable parameters
  
- `params`      : `DataFrame` of parameter values for each taxa/size-class combination
  

**Developer notes**

To force recomputation within a live session (e.g. after reloading with Revise):

```julia
ADRIA.invalidate_coral_spec_cache!()
```


If struct fields changed, also call `ADRIA.create_coral_struct()` afterwards.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/ecosystem/corals/Corals.jl#L209-L234" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.create_coral_instance' href='#ADRIA.create_coral_instance'><span class="jlbinding">ADRIA.create_coral_instance</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
create_coral_instance(bounds=(0.9, 1.1); overrides=Dict())
```


Construct a `Coral` instance with calibrated field values without redefining the struct. Use this instead of `create_coral_struct` when only an instance (not a struct redefinition) is needed.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/ecosystem/corals/Corals.jl#L449-L455" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.create_coral_struct' href='#ADRIA.create_coral_struct'><span class="jlbinding">ADRIA.create_coral_struct</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
create_coral_struct(bounds=(0.9, 1.1))
```


Generates Coral struct using the default parameter spec.

**Example**

```julia
# Define coral struct with auto-generated parameter ranges
# (default in ADRIA is ± 10%, triangular distribution with peak at 0.5)
create_coral_struct()
coral = Coral()

# Recreate coral spec ± 50% from nominal values
create_coral_struct((0.5, 1.5))
coral = Coral()
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/ecosystem/corals/Corals.jl#L424-L440" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.env_stats-Tuple{ADRIA.ResultSet, String, String}' href='#ADRIA.env_stats-Tuple{ADRIA.ResultSet, String, String}'><span class="jlbinding">ADRIA.env_stats</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
env_stats(rs::ResultSet, s_name::String, rcp::String)
env_stats(rs::ResultSet, s_name::String, rcp::String, scenario::Int)
env_stats(rs::ResultSet, s_name::String, stat::String, rcp::String, scenario::Int)
```


Extract statistics for a given environmental layer (&quot;DHW&quot; or &quot;wave&quot;)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/io/ResultSet.jl#L474-L480" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.loc_area-Tuple{ADRIA.ResultSet}' href='#ADRIA.loc_area-Tuple{ADRIA.ResultSet}'><span class="jlbinding">ADRIA.loc_area</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
loc_area(rs::ResultSet)::Vector{Float64}
```


Extract vector of a location&#39;s total area in its areal unit (m², km², etc).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/io/ResultSet.jl#L592-L596" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.loc_area-Tuple{Domain}' href='#ADRIA.loc_area-Tuple{Domain}'><span class="jlbinding">ADRIA.loc_area</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
loc_area(domain::Domain)::Vector{Float64}
```


Get location area for the given domain.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/Domain.jl#L314-L318" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.loc_coral_cover-Tuple{AbstractArray{Float64, 3}}' href='#ADRIA.loc_coral_cover-Tuple{AbstractArray{Float64, 3}}'><span class="jlbinding">ADRIA.loc_coral_cover</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
loc_coral_cover(C_cover_t::Array{Float64,3})::Vector{Float64}
```


Sum coral cover across all functional groups and size classes of a single timestep for each location.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/Domain.jl#L361-L365" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.loc_k_area-Tuple{ADRIA.ResultSet}' href='#ADRIA.loc_k_area-Tuple{ADRIA.ResultSet}'><span class="jlbinding">ADRIA.loc_k_area</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
loc_k_area(rs::ResultSet)::Vector{Float64}
```


Extract vector of a location&#39;s coral carrying capacity in terms of absolute area.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/io/ResultSet.jl#L573-L577" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.loc_k_area-Tuple{Domain}' href='#ADRIA.loc_k_area-Tuple{Domain}'><span class="jlbinding">ADRIA.loc_k_area</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
loc_k_area(domain::Domain)::Vector{Float64}
```


Get maximum coral cover area for the given domain in absolute area.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/Domain.jl#L325-L329" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.loc_recruits_cover-Tuple{Matrix{Float64}}' href='#ADRIA.loc_recruits_cover-Tuple{Matrix{Float64}}'><span class="jlbinding">ADRIA.loc_recruits_cover</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
loc_recruits_cover(recruits::Matrix{Float64})::Vector{Float64}
```


Absolute cover of recruits on each location.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/Domain.jl#L370-L374" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.run_scenario-Tuple{Domain, Int64, Union{DataFrames.DataFrameRow, AbstractVector}, Vector{Vector{CoralBlox.FunctionalGroup}}, NamedTuple}' href='#ADRIA.run_scenario-Tuple{Domain, Int64, Union{DataFrames.DataFrameRow, AbstractVector}, Vector{Vector{CoralBlox.FunctionalGroup}}, NamedTuple}'><span class="jlbinding">ADRIA.run_scenario</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
run_scenario(domain::Domain, idx::Int64, scenario::Union{AbstractVector,DataFrameRow}, functional_groups::Vector{Vector{FunctionalGroup}}, data_store::NamedTuple)::Nothing
run_scenario(domain::Domain, scenario::Union{AbstractVector,DataFrameRow})::NamedTuple
run_scenario(domain::Domain, scenario::Union{AbstractVector,DataFrameRow}, RCP::String)::NamedTuple
```


Run individual scenarios for a given domain, saving results to a Zarr data store. Results are stored in Zarr format at a pre-configured location. Sets up a new `cache` if not provided.

**Arguments**
- `domain` : Simulation domain (may be modified via `switch_RCPs!`).
  
- `idx` : Scenario index, to store results into `data_store`.
  
- `scenario` : Parameter row describing the scenario.
  
- `functional_groups` : Preallocated functional group buffers.
  
- `data_store` : Pre-opened store with arrays to write results into.
  

**Returns**

Nothing


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/scenario.jl#L570-L588" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.select-Tuple{ADRIA.ResultSet, String}' href='#ADRIA.select-Tuple{ADRIA.ResultSet, String}'><span class="jlbinding">ADRIA.select</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
select(r::ResultSet, op::String)
```


Hacky scenario filtering - to be replaced with more robust approach.

Only supports filtering by single attribute. Should be expanded to support filtering metric results too.

**Examples**

```julia
select(result, "guided .> 0.0")

# Above expands to:
# result.inputs.guided .> 0.0
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/io/ResultSet.jl#L521-L536" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.timesteps-Tuple{ADRIA.ResultSet}' href='#ADRIA.timesteps-Tuple{ADRIA.ResultSet}'><span class="jlbinding">ADRIA.timesteps</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
timesteps(rs::ResultSet)
```


Retrieve the time steps represented in the result set.

**Arguments**
- `rs` : ResultSet
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/io/ResultSet.jl#L548-L555" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.timesteps-Tuple{Domain}' href='#ADRIA.timesteps-Tuple{Domain}'><span class="jlbinding">ADRIA.timesteps</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Extract the time steps represented in the data package.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/Domain.jl#L388" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.timesteps-Tuple{YAXArrays.Cubes.YAXArray}' href='#ADRIA.timesteps-Tuple{YAXArrays.Cubes.YAXArray}'><span class="jlbinding">ADRIA.timesteps</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
timesteps(outcomes::YAXArray)::Vector{Int64}
```


Extract time step labels from a YAXArray. Returns an empty `Vector{Int64}` if the array has no `:timesteps` dimension.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/io/ResultSet.jl#L560-L565" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.ADRIADomain' href='#ADRIA.ADRIADomain'><span class="jlbinding">ADRIA.ADRIADomain</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
ADRIADomain{Σ,M,I,D,X,Y,Z}
```


Core ADRIA domain. Represents study area.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/ExtInterface/ADRIA/Domain.jl#L11-L15" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.Domain-Tuple{String, String, String, Vector, Vararg{String, 10}}' href='#ADRIA.Domain-Tuple{String, String, String, Vector, Vararg{String, 10}}'><span class="jlbinding">ADRIA.Domain</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
Domain(name::String, dpkg_path::String, rcp::String, timeframe::Vector, location_data_fn::String, location_id_col::String, cluster_id_col::String, k_area_col::String, area_col::String, init_coral_fn::String, conn_path::String, dhw_fn::String, wave_fn::String, cyclone_mortality_fn::String; calib_params_fn::String="")::ADRIADomain
```


Convenience constructor for Domain.

**Arguments**
- `name` : Name of domain
  
- `dpkg_path` : location of data package
  
- `rcp` : RCP scenario represented
  
- `timeframe` : Time steps represented
  
- `location_data_fn` : File name of spatial data used
  
- `location_id_col` : Column holding name of reef the location is associated with (non-unique)
  
- `cluster_id_col` : Column holding unique cluster names/ids
  
- `k_area_col` : Column holding habitable area proportion
  
- `area_col` : Column holding location area
  
- `init_coral_fn` : Name of file holding initial coral cover values
  
- `conn_path` : Path to directory holding connectivity data
  
- `dhw_fn` : Filename of DHW data cube in use
  
- `wave_fn` : Filename of wave data cube
  
- `cyclone_mortality_fn` : Filename of cyclone mortality data cube
  
- `calib_params_fn` : path to a CoralBlox calibration NetCDF. If empty or missing, ADRIA
  

default coral and growth acceleration parameters are used.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/ExtInterface/ADRIA/Domain.jl#L167-L189" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.Domain-Union{Tuple{T}, Tuple{String, String, ADRIA.EnvLayer, YAXArrays.Cubes.YAXArray{T, N, A} where {N, A<:AbstractArray{T, N}}, DataFrames.DataFrame, String, String, YAXArrays.Cubes.YAXArray, ADRIA.CoralGrowth, Vector{String}, Vector{String}, YAXArrays.Cubes.YAXArray, YAXArrays.Cubes.YAXArray, YAXArrays.Cubes.YAXArray}} where T<:Union{Float32, Float64}' href='#ADRIA.Domain-Union{Tuple{T}, Tuple{String, String, ADRIA.EnvLayer, YAXArrays.Cubes.YAXArray{T, N, A} where {N, A<:AbstractArray{T, N}}, DataFrames.DataFrame, String, String, YAXArrays.Cubes.YAXArray, ADRIA.CoralGrowth, Vector{String}, Vector{String}, YAXArrays.Cubes.YAXArray, YAXArrays.Cubes.YAXArray, YAXArrays.Cubes.YAXArray}} where T<:Union{Float32, Float64}'><span class="jlbinding">ADRIA.Domain</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Barrier function to create Domain struct without specifying Intervention, Criteria, Coral or SimConstant parameters, which are constructed internally.

Coral and growth acceleration parameters are taken from `calib_params_fn` when given, and default to ADRIA values otherwise.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/ExtInterface/ADRIA/Domain.jl#L74-L80" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='ADRIA.SimConstants' href='#ADRIA.SimConstants'><span class="jlbinding">ADRIA.SimConstants</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
SimConstants
```


Struct of simulation constants for ADRIA

**References**
1. Lough, J. M., Anderson, K. D., &amp; Hughes, T. P. (2018). Increasing thermal stress for tropical coral reefs: 1871-2017. Scientific Reports, 8(1), 6079. https://doi.org/10.1038/s41598-018-24530-9
  
2. Hughes, T. P., Kerry, J. T., Baird, A. H., Connolly, S. R.,   Dietzel, A., Eakin, C. M., Heron, S. F., Hoey, A. S.,   Hoogenboom, M. O., Liu, G., McWilliam, M. J., Pears, R. J.,   Pratchett, M. S., Skirving, W. J., Stella, J. S., &amp; Torda, G. (2018). Global warming transforms coral reef assemblages. Nature, 556(7702), 492-496. https://doi.org/10.1038/s41586-018-0041-2
  
3. Bozec, Y.-M., Rowell, D., Harrison, L., Gaskell, J., Hock, K.,   Callaghan, D., Gorton, R., Kovacs, E. M., Lyons, M., Mumby, P.,   &amp; Roelfsema, C. (2021). Baseline mapping to support reef restoration and   resilience-based management in the Whitsundays. https://doi.org/10.13140/RG.2.2.26976.20482
  
4. Bozec, Y.-M., Hock, K., Mason, R. A. B., Baird, M. E., Castro-Sanguino, C.,   Condie, S. A., Puotinen, M., Thompson, A., &amp; Mumby, P. J. (2022). Cumulative impacts across Australia&#39;s Great Barrier Reef: A mechanistic evaluation. Ecological Monographs, 92(1), e01494. https://doi.org/10.1002/ecm.1494
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/open-AIMS/ADRIA.jl/blob/894d691da86df90247def85cd9c7a50357ba6c3c/ADRIA/src/factors/const_params.jl#L1-L32" target="_blank" rel="noreferrer">source</a></Badge>

</details>

