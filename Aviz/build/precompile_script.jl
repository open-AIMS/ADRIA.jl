using CSV, DataFrames, GeoDataFrames, GeoInterface
using Statistics, Distributions
using Makie.GeometryBasics
using ADRIA, Makie

precompile(CSV.read, (String, DataFrame))
precompile(GeoDataFrames.read, (String,))

precompile(ADRIA.load_domain, (String, String))
precompile(ADRIA.load_results, (String,))

precompile(Figure, ())
precompile(Axis, (Figure,))
precompile(lines!, (Axis, Array))
precompile(scatter!, (Axis, Array))
precompile(series!, (Axis, Array))

precompile(hist!, (Axis, Array))
precompile(hist!, (Axis, Array, Int64))
precompile(violin!, (Axis, Array, Array))
