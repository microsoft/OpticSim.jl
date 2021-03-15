module GlassCat

using Polynomials
using Plots
using StringEncodings
using Unitful
using StaticArrays
using Base: @.
import Unitful: Length, Temperature, Quantity, Units

include("GlassTypes.jl")
include("constants.jl")
include("search.jl")
include("utilities.jl")
export plot_indices, index, polyfit_indices, absairindex, absorption, info, glassid, glassname, glassforid, isair, findglass, modelglass, glassfromMIL, GlassID

# import built glass cat source files via deps/deps.jl
const deps_path = joinpath(@__DIR__, "..", "..", "deps", "deps.jl")
include(deps_path)

end # module
export GlassCat
