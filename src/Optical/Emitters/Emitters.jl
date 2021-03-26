#=
Emitters are defined by Pixels and SpatialLayouts. An emitter has a spectrum, and an optical power distribution over the hemisphere. 
These are intrinsic physical properties of the emitter.
=#

module Emitters
export Geometry, Spectrum, Directions, Origins, AngularPower, Sources

import ..OpticSim
# import Makie
# import Makie.AbstractPlotting
# import Makie.AbstractPlotting.MakieLayout

# using LinearAlgebra
# # import  Distributions
# using DataFrames
# using StaticArrays
# using Unitful

struct MissingImplementationException <: Exception
    msg::String
end

# time measuring macro with a caption
macro timeCap(caption, ex)
    quote
        print("[Timing] [$($caption)] ")
        @time $(esc(ex))
    end
end
export @timeCap

# defining name placeholders to override in nested modules
generate() = 0 
export generate
visual_size() = 0
export visual_size
apply() = 0
export apply

include("Geometry.jl")
include("Spectrum.jl")
include("Directions.jl")
include("Origins.jl")
include("AngularPower.jl")
include("Sources.jl")
include("Visualization.jl")

end # module Emitters
