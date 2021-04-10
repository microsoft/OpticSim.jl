#=
Emitters are defined by Pixels and SpatialLayouts. An emitter has a spectrum, and an optical power distribution over the hemisphere. 
These are intrinsic physical properties of the emitter.
=#

module Emitters

using ...OpticSim, ..Geometry
using LinearAlgebra

# defining name placeholders to override in nested modules
generate() = 0 
export generate
visual_size() = 0
export visual_size
apply() = 0
export apply

# define some utility functions 
right(t::OpticSim.Geometry.Transform{T}) where {T<:Real} = normalize(Vec3(t[1,1], t[2,1], t[3,1]))
up(t::OpticSim.Geometry.Transform{T}) where {T<:Real} = normalize(Vec3(t[1,2], t[2,2], t[3,2]))
forward(t::OpticSim.Geometry.Transform{T}) where {T<:Real} = normalize(Vec3(t[1,3], t[2,3], t[3,3]))
OpticSim.origin(t::OpticSim.Geometry.Transform{T}) where {T<:Real} = Vec3(t[1,4], t[2,4], t[3,4])
export right, up, forward

include("Spectrum.jl")
include("Directions.jl")
include("Origins.jl")
include("AngularPower.jl")
include("Sources.jl")

# exporting stuff from the Emitters module
export Spectrum, Directions, Origins, AngularPower, Sources

end # module Emitters

export Emitters
