module Repeat

export basis

using Base: offset_if_vec
using StaticArrays:SVector,MVector

include("Lattice.jl")
include("HexagonalLattice.jl")
include("RectangularLattice.jl")
include("Array.jl")
include("Projections.jl")

end #module
export Repeat