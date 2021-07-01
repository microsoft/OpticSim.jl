module Repeat

using StaticArrays:SVector,MVector

include("Lattice.jl")
include("HexagonalLattice.jl")
include("RectangularLattice.jl")
include("Array.jl")
include("Projections.jl")

end #module
export Repeat