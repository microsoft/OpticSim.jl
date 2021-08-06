
module HeadEye

using ..OpticSim
using ..OpticSim.Geometry
using ..OpticSim.Repeat
using ..OpticSim.Emitters

using Colors
using StaticArrays, Statistics
using LinearAlgebra
import Makie, GeometryBasics, FileIO, MeshIO


include("Misc.jl")
include("Arrangement.jl")
include("HeadEyeClasses.jl")
include("Examples.jl")


end # module HeadEye
export HeadEye