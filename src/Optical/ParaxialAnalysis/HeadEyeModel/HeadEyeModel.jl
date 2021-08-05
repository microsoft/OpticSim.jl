
module HeadEye

using ....OpticSim
using ....OpticSim.Geometry
using ....OpticSim.Repeat
using ....OpticSim.Emitters

using Colors
using StaticArrays
using LinearAlgebra
import Makie, GeometryBasics, FileIO, MeshIO


include("Misc.jl")
include("Arrangement.jl")
include("HeadEyeClasses.jl")


end # module HeadEye
export HeadEye