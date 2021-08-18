# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

module HeadEye

using ..OpticSim
using ..OpticSim.Geometry
using ..OpticSim.Repeat
using ..OpticSim.Emitters
using ..OpticSim.Repeat

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