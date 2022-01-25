# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

module Vis


using ..OpticSim
using ..OpticSim: euclideancontrolpoints, evalcsg, vertex, makiemesh, detector, centroid, lower, upper, intervals, α
using ..OpticSim.Geometry
using ..OpticSim.Repeat

using Unitful
using ImageView
using Images
using ColorTypes
using ColorSchemes
using StaticArrays
using LinearAlgebra
import Makie
import GeometryBasics
import Plots
import Luxor
using FileIO

include("Visualization.jl")
include("Emitters.jl")
include("VisRepeatingStructures.jl")

end # module Vis
export Vis
