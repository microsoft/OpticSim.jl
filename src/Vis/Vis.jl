# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

module Vis

using Base: byte_string_classify
using ..OpticSim
using ..OpticSim: euclideancontrolpoints, evalcsg, vertex, makiemesh, detector, centroid, lower, upper, intervals, α
using ..OpticSim.Geometry

using Unitful
using ImageView
using Images
using ColorTypes
using ColorSchemes
using StaticArrays
using LinearAlgebra
import Makie
import Makie.AbstractPlotting
import Makie.AbstractPlotting.MakieLayout
import Plots
import Luxor
using FileIO

# If using precompiled system image (which we always are) you have to run AbstractPlotting.__init__() after loading Makie
# during precompilation, the display stack gets shuffled around such that the Makie display does not take priority.
# See https://discourse.julialang.org/t/makie-doesnt-display-plot-when-using-a-custom-julia-sysimage/38515.
__init__() = AbstractPlotting.__init__()

include("Visualization.jl")
include("Emitters.jl")
include("VisRepeatingStructures.jl")

end # module Vis
export Vis
