# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

"""Contains example usage of the features in the OpticSim.jl package."""
module Examples
using ..OpticSim
using ..OpticSim.Vis
using ..OpticSim.Geometry
using ..OpticSim.Emitters
using ..OpticSim.Emitters.Spectrum
using ..OpticSim.Emitters.Origins
using ..OpticSim.Emitters.Directions
using ..OpticSim.Emitters.Sources
using ..OpticSim.GlassCat

using StaticArrays
using DataFrames
using Images
using Unitful
using Plots
using LinearAlgebra

include("docs.jl")
include("other.jl")

end #module Examples
export Examples
