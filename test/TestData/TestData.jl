# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

"""
Contains all complex data used for testing and benchmarking.
"""
module TestData

using OpticSim
using OpticSim: tobeziersegments # try to use only exported functions so this list should stay short
using OpticSim.Geometry
using OpticSim.Emitters
using OpticSim.GlassCat

using StaticArrays
using DataFrames: DataFrame
using Unitful
using Images
using Base: @.
using LinearAlgebra

include("curves.jl")
include("surfaces.jl")
include("lenses.jl")
include("systems.jl")
include("other.jl")

end # module TestData
export TestData
