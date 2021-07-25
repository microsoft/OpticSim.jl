# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

module ParaxialAnalysis

using StaticArrays: promote_tuple_eltype
using ..OpticSim
using ..OpticSim.Geometry
using StaticArrays
using LinearAlgebra
using LazySets

include("LensletAssembly.jl")
include("SphericalPolygon.jl")
include("BeamEnergy.jl")

end #module
export ParaxialAnalysis