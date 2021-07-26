# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

module ParaxialAnalysis

using StaticArrays: promote_tuple_eltype
using ..OpticSim
using ..OpticSim.Geometry
using ..OpticSim.Repeat
using StaticArrays
using LinearAlgebra
using LazySets

include("LensletAssembly.jl")
include("SphericalPolygon.jl")
include("BeamEnergy.jl")
include("paraxial_analysis_test.jl")

end #module
export ParaxialAnalysis