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
import LazySets
using Unitful
using Unitful.DefaultSymbols

include("LensletAssembly.jl")
include("HeadEyeModel/HeadEyeModel.jl")

end #module
export ParaxialAnalysis