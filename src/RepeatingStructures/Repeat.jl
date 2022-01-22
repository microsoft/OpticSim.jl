# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

module Repeat

export basismatrix

using Base: offset_if_vec
using StaticArrays:SVector,MVector,SMatrix,MMatrix
using DataFrames:DataFrame
using Colors
import LazySets
using LinearAlgebra:norm
import ..OpticSim #only LensletAssembly uses OpticSim. This doesn't seem like a great idea. Probably should move LensletAssembly somewhere else or at least remove the dependency.
import Unitful

include("Lattice.jl")
include("HexagonalLattice.jl")
include("RectangularLattice.jl")
include("Array.jl")
include("Cluster.jl")
include("Multilens/Multilens.jl")
include("LensletAssembly.jl")


end #module
export Repeat