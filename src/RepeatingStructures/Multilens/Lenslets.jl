# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.
module Lenslets
using LinearAlgebra
import Unitful
using Unitful:uconvert,ustrip
using Unitful.DefaultSymbols:mm,μm,nm #Unitful.DefaultSymbols exports a variable T which can cause major confusion with type declarations since T is a common parametric type symbol. Import only what is needed.
using Colors
using StaticArrays
import DataFrames
import ..Repeat
import SpecialFunctions
import Plots
import ...OpticSim
import ...OpticSim.Geometry
import ...OpticSim.Repeat:region,HexBasis1,HexBasis3,LatticeBasis,LatticeCluster,clustersize,ClusterWithProperties

include("HexClusters.jl")
include("HexTilings.jl")
include("Analysis.jl")
include("DisplayGeneration.jl")
include("LensletAssignment.jl")


end # module
export Lenslets