# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.
module Multilens

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
import ...OpticSim:plane_from_points,centroid,pointonplane,shape,focallength
import ...OpticSim.Geometry
import ...OpticSim.Repeat:region,HexBasis1,HexBasis3,LatticeBasis,LatticeCluster,clustersize,ClusterWithProperties,cluster_coordinates_from_tile_coordinates, AbstractLatticeCluster
using ...OpticSim.Data

include("HexClusters.jl")
include("HexTilings.jl")
include("Analysis.jl")
include("DisplayGeneration.jl")
include("LensletAssignment.jl")
include("Example.jl")


end # module
export Multilens