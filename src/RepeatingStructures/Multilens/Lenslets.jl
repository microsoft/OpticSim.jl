# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.
module Lenslets
using LinearAlgebra
using Unitful
using Unitful.DefaultSymbols:mm,μm #Unitful.DefaultSymbols exports a variable T which can cause major confusion with type declarations since T is a common parametric type symbol. Import only what is needed.
using Colors
using StaticArrays
import DataFrames
import ..Repeat
import SpecialFunctions
import Plots
import ...OpticSim

include("HexClusters.jl")
include("HexTilings.jl")
include("Analysis.jl")
include("DisplayGeneration.jl")

end # module
export Lenslets