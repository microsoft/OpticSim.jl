# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

"""
Optimization interface consists of two functions `optimizationvariables` and `updateoptimizationvariables`.
`optimizationvariables` packs variables to be optimized into a vector.
`updateoptimizationvariables` receives a vector of variables and creates a new optical system with the variable values.
"""
module Optimization

using ..OpticSim
using ..OpticSim: detectorsize, temperature, pressure
using ..OpticSim.Emitters

using ForwardDiff
using DataFrames
using Optim
using JuMP
using Ipopt
using Zygote
using NLopt

include("OptimizationVariables.jl")
include("Constraints.jl")
include("TestCases.jl")

end #module
export Optimization
