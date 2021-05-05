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

# include("UI/AxisSymmetric.jl") # old version
include("UI/BlinkUtils.jl")
include("UI/AxisSymmetricOptimization.jl")
include("UI/AxisSymmetricOptimizationUI.jl")

end #module
export Optimization

#Notes 
# User defined objective or constraint functions must be registered with JuMP before use. If we intend to automatically set up and run the optimization in JuMP then we will need an interface that includes a list of objective functions
# and a list of constraint functions. We also need to define an interface for the end programmer for defining constructor functions that build the optical system given some set of inputs and that extract the optimizable variables
# from the system. These two must be consistent. 

struct OptimizationVariable{T} <: Number #maybe or maybe not
    value::T
    optimize::Bool
end

#interface for user constructor function takes OptimizationVariable args instead of float args. User provides a vector guess and a function that takes the vector and peels off the elements to be used in the actual function call.
# makesystem(a::OptimizationVariable{T}...) where{T<:Number}