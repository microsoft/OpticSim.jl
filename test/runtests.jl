# MIT License

# Copyright (c) Microsoft Corporation.

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE

using Test
using FiniteDifferences
using StaticArrays
using LinearAlgebra
using Suppressor
using Random
using Unitful
using Plots

using OpticSim
# curve/surface imports
using OpticSim: findspan, makemesh, knotstoinsert, coefficients, inside, quadraticroots, tobeziersegments, evalcsg, makiemesh
# interval imports
using OpticSim: α, halfspaceintersection, positivehalfspace, lower, upper, EmptyInterval, rayorigininterval, intervalcomplement, intervalintersection, intervalunion, RayOrigin, Infinity, Intersection
include("TestData/TestData.jl")
# bounding box imports
using OpticSim: doesintersect
# RBT imports
using OpticSim: rotmat, rotmatd
# lens imports
using OpticSim: reflectedray, snell, refractedray, trace, intersection, pathlength, power
#Bounding volume hierarchy imports
# using OpticSim: partition!

const COMP_TOLERANCE = 25 * eps(Float64)
const RTOLERANCE = 1e-10
const ATOLERANCE = 10 * eps(Float64)
const SEED = 12312487
const ALL_TESTS = isempty(ARGS) || "all" in ARGS || "All" in ARGS || "ALL" in ARGS

"""Optional testset: if the name of the testset is passed as an arg then it will be executed. If none are specified, or "all" is given as an arg then all optional testsets will be executed """
macro otestset(name, expr)
    quote
        if ALL_TESTS || $name in ARGS
            @testset $name begin
                $expr
            end
        end
    end
end

"""Evalute all functions not requiring arguments in a given module and test they don't throw anything"""
macro test_all_no_arg_functions(m)
    quote
        for n in names($m, all = true)
            # get all the valid functions
            if occursin("#", string(n)) || string(n) == string(nameof($m))
                continue
            end
            f = Core.eval($m, n)
            # get the methods of this function
            for meth in methods(f)
                # if this method has no args then try and evaluate it
                if meth.nargs == 1
                    # suppress STDOUT to prevent loads of stuff spamming the console
                    @test (@suppress_out begin
                        f()
                    end;
                    true)
                end
            end
        end
    end
end

include("testsets/JuliaLang.jl")
# include("testsets/BVH.jl")
include("testsets/TestData.jl")
# include("testsets/Examples.jl")  # slow
include("testsets/General.jl")
include("testsets/SurfaceDefs.jl")
include("testsets/Intersection.jl")
include("testsets/Lenses.jl")
include("testsets/OpticalSystem.jl")
include("testsets/Emitters.jl")  # TODO
include("testsets/Comparison.jl")
include("testsets/Visualization.jl")

include("Benchmarks/Benchmarks.jl")
include("testsets/Allocations.jl")

include("testsets/GlassCat.jl")
