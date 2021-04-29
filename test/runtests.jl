# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

using Test
using FiniteDifferences
using StaticArrays
using LinearAlgebra
using Suppressor
using Random
using Unitful
using Plots

using OpticSim
# Geometry Module
using OpticSim.Geometry
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
const TESTSET_DIR = "testsets"

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

include("Benchmarks/Benchmarks.jl")

alltestsets = [
    "JuliaLang",
    # "BVH",
    # "Examples", # slow
    "General",
    "SurfaceDefs",
    "Intersection",
    "Lenses",
    "OpticalSystem",
    "Emitters", # TODO
    "Comparison",
    "Visualization",
    "Allocations",
    "GlassCat",
    "Transform"
]

runtestsets = ALL_TESTS ? alltestsets : intersect(alltestsets, ARGS)
include.([joinpath(TESTSET_DIR, "$(testset).jl") for testset in runtestsets])
