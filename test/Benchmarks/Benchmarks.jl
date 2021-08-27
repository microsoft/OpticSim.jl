# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

module Benchmarks

using BenchmarkTools
using Unitful
using StaticArrays

using OpticSim
using OpticSim: Sphere, Ray, replprint, trace, Cylinder, AcceleratedParametricSurface, newton, Infinity, RayOrigin, NullInterface, IntervalPoint, intervalintersection, LensTrace, FresnelInterface, wavelength, mᵢandmₜ, refractedray, direction, origin, pathlength, LensAssembly, NoPower, power, snell, fresnel, reflectedray, BoundingBox
using OpticSim.Examples

include("../TestData/TestData.jl")

rayz() = Ray([0.0, 0.0, 10.0], [0.0, 0.0, -1.0])
perturbrayz() = Ray([0.0, 0.0, 10.0], [0.001, 0.001, -1.0])
rayx() = Ray([10.0, 0.0, 0.0], [-1.0, 0.0, 0.0])
bezierray() = Ray([0.5, 0.5, 10.0], [0.0, 0.0, -1.0])
optray() = OpticalRay(origin(rayz()), direction(rayz()), 1.0, 0.55)
perturboptray() = OpticalRay(origin(perturbrayz()), direction(perturbrayz()), 1.0, 0.55)

#############################
### BENCHMARK DEFINITIONS ###
#############################

# benchmarks are listed in the form: benchmark() = function, (args...)
# e.g. surfaceintersection, (surface, ray)

double_convex() = trace, (TestData.doubleconvex(), optray())
double_concave() = trace, (TestData.doubleconcave(), optray())
zernike_lens() = trace, (TestData.zernikesystem(), perturboptray())
aspheric_lens() = trace, (TestData.asphericsystem(), perturboptray())
cooke_triplet() = trace, (TestData.cooketriplet(), optray())
chebyshev_lens() = trace, (TestData.chebyshevsystem(), perturboptray())
gridsag_lens() = trace, (TestData.gridsagsystem(), perturboptray())
multi_hoe() = trace, (TestData.multiHOE(), OpticalRay(SVector(0.0, 3.0, 3.0), SVector(0.0, -1.0, -1.0), 1.0, 0.55))

system_benchmarks() = double_concave, double_convex, aspheric_lens, zernike_lens, chebyshev_lens, gridsag_lens, cooke_triplet, multi_hoe

triangle() = surfaceintersection, (Triangle(SVector(-1.0, -1.0, 0.0), SVector(1.0, 0.0, 0.0), SVector(0.0, 1.0, 0.0)), rayz())
sphericalcap() = surfaceintersection, (SphericalCap(1.0, π / 3), rayz())
sphere_outside() = surfaceintersection, (Sphere(0.5), rayz())
sphere_inside() = surfaceintersection, (Sphere(0.5), Ray([0.0, 0.0, 0.0], [0.0, 0.0, -1.0]))
plane() = surfaceintersection, (Plane(0.0, 0.0, 1.0, 0.0, 0.0, 0.0), rayz())
rectangle() = surfaceintersection, (Rectangle(1.0, 1.0, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, 0.0)), rayz())
ellipse() = surfaceintersection, (Ellipse(1.0, 1.0, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, 0.0)), rayz())
cylinder_parallel() = surfaceintersection, (Cylinder(1.0), rayz())
cylinder_perp() = surfaceintersection, (Cylinder(1.0), rayx())
bbox_parallel() = surfaceintersection, (BoundingBox(-1.0, 1.0, -1.0, 1.0, -5.0, 5.0), rayz())
bbox_perp() = surfaceintersection, (BoundingBox(-1.0, 1.0, -1.0, 1.0, -5.0, 5.0), rayx())

simple_surface_benchmarks() = plane, rectangle, ellipse, triangle, sphere_outside, sphere_inside, sphericalcap, cylinder_parallel, cylinder_perp, bbox_parallel, bbox_perp

zernike_surface1() = surfaceintersection, (AcceleratedParametricSurface(TestData.zernikesurface1(), 30), perturbrayz())
zernike_surface2() = surfaceintersection, (AcceleratedParametricSurface(TestData.zernikesurface3(), 30), perturbrayz())
qtype_surface() = surfaceintersection, (AcceleratedParametricSurface(TestData.qtypesurface1(), 30), perturbrayz())
bezier_surface() = surfaceintersection, (AcceleratedParametricSurface(TestData.beziersurface(), 30), bezierray())
chebyshev_surface1() = surfaceintersection, (AcceleratedParametricSurface(TestData.chebyshevsurface1(), 30), perturbrayz())
chebyshev_surface2() = surfaceintersection, (AcceleratedParametricSurface(TestData.chebyshevsurface2(), 30), perturbrayz())
gridsag_surface() = surfaceintersection, (AcceleratedParametricSurface(TestData.gridsagsurfacebicubic(), 30), perturbrayz())

accel_surface_benchmarks() = zernike_surface1, zernike_surface2, qtype_surface, bezier_surface, chebyshev_surface1, chebyshev_surface2, gridsag_surface

all_benchmarks() = (simple_surface_benchmarks()..., accel_surface_benchmarks()..., system_benchmarks()...)

#############################

function runbenchmark(b; kwargs...)
    f, args = b()
    qkwargs = [:($(keys(kwargs)[i]) = $(values(kwargs)[i])) for i in 1:length(kwargs)]
    if f === trace
        return @eval (@benchmark($f(($args)..., test = true), $(qkwargs...)))
    else
        return @eval (@benchmark($f(($args)...), setup = (OpticSim.emptyintervalpool!()), $(qkwargs...)))
    end
end

timings(f::Function) = timings(f()...)
function timings(functions::Vararg{Function})
    for func in functions
        println("timing $func")
        replprint(runbenchmark(func))
        println("\n\n***************\n")
    end
end

function runforpipeline(ismaster::Bool)
    io = open(ismaster ? "timings_master.csv" : "timings_pr.csv", "w")
    write(io, "Function, Memory, Allocs, Min Time, Mean Time, Max Time\n")
    for f in all_benchmarks()
        b = runbenchmark(f)
        write(io, [
            "$f",
            "$(BenchmarkTools.prettymemory(b.memory))",
            "$(b.allocs)",
            "$(BenchmarkTools.prettytime(minimum(b.times)))",
            "$(BenchmarkTools.prettytime(BenchmarkTools.mean(b.times)))",
            "$(BenchmarkTools.prettytime(maximum(b.times)))\n"
        ].join(", ")
        )
    end
    close(io)
end

########################################################

end #module

export Benchmarks
