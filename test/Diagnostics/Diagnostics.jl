# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

# Programs used to visualize output, profile code or perform debugging tasks, as opposed to unit testing
module Diagnostics

using OpticSim
using OpticSim: replprint
using OpticSim.Vis
using OpticSim.Optimization
using LinearAlgebra
using StaticArrays

import Unitful
import Plots

include("../TestData/TestData.jl")

function testbigfloat()
    λ = 550 * Unitful.u"nm"
    temp = 20 * Unitful.u"°C"

    TOLERANCE = 1e-8

    function printintersections(a::Type{T}) where {T<:Real}
        println("* Intersections for type $T *")
        r1 = OpticalRay(Vector{T}([0.0, 0.0, 1.0]), Vector{T}([0.0, 0.0, -1.0]), one(T), T(0.55))
        r2 = OpticalRay(Vector{T}([2.0, 2.0, 1.0]), Vector{T}([0.0, 0.0, -1.0]), one(T), T(0.55))
        r3 = OpticalRay(Vector{T}([5.0, 5.0, 1.0]), Vector{T}([0.0, 0.0, -1.0]), one(T), T(0.55))
        r4 = OpticalRay(Vector{T}([0.0, -5.0, 1.0]), Vector{T}([0.0, sind(5.0), -cosd(5.0)]), one(T), T(0.55))
        r5 = OpticalRay(Vector{T}([-5.0, -5.0, 1.0]), Vector{T}([sind(5.0), sind(2.0), -cosd(2.0) * cosd(5.0)]), one(T), T(0.55))

        a = TestData.doubleconvex(T, temperature = temp)
        for r in (r1, r2, r3, r4, r5)
            println(point(OpticSim.intersection(OpticSim.trace(a, r, test = true))))
        end
        println()
    end

    printintersections(Float64)
    printintersections(BigFloat)

    a = TestData.doubleconvex(BigFloat)
    r2 = OpticalRay(Vector{BigFloat}([2.0, 2.0, 1.0]), Vector{BigFloat}([0.0, 0.0, -1.0]), one(BigFloat), BigFloat(0.55))
    r3 = OpticalRay(Vector{BigFloat}([5.0, 5.0, 1.0]), Vector{BigFloat}([0.0, 0.0, -1.0]), one(BigFloat), BigFloat(0.55))
    r4 = OpticalRay(Vector{BigFloat}([0.0, -5.0, 1.0]), Vector{BigFloat}([0.0, sind(5.0), -cosd(5.0)]), one(BigFloat), BigFloat(0.55))
    r5 = OpticalRay(Vector{BigFloat}([-5.0, -5.0, 1.0]), Vector{BigFloat}([sind(5.0), sind(2.0), -cosd(2.0) * cosd(5.0)]), one(BigFloat), BigFloat(0.55))
    @assert isapprox(point(OpticSim.intersection(OpticSim.trace(a, r2, test = true))), [-0.06191521590711035, -0.06191521590711035, -67.8], atol = TOLERANCE)
    @assert isapprox(point(OpticSim.intersection(OpticSim.trace(a, r3, test = true))), [-0.2491105067897657, -0.2491105067897657, -67.8], atol = TOLERANCE)
    @assert isapprox(point(OpticSim.intersection(OpticSim.trace(a, r4, test = true))), [0.0, 5.639876913179362, -67.8], atol = TOLERANCE)
    @assert isapprox(point(OpticSim.intersection(OpticSim.trace(a, r5, test = true))), [5.75170097290395, 2.504152441922817, -67.8], atol = TOLERANCE)
end

function vistest(sys::AbstractOpticalSystem{Float64}; kwargs...)
    λ = 0.550
    r1 = OpticalRay([0.0, 0.0, 1.0], [0.0, 0.0, -1.0], 1.0, λ)
    r2 = OpticalRay([2.0, 2.0, 1.0], [0.0, 0.0, -1.0], 1.0, λ)
    r3 = OpticalRay([5.0, 5.0, 1.0], [0.0, 0.0, -1.0], 1.0, λ)
    r4 = OpticalRay([0.0, -5.0, 1.0], [0.0, 0.08715574274765818, -0.9961946980917454], 1.0, λ)
    r5 = OpticalRay([-5.0, -5.0, 1.0], [0.08715574274765818, -0.01738599476176408, -0.9960429728140486], 1.0, λ)
    raygen = RayListSource(r1, r2, r3, r4, r5)
    Vis.drawtracerays(sys, raygenerator = raygen, trackallrays = true, test = true; kwargs...)
    for (i, r) in enumerate(raygen)
        t = OpticSim.trace(sys, r, test = true)
        if t !== nothing
            println("=== RAY $i ===")
            println(OpticSim.point(t))
            println(OpticSim.pathlength(t))
            println()
        end
    end
end

function visualizerefraction()
    nₛ = SVector{3,Float64}([0.0, 0.0, 1.0])
    c = 1 / sqrt(2.0)
    r = Ray([0.0, c, c], [0.0, -c, -c])
    rray = OpticSim.refractedray(1.0, 1.5, nₛ, direction(r))
    Vis.draw(Ray([0.0, 0.0, 0.0], Array{Float64,1}(nₛ)), color = :black)
    Vis.draw!(r, color = :green)
    Vis.draw!(Ray([0.0, 0.0, 0.0], Array{Float64,1}(rray)), color = :red)
    Vis.draw!(Ray([0.0, 0.0, 0.0], Array{Float64,1}(-nₛ)), color = :yellow)

    temp = OpticSim.refractedray(1.0, 1.5, nₛ, direction(r))
    temp = [0.0, temp[2], -temp[3]]
    r = Ray([0.0, -temp[2], -temp[3]], temp)
    rray = OpticSim.refractedray(1.5, 1.0, nₛ, direction(r))
    Vis.draw!(Ray([0.0, 0.0, 0.0], Array{Float64,1}(nₛ)), color = :black)
    Vis.draw!(r, color = :green)
    Vis.draw!(Ray([0.0, 0.0, 0.0], Array{Float64,1}(rray)), color = :red)
    Vis.draw!(Ray([0.0, 0.0, 0.0], Array{Float64,1}(-nₛ)), color = :yellow)
end

function plotreflectedvsrefractedpower()
    lens = OpticSim.SphericalLens(OpticSim.GlassCat.SCHOTT.BAK50, 0.0, Inf64, Inf64, 5.0, 10.0)
    reflectpow = Array{Float64,1}(undef, 0)
    refractpow = Array{Float64,1}(undef, 0)
    green = 500 * Unitful.u"nm"
    glass = OpticSim.GlassCat.SCHOTT.BAK50

    incidentindex = OpticSim.GlassCat.index(glass, green)
    transmittedindex = OpticSim.GlassCat.index(OpticSim.GlassCat.Air, green)

    for θ in 0.0:0.01:(π / 2.0)
        origin = [0.0, tan(θ), 1.0]
        dir = -origin
        r = Ray(origin, dir)

        intsct = OpticSim.closestintersection(OpticSim.surfaceintersection(lens(), r))

        nml = SVector{3,Float64}(0.0, 0.0, 1.0)

        rdir = direction(r)
        (nᵢ, nₜ) = OpticSim.mᵢandmₜ(incidentindex, transmittedindex, nml, r)
        reflected = OpticSim.reflectedray(nml, rdir)
        refracted = OpticSim.refractedray(nᵢ, nₜ, nml, rdir)

        (sinθᵢ, sinθₜ) = OpticSim.snell(nml, direction(r), nᵢ, nₜ)
        (powᵣ, powₜ) = OpticSim.fresnel(nᵢ, nₜ, sinθᵢ, sinθₜ)
        push!(reflectpow, powᵣ)
        push!(refractpow, powₜ)
    end

    Plots.plot((0.0:0.01:(π / 2.0)), [reflectpow, refractpow], label = ["reflected" "refracted"])
end

end #module
export Diagnostics
