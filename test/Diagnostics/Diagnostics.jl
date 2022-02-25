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
    raygen = Emitters.Sources.RayListSource(r1, r2, r3, r4, r5)
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
    lens = OpticSim.SphericalLens(OpticSim.Examples.Examples_BAK50, 0.0, Inf64, Inf64, 5.0, 10.0)
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

### OPTIMIZATION TESTING

using FiniteDifferences
using ForwardDiff
using DataFrames: DataFrame
using Optim
using JuMP
using Ipopt
using Zygote
using NLopt

doubleconvexprescription() = DataFrame(SurfaceType = ["Object", "Standard", "Standard", "Image"], Radius = [(Inf64), 60.0, -60.0, (Inf64)], Thickness = [(Inf64), (10.0), (77.8), missing], Material = [OpticSim.GlassCat.Air, OpticSim.Examples.Examples_N_BK7, OpticSim.GlassCat.Air, missing], SemiDiameter = [(Inf64), (9.0), (9.0), (15.0)])

function doubleconvex(a::AbstractVector{T}; detpix::Int = 100) where {T<:Real}
    frontradius = a[1]
    rearradius = a[2]
    #! format: off
    AxisymmetricOpticalSystem{T}(DataFrame(
        SurfaceType = ["Object", "Standard", "Standard", "Image"],
        Radius = [T(Inf64), frontradius, rearradius, T(Inf64)],
        Conic = [missing, -1.0, 1.0, missing],
        Thickness = [T(Inf64), T(10.0), T(77.8), missing],
        Material = [OpticSim.GlassCat.Air, OpticSim.Examples.Examples_N_BK7, OpticSim.GlassCat.Air, missing],
        SemiDiameter = [T(Inf64), T(9.0), T(9.0), T(15.0)],
    ), detpix, detpix, T, temperature = OpticSim.GlassCat.TEMP_REF_UNITFUL, pressure = OpticSim.GlassCat.PRESSURE_REF)
    #! format: on
end

function RMS_spot_size(a::AbstractVector{T}, b::AxisymmetricOpticalSystem{T}, samples::Int = 3) where {T}
    # RMSE spot size
    lens = Optimization.updateoptimizationvariables(b, a)
    field = HexapolarField(lens, collimated = true, samples = samples)
    error = zero(T)
    hits = 0
    for r in field
        traceres = OpticSim.trace(lens, r, test = true)
        if traceres !== nothing
            hitpoint = point(traceres)
            if abs(hitpoint[1]) > eps(T) && abs(hitpoint[2]) > eps(T)
                dist_to_axis = hitpoint[1]^2 + hitpoint[2]^2
                error += dist_to_axis
            end
            hits += 1
        end
    end
    if hits > 0
        error = sqrt(error / hits)
    end
    return error
end

function testfinitedifferences()
    lens = Examples.doubleconvex(60.0, -60.0)
    start = Optimization.optimizationvariables(lens)

    u = 60.0
    println("time to compute objective function $(@time RMS_spot_size([u,-60.0],lens))")
    fdm = central_fdm(3, 1)
    f1(u) = RMS_spot_size([u, -60.0], lens)
    fu = 0.0
    for i in 1:100
        fu = fdm(f1, 60.0)
    end
    return fu
end
export testfinitedifferences

function testoptimization(; lens = Examples.cooketriplet(), constrained = false, algo = constrained ? IPNewton() : LBFGS(), chunk_size = 1, samples = 3)
    @info "NaN safe: $(ForwardDiff.NANSAFE_MODE_ENABLED)"

    start, lower, upper = Optimization.optimizationvariables(lens)
    optimobjective = (arg) -> RMS_spot_size(arg, lens, samples)

    if chunk_size == 0
        gcfg = ForwardDiff.GradientConfig(optimobjective, start)
        hcfg = ForwardDiff.HessianConfig(optimobjective, start)
    else
        gcfg = ForwardDiff.GradientConfig(optimobjective, start, ForwardDiff.Chunk{chunk_size}())
        hcfg = ForwardDiff.HessianConfig(optimobjective, start, ForwardDiff.Chunk{chunk_size}())
        @info "$(length(start)) optimization variables, chunk size = $chunk_size"
    end

    g! = (G, x) -> ForwardDiff.gradient!(G, optimobjective, x, gcfg)
    h! = (H, x) -> ForwardDiff.hessian!(H, optimobjective, x, hcfg)

    @info "Starting objective function $(optimobjective(start))"
    if constrained && !all(isinf.(abs.(lower))) && !all(isinf.(upper))
        if !(algo isa Optim.ConstrainedOptimizer)
            @error "This algorithm can't handle constraints, use IPNewton"
            return
        end
        @info "Constraints: $lower, $upper"
        df = TwiceDifferentiable(optimobjective, g!, h!, start)
        dfc = TwiceDifferentiableConstraints(lower, upper)
        res = optimize(df, dfc, start, algo, Optim.Options(show_trace = true, iterations = 100, allow_f_increases = true))
    else
        if constrained
            @warn "No constraints to apply, using unconstrained optimization"
        end
        res = optimize(optimobjective, g!, h!, start, algo, Optim.Options(show_trace = true, iterations = 100, allow_f_increases = true))
    end

    # println("== Computing gradient")
    # println("First call:")
    # @time g!(zeros(length(start)), start)
    # println("Second call:")
    # @time g!(zeros(length(start)), start)

    # println("== Computing hessian")
    # println("First call:")
    # @time h!(zeros(length(start), length(start)), start)
    # println("Second call:")
    # @time h!(zeros(length(start), length(start)), start)

    # println()

    final = Optim.minimizer(res)
    newlens = Optimization.updateoptimizationvariables(lens, final)
    @info "Result: $final"

    field = HexapolarField(newlens, collimated = true, samples = samples)
    Vis.drawtraceimage(newlens, raygenerator = field)
    Vis.drawtracerays(newlens, raygenerator = field, trackallrays = true, test = true)
end

function testnlopt()
    lens = Examples.doubleconvex(60.0, -60.0)
    # lens = Examples.doubleconvex()
    # lens = Examples.telephoto(6,.5)
    # Vis.drawtraceimage(lens)
    start = OpticSim.optimizationvariables(lens)

    optimobjective = (arg) -> RMS_spot_size(arg, lens)
    println("starting objective function $(optimobjective(start))")

    model = Model(NLopt.optimizer)
    @variable(model, rad1, start = 60.0)
    @variable(model, rad2, start = -60.0)
    register(model, :f, 2, f, autodiff = true)
    @NLconstraint(model, 30.0 <= rad1 <= 85.0)
    @NLconstraint(model, -100.0 <= rad1 <= -55.0)
    @NLobjective(model, Min, f(rad1, rad2))
    optimize!(model)

    println("rad1 $(value(rad1)) rad2 $(value(rad2))")

end

function testJuMP()
    function innerfunc(radius1::T, radius2::T) where {T<:Real}

        lens = Examples.doubleconvex(radius1, radius2)

        field = HexapolarField(lens, collimated = true, samples = 10)
        error = zero(T)
        hits = 0
        for r in field
            traceres = OpticSim.trace(lens, r, test = true)

            if !(nothing === traceres)
                hitpoint = point(traceres)
                if abs(hitpoint[1]) > eps(T) && abs(hitpoint[2]) > eps(T)
                    dist_to_axis = sqrt(hitpoint[1]^2 + hitpoint[2]^2)
                    error += dist_to_axis
                end
                hits += 1
            end
        end

        # if isa(radius1,ForwardDiff.Dual)
        #     println("radius1 $(realpart(radius1)) radius2 $(realpart(radius2))")
        # end

        error /= hits
        # println("error in my code $error")

        # Vis.drawtracerays(doubleconvex(radiu1,radius2),trackallrays = true, test = true)

        return error
    end

    function f(radius1::T, radius2::T) where {T<:Real}
        innerfunc(radius1, radius2)
        # try
        #     innerfunc(radius1, radius2)
        # catch err
        #     return zero(T)
        # end
    end

    # Vis.drawtracerays(doubleconvex(60.0,-60.0), trackallrays = true, test = true)

    model = Model(Ipopt.Optimizer)
    set_time_limit_sec(model, 10.0)
    register(model, :f, 2, f, autodiff = true)
    @variable(model, 10.0 <= rad1 <= 85.0, start = 60.0)
    @variable(model, -150 <= rad2 <= 20.0, start = -10.0)
    @NLobjective(model, Min, f(rad1, rad2))
    optimize!(model)

    println("rad1 $(value(rad1)) rad2 $(value(rad2))")
    Vis.drawtracerays(Examples.doubleconvex(value(rad1), value(rad2)), trackallrays = true, test = true)
end

function simpletest()
    lens = Examples.cooketriplet()
    # lens = Examples.doubleconvex(60.0,-60.0)
    testgenericoptimization(lens, RMS_spot_size)
end

function RMS_spot_size(lens, x::T...) where {T}
    lens = Optimization.updateoptimizationvariables(lens, collect(x))
    field = HexapolarField(lens, collimated = true, samples = 10)
    error = zero(T)
    hits = 0
    for r in field
        traceres = OpticSim.trace(lens, r, test = true)

        if !(nothing === traceres)
            hitpoint = point(traceres)
            if abs(hitpoint[1]) > eps(T) && abs(hitpoint[2]) > eps(T)
                dist_to_axis = sqrt(hitpoint[1]^2 + hitpoint[2]^2)
                error += dist_to_axis
            end
            hits += 1
        end
    end

    error /= hits

    return error
end

function testZygoteGradient(lens::S, objective::Function) where {S<:AbstractOpticalSystem}
    vars = Optimization.optimizationvariables(lens)
    gradient = objective'(vars)
    println(gradient)
end

function testgenericoptimization(lens::S, objective::Function) where {S<:AbstractOpticalSystem}
    vars = Optimization.optimizationvariables(lens)
    function newobjective(a)
        objective(lens, a...)
    end

    numvars = length(vars)
    model = Model(Ipopt.Optimizer)
    set_time_limit_sec(model, 10.0)
    register(model, :newobjective, numvars, newobjective, autodiff = true)
    # @variable(model, bounds[i][1] <= x[i=1:numvars] <= bounds[i][2],start = vars)
    @variable(model, x[i = 1:numvars], start = vars[i])

    @NLobjective(model, Min, newobjective(x...))
    # d = NLPEvaluator(model)
    # MOI.initialize(d, [:Grad])
    # println("objective $(MOI.eval_objective(d, vars))")

    # ∇f = zeros(length(vars))
    # println(vars)
    # MOI.eval_objective_gradient(d, ∇f, vars)
    # println("gradient $∇f")
    # println("starting forward diff")
    # time = @time ForwardDiff.gradient(newobjective,vars)
    # println("forward diff time1 $time")
    # time = @time ForwardDiff.gradient(newobjective,vars)
    # println("forward diff time2 $time")
    println("starting reverse diff")
    time = @time ReverseDiff.gradient(newobjective, vars)
    println("reverse diff time1 $time")
    time = @time ReverseDiff.gradient(newobjective, vars)
    println("reverse diff time2 $time")

    println("ending forward diff")
    # fdm = central_fdm(2,1)
    # replprint( @benchmark(FiniteDifferences.grad($fdm,$newobjective,$vars...)))
    # println("finite difference gradient $(FiniteDifferences.grad(fdm,newobjective,vars...))")
    # optimize!(model)
    lens = Optimization.updateoptimizationvariables(lens, value.(x))

    Vis.drawtracerays(lens, trackallrays = true, test = true)
end

function testIpopt()

    model = Model(Ipopt.Optimizer)

    @variable(model, x1, start = 1)
    @variable(model, x2, start = 5)
    @variable(model, x3, start = 5)
    @variable(model, x4, start = 1)

    f(x1, x2, x3, x4) = x1 * x4 * (x1 + x2 + x3) + x3

    register(model, :f, 4, f, autodiff = true)
    @NLobjective(model, Min, f(x1, x2, x3, x4))
    @NLconstraint(model, x1 * x2 * x3 * x4 >= 25)
    @NLconstraint(model, x1^2 + x2^2 + x3^2 + x4^2 == 40)

    optimize!(model)

    println("x1 $(value(x1)) x2 $(value(x2)) x3 $(value(x3)) x4 $(value(x4))")
end

function testtypesignature()
    a = OpticSim.Sphere(1.0)
    b = (a ∪ a) ∩ (a ∪ a)
end

function testoptimizationvariables()
    lens = Examples.cooketriplet()
    vars = OpticSim.optimizationvariables(lens)
    vars .= [Float64(i) for i in 1:length(vars)]
    newlens = OpticSim.updateoptimizationvariables(lens, vars)
    println("original lens")
    show(lens)
    println("updated lens")
    show(newlens)
end

end #module
export Diagnostics
