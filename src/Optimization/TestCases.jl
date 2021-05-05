# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

doubleconvexprescription() = DataFrame(Surface = [:Object, 1, 2, :Image], Radius = [(Inf64), 60.0, -60.0, (Inf64)], Thickness = [(Inf64), (10.0), (77.8), missing], Material = [OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, missing], SemiDiameter = [(Inf64), (9.0), (9.0), (15.0)])

function doubleconvex(a::AbstractVector{T}; detpix::Int = 100) where {T<:Real}
    frontradius = a[1]
    rearradius = a[2]
    #! format: off
    AxisymmetricOpticalSystem{T}(DataFrame(
        Surface = [:Object, 1, 2, :Image],
        Radius = [T(Inf), frontradius, rearradius, T(Inf)],
        Conic = [missing, -1.0, 1.0, missing],
        Thickness = [T(Inf), T(10.0), T(77.8), missing],
        Material = [OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, missing],
        SemiDiameter = [T(Inf), T(9.0), T(9.0), T(15.0)],
    ), detpix, detpix, T, temperature = OpticSim.GlassCat.TEMP_REF_UNITFUL, pressure = OpticSim.GlassCat.PRESSURE_REF)
    #! format: on
end



# function testfinitedifferences()
#     lens = Examples.doubleconvex(60.0, -60.0)
#     start = Optimization.optimizationvariables(lens)

#     u = 60.0
#     println("time to compute objective function $(@time RMS_spot_size([u,-60.0],lens))")
#     fdm = central_fdm(3, 1)
#     f1(u) = RMS_spot_size([u, -60.0], lens)
#     fu = 0.0
#     for i in 1:100
#         fu = fdm(f1, 60.0)
#     end
#     return fu
# end
# export testfinitedifferences

function testOptim(; lens = Examples.cooketriplet(), constrained = false, algo = constrained ? IPNewton() : LBFGS(), chunk_size = 1, samples = 3)
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
        res = Optim.optimize(df, dfc, start, algo, Optim.Options(show_trace = true, iterations = 100, allow_f_increases = true))
    else
        if constrained
            @warn "No constraints to apply, using unconstrained optimization"
        end
        res = Optim.optimize(optimobjective, g!, h!, start, algo, Optim.Options(show_trace = true, iterations = 100, allow_f_increases = true))
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

#This example doesn't work. Leaving it here on the off chance we decide to revisit nlopt as an optimization platform
# function testnlopt()
#     lens = Examples.doubleconvex(60.0, -60.0)
#     # lens = Examples.doubleconvex()
#     # lens = Examples.telephoto(6,.5)
#     # Vis.drawtraceimage(lens)
#     start = Optimization.optimizationvariables(lens)

#     optimobjective = (arg) -> RMS_spot_size(arg, lens)
#     println("starting objective function $(optimobjective(start))")

#     model = Model(NLopt.optimizer)
#     @variable(model, rad1, start = 60.0)
#     @variable(model, rad2, start = -60.0)
#     register(model, :f, 2, f, autodiff = true)
#     @NLconstraint(model, 30.0 <= rad1 <= 85.0)
#     @NLconstraint(model, -100.0 <= rad1 <= -55.0)
#     @NLobjective(model, Min, f(rad1, rad2))
#     NLOpt.optimize!(model)

#     println("rad1 $(value(rad1)) rad2 $(value(rad2))")

# end

function testJuMP()
    function innerfunc(radius1::T, radius2::T) where {T<:Real}

        lens = Examples.doubleconvex(radius1, radius2)

        field = Sources.Source(transform = Geometry.translation(0.0,0.0,10.0), origins = Origins.Hexapolar(5,5.0,5.0),directions = Directions.Constant(0.0,0.0,-1.0))
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
    set_time_limit_sec(model, 100.0)
    register(model, :f, 2, f, autodiff = true)
    @variable(model, 10.0 <= rad1 <= 85.0, start = 60.0)
    @variable(model, -150 <= rad2 <= 20.0, start = -10.0)
    @NLobjective(model, Min, f(rad1, rad2))
    JuMP.optimize!(model)

    println("rad1 $(value(rad1)) rad2 $(value(rad2))")
    Vis.drawtracerays(Examples.doubleconvex(value(rad1), value(rad2)), trackallrays = true, test = true)
end

function simpletest()
    lens = Examples.cooketriplet()
    # lens = Examples.doubleconvex(60.0,-60.0)
    testgenericoptimization(lens, RMS_spot_size)
end



function testZygoteGradient(lens::S, objective::Function) where {S<:OpticalSystem}
    vars = Optimization.optimizationvariables(lens)
    gradient = objective'(vars)
    println(gradient)
end

function testgenericoptimization(lens::S, objective::Function) where {S<:OpticalSystem}
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
    println("here")
    a = OpticSim.Sphere(1.0)
    b = csgintersection(csgunion(a, a), csgunion(a, a))
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