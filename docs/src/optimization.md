# Optimization

We are currently in the early stages of implementing optimization for lens surfaces.

We have had success using [Optim.jl](https://julianlsolvers.github.io/Optim.jl/stable/) and [ForwardDiff.jl](http://www.juliadiff.org/ForwardDiff.jl/stable/) among other packages, though ForwardDiff.jl gradient and hessian pre-compilation can be very slow for complex systems.

Other optimization packages _should_ integrate easily with the system too: [JuMP.jl](https://github.com/jump-dev/JuMP.jl), [Ipopt.jl](https://github.com/jump-dev/Ipopt.jl) and [NLOpt.jl](https://github.com/JuliaOpt/NLopt.jl) are some options.

Merit functions must currently be implemented by the user, this is quite straight-forward as [`trace`](@ref) returns the vast majority of information that could be needed in the form of a [`LensTrace`](@ref) which hits the detector (or `nothing` if the ray doesn't reach the detector).

There are some helper functions implemented for [`AxisymmetricOpticalSystem`](@ref)s which can make optimization of basic systems much easier:

```@docs
OpticSim.Optimization
OpticSim.Optimization.optimizationvariables
OpticSim.Optimization.updateoptimizationvariables
```

It is of course possible to write your own optimization loop for more complex (i.e. non-AxisymmetricOpticalSystem) systems and this _should_ work without issue.

## Example

```julia
using OpticSim
using ForwardDiff
using Optim

function objective(a::AbstractVector{T}, b::AxisymmetricOpticalSystem{T}, samples::Int = 3) where {T}
    # RMSE spot size
    system = Optimization.updateoptimizationvariables(b, a)
    # distribute rays evenly across entrance pupil using HexapolarField. The latest version of OpticSim no longer supports HexapolarField. 
    field = HexapolarField(system, collimated = true, samples = samples)
    error = zero(T)
    hits = 0
    for r in field
        traceres = trace(system, r, test = true)
        if traceres !== nothing # ignore rays which miss
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
    # if hits == 0 returns 0 - not ideal!
    return error
end

start, lower, upper = Optimization.optimizationvariables(system)
optimobjective = arg -> objective(arg, system)
gcfg = ForwardDiff.GradientConfig(optimobjective, start, ForwardDiff.Chunk{1}()) # speed up ForwardDiff significantly
hcfg = ForwardDiff.HessianConfig(optimobjective, start, ForwardDiff.Chunk{1}())
g! = (G, x) -> ForwardDiff.gradient!(G, optimobjective, x, gcfg)
h! = (H, x) -> ForwardDiff.hessian!(H, optimobjective, x, hcfg)
df = TwiceDifferentiable(optimobjective, g!, h!, start)
dfc = TwiceDifferentiableConstraints(lower, upper) # constrain the optimization to avoid e.g. thickness < 0
res = optimize(df, dfc, start, algo, Optim.Options(show_trace = true, iterations = 100, allow_f_increases = true))
final = Optim.minimizer(res)
new_system = Optimization.updateoptimizationvariables(system, final)
```
