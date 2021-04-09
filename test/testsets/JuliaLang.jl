# FIXING ERRONEOUS UNBOUND TYPE ERRORS THAT OCCUR WITH VARARG
#=
The set will contain something like this:
    MultiHologramInterface(interfaces::Vararg{HologramInterface{T}, N}) where {T<:Real, N}
To get the signature run:
    methods(OpticSim.MultiHologramInterface).ms[I].sig            where I is typically 1 or 2 depending on the number of methods
This gives:
    Tuple{Type{MultiHologramInterface}, Vararg{HologramInterface{T}, N}} where N where T<:Real
We then have to use which to find the method:
    which(OpticSim.MultiHologramInterface, Tuple{Vararg{HologramInterface{T}, N}} where N where T<:Real)
and pop it from the set
=#

@testset "JuliaLang" begin
    # here we ignore the specific methods which we know are ok but are still failing
    # for some weird reason the unbound args check seems to fail for some (seemingly random) methods with Vararg arguments
    methods_to_ignore = Dict(
        OpticSim => VERSION >= v"1.6.0-DEV" ? [
            (OpticSim.LensAssembly, Tuple{Vararg{Union{CSGTree{T},LensAssembly{T},OpticSim.Surface{T}},N} where N} where {T<:Real}),
            (OpticSim.RayListSource, Tuple{Vararg{OpticalRay{T,3},N} where N} where {T<:Real}),
            (OpticSim.OpticalSourceGroup, Tuple{Vararg{OpticalRayGenerator{T},N} where N} where {T<:Real}),
            (OpticSim.PixelSource, Tuple{Vararg{P,N} where N} where {P<:OpticalRayGenerator{T}} where {T<:Real}),
            (OpticSim.MultiHologramInterface, Tuple{Vararg{HologramInterface{T},N}} where {N} where {T<:Real}),
        ] : [],
        OpticSim.Zernike => [],
        OpticSim.QType => [],
        OpticSim.Vis => [
            (OpticSim.Vis.drawcurves, Tuple{Vararg{Spline{P,S,N,M},N1} where N1} where {M} where {N} where {S} where {P}),
            (OpticSim.Vis.draw, Tuple{Vararg{S,N} where N} where {S<:Union{OpticSim.Surface{T},OpticSim.TriangleMesh{T}}} where {T<:Real}),
            (OpticSim.Vis.draw!, Tuple{OpticSim.Vis.MakieLayout.LScene,Vararg{S,N} where N} where {S<:Union{OpticSim.Surface{T},OpticSim.TriangleMesh{T}}} where {T<:Real})
        ],
        OpticSim.Examples => [],
        OpticSim.Chebyshev => [],
        OpticSim.GlassCat => [],
        OpticSim.Optimization => [],
    )

    for (mod, ignore_list) in methods_to_ignore
        unbound = setdiff(
            detect_unbound_args(mod),
            [which(method, signature) for (method, signature) in ignore_list]
        )
        # also ignore any generate methods created due to default or keyword arguments
        filter!(m -> !occursin("#", string(m.name)), unbound)
        @test isempty(unbound)

        ambiguous = detect_ambiguities(mod)
        @test isempty(ambiguous)
    end
end # testset JuliaLang
