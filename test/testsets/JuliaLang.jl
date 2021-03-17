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

@otestset "JuliaLang" begin
    # ensure there aren't any ambiguities or unbound args
    let unbound = Set{Method}(detect_unbound_args(OpticSim))
        # Pretty hacky way to ignore specific methods, for some weird reason the unbound args check seems to fail for some (seemingly random) methods with Vararg arguments
        # here we ignore the specific methods which we know are ok but are still failing
        if VERSION >= v"1.6.0-DEV"
            pop!(unbound, which(OpticSim.LensAssembly, Tuple{Vararg{Union{CSGTree{T},LensAssembly{T},Surface{T}},N} where N} where {T<:Real}))
            pop!(unbound, which(OpticSim.RayListSource, Tuple{Vararg{OpticalRay{T,3},N} where N} where {T<:Real}))
            pop!(unbound, which(OpticSim.OpticalSourceGroup, Tuple{Vararg{OpticalRayGenerator{T},N} where N} where {T<:Real}))
            pop!(unbound, which(OpticSim.PixelSource, Tuple{Vararg{P,N} where N} where {P<:OpticalRayGenerator{T}} where {T<:Real}))
            pop!(unbound, which(OpticSim.MultiHologramInterface, Tuple{Vararg{HologramInterface{T},N}} where {N} where {T<:Real}))
        end
        # and ignore any generate methods created due to default or keyword arguments
        for m in unbound
            if occursin("#", string(m.name))
                pop!(unbound, m)
            end
        end
        @test unbound == Set()
    end
    let ambiguous = Set{Method}(detect_ambiguities(OpticSim))
        @test ambiguous == Set()
    end
    let unbound = Set{Method}(detect_unbound_args(OpticSim.Zernike))
        @test unbound == Set()
    end
    let ambiguous = Set{Method}(detect_ambiguities(OpticSim.Zernike))
        @test ambiguous == Set()
    end
    let unbound = Set{Method}(detect_unbound_args(OpticSim.QType))
        @test unbound == Set()
    end
    let ambiguous = Set{Method}(detect_ambiguities(OpticSim.QType))
        @test ambiguous == Set()
    end
    let unbound = Set{Method}(detect_unbound_args(OpticSim.Vis))
        # Pretty hacky way to ignore specific methods, for some weird reason the unbound args check seems to fail for some (seemingly random) methods with Vararg arguments
        # here we ignore the specific methods which we know are ok but are still failing
        pop!(unbound, which(OpticSim.Vis.drawcurves, Tuple{Vararg{Spline{P,S,N,M},N1} where N1} where {M} where {N} where {S} where {P}))
        pop!(unbound, which(OpticSim.Vis.draw, Tuple{Vararg{S,N} where N} where {S<:Union{OpticSim.Surface{T},OpticSim.TriangleMesh{T}}} where {T<:Real}))
        pop!(unbound, which(OpticSim.Vis.draw!, Tuple{OpticSim.Vis.MakieLayout.LScene,Vararg{S,N} where N} where {S<:Union{OpticSim.Surface{T},OpticSim.TriangleMesh{T}}} where {T<:Real}))
        # and ignore any generate methods created due to default or keyword arguments
        for m in unbound
            if occursin("#", string(m.name))
                pop!(unbound, m)
            end
        end
        @test unbound == Set()
    end
    let ambiguous = Set{Method}(detect_ambiguities(OpticSim.Vis))
        @test ambiguous == Set()
    end
    let unbound = Set{Method}(detect_unbound_args(OpticSim.TestData))
        @test unbound == Set()
    end
    let ambiguous = Set{Method}(detect_ambiguities(OpticSim.TestData))
        @test ambiguous == Set()
    end
end # testset JuliaLang
