# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

export LensAssembly, elements, access
export LensTrace, intersection

"""
    LensAssembly{T<:Real}

Structure which contains the elements of the optical system, these can be [`CSGTree`](@ref) or [`Surface`](@ref) objects.

In order to prevent type ambiguities bespoke structs are created for each possible number of elements e.g. `LensAssembly3`.
These are parameterized by the types of the elements to prevent ambiguities.
Basic surface types such as [`Rectangle`](@ref) (which can occur in large numbers) are stored independently in `Vector`s, so type paramters are only needed for CSG objects.

Each struct looks like this:
```julia
struct LensAssemblyN{T,T1,T2,...,TN} <: LensAssembly{T}
    axis::SVector{3,T}
    planarshapes::Vector{PlanarShape{T}}
    paraxials::Vector{ParaxialLens{T}}
    E1::T1
    E2::T2
    ...
    EN::TN
end
```

Where `Ti <: Union{Surface{T},CSGTree{T}}`.

To create a LensAssembly object the following functions can be used:
```julia
LensAssembly(elements::Vararg{Union{Surface{T},CSGTree{T},LensAssembly{T}}}; axis = SVector(0.0, 0.0, 1.0)) where {T<:Real}
```
"""
abstract type LensAssembly{T<:Real} end

LensAssemblyElement{T} = Union{Surface{T}, CSGTree{T}, LensAssembly{T}}

Base.show(io::IO, ::Type{<:LensAssembly}) = print(io, "LensAssembly")
Base.show(io::IO, a::LensAssembly{T}) where {T<:Real} = print(io, "LensAssembly{$T}($(length(elements(a))) elements)")

# in order to prevent type ambiguities lens assembly must be paramaterised by its elements, but we want lens assembly to
# have an arbitrary number of elements as such we have to use macros to create lens assembly structs of a given size,
# and methods to handle them in a type unambiguous way we may often end up with lots of simple surfaces (Rectangle,
# Ellipse) in the assembly, and we can store these more efficiently in Vectors
macro lensassembly_constructor(N)
    # generate the type definition for a lens assembly of size N
    typesig = [:($(Symbol("T$(i)")) <: Union{Surface{T},CSGTree{T}}) for i in 1:N]
    fields = [:($(Symbol("E$(i)"))::$(Symbol("T$(i)"))) for i in 1:N]
    name = Symbol("LensAssembly$N")
    return esc(quote
        struct $(name){T,$(typesig...)} <: LensAssembly{T}
            axis::SVector{3,T}
            planarshapes::Vector{PlanarShape{T}}
            paraxials::Vector{ParaxialLens{T}}
            $(fields...)
        end
    end)
end

macro lensassembly_intersection(N)
    # generates the _closestintersection() function for a lens assembly of size N, this should never be called directly,
    # instead call closestintersection() in all cases
    type = :($(Symbol("LensAssembly$(N)")))
    fieldnames = [Symbol("E$(i)") for i in 1:N]
    checks = [quote
        intvl = surfaceintersection(obj.$(fieldnames[i]), r)
        if !(intvl isa EmptyInterval)
            intsct = closestintersection(intvl::Union{Interval{T},DisjointUnion{T}})
            if intsct !== nothing
                αcurr = α(intsct)
                if αcurr < αmin
                    αmin = αcurr
                    closest = intsct
                end
            end
        end
    end for i in 1:N]
    return esc(quote
        function _closestintersection(obj::$(type){T}, r::AbstractRay{T,N})::Union{Nothing,Intersection{T,N}} where {T<:Real,N}
            αmin = typemax(T)
            closest = nothing
            for element in obj.planarshapes
                intvl = surfaceintersection(element, r)
                if !(intvl isa EmptyInterval)
                    intsct = closestintersection(intvl::Union{Interval{T},DisjointUnion{T}})
                    if intsct !== nothing
                        αcurr = α(intsct)
                        if αcurr < αmin
                            αmin = αcurr
                            closest = intsct
                        end
                    end
                end
            end
            for element in obj.paraxials
                intvl = surfaceintersection(element, r)
                if !(intvl isa EmptyInterval)
                    intsct = closestintersection(intvl::Union{Interval{T},DisjointUnion{T}})
                    if intsct !== nothing
                        αcurr = α(intsct)
                        if αcurr < αmin
                            αmin = αcurr
                            closest = intsct
                        end
                    end
                end
            end
            $(checks...)
            return closest
        end
    end)
end

macro lensassembly_elements(N)
    # generates the _elements() function for a lens assembly of size N, this should never be called directly, instead call elements() in all cases
    type = :($(Symbol("LensAssembly$(N)")))
    fieldnames = [Symbol("E$(i)") for i in 1:N]
    tuple = [:(obj.$(fieldnames[i])) for i in 1:N]
    return esc(quote
        function _typed_elements(obj::$(type){T}) where {T<:Real}
            return tuple($(tuple...))
        end
    end)
end

function compile_lensassembly(N)
    if !isdefined(OpticSim, :($(Symbol("LensAssembly$N"))))
        @eval @lensassembly_constructor($N)
        @eval @lensassembly_intersection($N)
        @eval @lensassembly_elements($N)
    end
end

# pregenerate a small number at compile time
compile_lensassembly.(1:PREGENERATED_LENS_ASSEMBLY_SIZE)

function LensAssembly(elements::Vararg{LensAssemblyElement{T}}; axis::SVector{3,T} = SVector{3,T}(0.0, 0.0, 1.0)) where {T<:Real}
    # make the actual object
    typed_elements = Vector{Union{Surface{T}, CSGTree{T}}}(undef, 0)
    planarshapes = Vector{PlanarShape{T}}(undef, 0)
    paraxials = Vector{ParaxialLens{T}}(undef, 0)
    for e in elements
        if e isa PlanarShape{T}
            push!(planarshapes, e)
        elseif e isa ParaxialLens{T}
            push!(paraxials, e)
        elseif e isa Surface{T} || e isa CSGTree{T}
            push!(typed_elements, e)
        elseif e isa LensAssembly{T}
            # unpack the lens assembly
            planarshapes = vcat(planarshapes, e.planarshapes)
            paraxials = vcat(paraxials, e.paraxials)
            append!(typed_elements, typed_elements(e))
        else
            error("Unrecognized LensAssembly element type $typeof(e)")
        end
    end

    # make the methods for this N
    N = length(typed_elements)
    compile_lensassembly(N)

    type = :($(Symbol("LensAssembly$(N)")))
    eltypes = typeof.(typed_elements)
    naxis = normalize(axis)
    return eval(:($(type){$T,$(eltypes...)}($naxis, $planarshapes, $paraxials, $(typed_elements...))))
end

function typed_elements(ass::LensAssembly{T}) where {T<:Real}
    # under certain circumstances we run into the world age problem (see https://discourse.julialang.org/t/how-to-bypass-the-world-age-problem/7012)
    # to avoid this we need to use invokelatest(), but this is type ambiguous and ruins performance, and in most cases isn't necessary
    # so we can use try/catch to only call it when we need (and warn appropriately), retaining performance in almost all cases
    # NEVER call _typed_elements directly, always use this method
    try
        _typed_elements(ass)
    catch MethodError
        @warn "First time call to LensAssembly of previously unused size, performance will be worse than normal." maxlog = 1
        Base.invokelatest(_typed_elements, ass)
    end
end

elements(ass::LensAssembly{T}) where {T<:Real} = [ass.planarshapes..., ass.paraxials..., typed_elements(ass)...]

function closestintersection(ass::LensAssembly{T}, r::AbstractRay{T,N})::Union{Nothing,Intersection{T,N}} where {T<:Real,N}
    # under certain circumstances we run into the world age problem (see https://discourse.julialang.org/t/how-to-bypass-the-world-age-problem/7012)
    # to avoid this we need to use invokelatest(), but this is type ambiguous and ruins performance, and in most cases isn't necessary
    # so we can use try/catch to only call it when we need (and warn appropriately), retaining performance in almost all cases
    # NEVER call _closestintersection directly, always use this method
    try
        _closestintersection(ass, r)
    catch MethodError
        @warn "First time call to LensAssembly of previously unused size, performance will be worse than normal." maxlog = 1
        Base.invokelatest(_closestintersection, ass, r)
    end
end

function BoundingBox(la::LensAssembly{T}) where {T<:Real}
    es = elements(la)
    bbox = BoundingBox(es[1])
    for b in es[2:end]
        if b isa ParametricSurface{T} || b isa CSGTree{T}
            bbox = union(bbox, b)
        end
    end
    return bbox
end

axis(a::LensAssembly{T}) where {T<:Real} = a.axis

#############################################################################

"""
    LensTrace{T<:Real,N}

Contains an intersection point and the ray segment leading to it from within an optical trace.
The ray carries the path length, power, wavelength, number of intersections and source number, all of which are accessible directly on this class too.


Has the following accessor methods:
```julia
ray(a::LensTrace{T,N}) -> OpticalRay{T,N}
intersection(a::LensTrace{T,N}) -> Intersection{T,N}
power(a::LensTrace{T,N}) -> T
wavelength(a::LensTrace{T,N}) -> T
pathlength(a::LensTrace{T,N}) -> T
point(a::LensTrace{T,N}) -> SVector{N,T}
uv(a::LensTrace{T,N}) -> SVector{2,T}
sourcenum(a::LensTrace{T,N}) -> Int
nhits(a::LensTrace{T,N}) -> Int
```
"""
struct LensTrace{T<:Real,N}
    ray::OpticalRay{T,N}
    intersection::Intersection{T,N}
end

ray(a::LensTrace{T,N}) where {T<:Real,N} = a.ray
intersection(a::LensTrace{T,N}) where {T<:Real,N} = a.intersection
power(a::LensTrace{T,N}) where {T<:Real,N} = power(ray(a))
wavelength(a::LensTrace{T,N}) where {T<:Real,N} = wavelength(ray(a))
pathlength(a::LensTrace{T,N}) where {T<:Real,N} = pathlength(ray(a))
point(a::LensTrace{T,N}) where {T<:Real,N} = point(intersection(a))
uv(a::LensTrace{T,N}) where {T<:Real,N} = uv(intersection(a))
sourcenum(a::LensTrace{T,N}) where {T<:Real,N} = sourcenum(ray(a))
nhits(a::LensTrace{T,N}) where {T<:Real,N} = nhits(ray(a))

function Base.print(io::IO, a::LensTrace{T,N}) where {T,N}
    println(io, "Ray\n$(ray(a))")
    println(io, "Lensintersection\n$(intersection(a))")
end

struct NoPower end
const nopower = NoPower()

# isnopower(::NoPower) = true
# isnopower(::Any) = false

#############################################################################

"""
    trace(assembly::LensAssembly{T}, r::OpticalRay{T}, temperature::T = 20.0, pressure::T = 1.0; trackrays = nothing, test = false)

Returns the ray as it exits the assembly in the form of a [`LensTrace`](@ref) object if it hits any element in the assembly, otherwise `nothing`.
Recursive rays are offset by a small amount (`RAY_OFFSET`) to prevent it from immediately reintersecting the same lens element.

`trackrays` can be passed an empty vector to accumulate the `LensTrace` objects at each intersection of `ray` with a surface in the assembly.
"""
function trace(assembly::LensAssembly{T}, r::OpticalRay{T,N}, temperature::T = T(OpticSim.GlassCat.TEMP_REF), pressure::T = T(OpticSim.GlassCat.PRESSURE_REF); trackrays::Union{Nothing,Vector{LensTrace{T,N}}} = nothing, test::Bool = false, recursion::Int = 0)::Union{Nothing,NoPower,LensTrace{T,N}} where {T<:Real,N}
    if power(r) < POWER_THRESHOLD || recursion > TRACE_RECURSION_LIMIT
        return nopower
    end
    intsct = closestintersection(assembly, r)
    if intsct === nothing
        return nothing
    else
        surfintsct = point(intsct)
        nml = normal(intsct)
        opticalinterface = interface(intsct)
        λ = wavelength(r)

        if VERSION < v"1.6.0-DEV"
            # TODO REMOVE
            temp = @unionsplit OpticalInterface T opticalinterface processintersection(opticalinterface, surfintsct, nml, r, temperature, pressure, test, recursion == 0)
        else
            temp = processintersection(opticalinterface, surfintsct, nml, r, temperature, pressure, test, recursion == 0)
        end

        if temp === nothing
            return nopower
        else
            raydirection, raypower, raypathlength = temp
            if trackrays !== nothing
                # we want the power that is hitting the surface, i.e. power on incident ray, but the path length at this intersection
                # i.e. with the path length for this intersection added and power modulated by absorption only
                push!(trackrays, LensTrace(OpticalRay(ray(r), raypower, λ, opl = raypathlength, nhits = nhits(r), sourcenum = sourcenum(r), sourcepower = sourcepower(r)), intsct))
            end
            offsetray = OpticalRay(surfintsct + RAY_OFFSET * raydirection, raydirection, raypower, λ, opl = raypathlength, nhits = nhits(r) + 1, sourcenum = sourcenum(r), sourcepower = sourcepower(r))
            res = trace(assembly, offsetray, temperature, pressure, trackrays = trackrays, test = test, recursion = recursion + 1)
            if res === nothing
                return LensTrace(OpticalRay(surfintsct, raydirection, raypower, λ, opl = raypathlength, nhits = nhits(r), sourcenum = sourcenum(r), sourcepower = sourcepower(r)), intsct)
            else
                return res
            end
        end
    end
end

# """approximate focal length of simple double convex lens. For plano convex set one of the radii to Inf"""
focallength(index::Float64, rf::Float64, rb::Float64) = 1 / ((index - 1) * (1 / rf - 1 / rb))
