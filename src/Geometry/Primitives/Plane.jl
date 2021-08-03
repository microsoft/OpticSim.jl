# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

"""
    Plane{T,N} <: ParametricSurface{T,N}

Infinite planar surface where the positive normal side is outside the surface.

By default this will not create any geometry for visualization, the optional `vishalfsizeu` and `vishalfsizev` arguments can be used to draw the plane as a rectangle for visualization **note that this does not fully represent the surface**.
In this case, the rotation of the rectangle around the normal to the plane is defined by `visvec` - `surfacenormal×visvec` is taken as the vector along the u axis.

```julia
Plane(surfacenormal::SVector{N,T}, pointonplane::SVector{N,T}; interface::NullOrFresnel{T} = nullinterface(T), vishalfsizeu::T = 0.0, vishalfsizev::T = 0.0, visvec::SVector{N,T} = [0.0, 1.0, 0.0])
Plane(nx::T, ny::T, nz::T, x::T, y::T, z::T; interface::NullOrFresnel{T} = nullinterface(T), vishalfsizeu::T = 0.0, vishalfsizev::T = 0.0, visvec::SVector{N,T} = [0.0, 1.0, 0.0])
```
"""
struct Plane{T,N} <: ParametricSurface{T,N}
    normal::SVector{N,T}
    d::T
    pointonplane::SVector{N,T}
    surface_id::UUID
    interface::NullOrFresnel{T}
    # below only for visualization purposes
    vishalfsizeu::T
    vishalfsizev::T
    visuvec::SVector{N,T}
    visvvec::SVector{N,T}

    function Plane(surfacenormal::AbstractArray{T,1}, pointonplane::AbstractArray{T,1}; interface::NullOrFresnel{T} = NullInterface(T), vishalfsizeu::T = zero(T), vishalfsizev::T = zero(T), visvec::AbstractArray{T,1} = [0.0, 1.0, 0.0]) where {T}
        @assert length(surfacenormal) == length(pointonplane) == length(visvec)
        N = length(surfacenormal)
        return Plane(SVector{N,T}(surfacenormal), SVector{N,T}(pointonplane), interface = interface, vishalfsizeu = vishalfsizeu, vishalfsizev = vishalfsizev, visvec = SVector{N,T}(visvec))
    end

    function Plane(surfacenormal::SVector{N,T}, pointonplane::SVector{N,T}; interface::NullOrFresnel{T} = NullInterface(T), vishalfsizeu::T = zero(T), vishalfsizev::T = zero(T), visvec::SVector{N,T} = SVector{3,T}(0.0, 1.0, 0.0)) where {T<:Real,N}
        norml = normalize(surfacenormal)
        d = dot(norml, pointonplane)
        if abs(dot(visvec, norml)) == one(T)
            visvec = SVector{3,T}(1.0, 0.0, 0.0)
        end
        uvec = normalize(cross(normalize(visvec), norml))
        vvec = normalize(cross(norml, uvec))
        return new{T,N}(norml, d, pointonplane, uuid1(Random.GLOBAL_RNG), interface, vishalfsizeu, vishalfsizev, uvec, vvec)
    end

    function Plane(nx::T, ny::T, nz::T, x::T, y::T, z::T; interface::NullOrFresnel{T} = NullInterface(T), vishalfsizeu::T = zero(T), vishalfsizev::T = zero(T), visvec = SVector{3,T}(0.0, 1.0, 0.0)) where {T<:Real}
        return Plane(SVector{3,T}(nx, ny, nz), SVector{3,T}(x, y, z), interface = interface, vishalfsizeu = vishalfsizeu, vishalfsizev = vishalfsizev, visvec = visvec)
    end
end
export Plane

Base.show(io::IO, a::Plane{T}) where {T<:Real} = print(io, "Plane{$T}($(a.pointonplane), $(normal(a)), $(interface(a)))")

surface_id(a::Plane{T,N}) where {T<:Real,N} = a.surface_id
interface(a::Plane{T,N}) where {T<:Real,N} = a.interface
normal(pln::Plane{T,N}) where {T<:Real,N} = pln.normal

inside(pln::Plane{T,3}, p::SVector{3,T}) where {T<:Real} = dot(normal(pln), p) - pln.d < zero(T)

onsurface(pln::Plane{T,3}, p::SVector{3,T}) where {T<:Real} = samepoint(dot(normal(pln), p), pln.d)

distancefromplane(p::Plane{T,N}, point::SVector{N,T}) where {N,T<:Real} = dot(normal(p), point) - p.d

uvrange(::Type{Plane{T,N}}) where {T<:Real,N} = ((-one(T), one(T)), (-one(T), one(T)))
point(p::Plane{T}, u::T, v::T) where {T<:Real} = p.pointonplane + p.vishalfsizeu * u * p.visuvec + p.vishalfsizev * v * p.visvvec
normal(p::Plane{T}, ::T, ::T) where {T<:Real} = normal(p)

function surfaceintersection(pln::Plane{T,N}, r::AbstractRay{T,N}) where {T<:Real,N}
    n̂ = normal(pln)
    d = direction(r)
    o = origin(r)
    nd = dot(n̂, d)
    if samepoint(nd, zero(T))
        # ray and plane are parallel
        if inside(pln, o) || onsurface(pln, o)
            return rayorigininterval(Infinity(T))
        else
            # no intersection only if the ray is strictly outside of the palne
            return EmptyInterval(T)
        end
    end
    t = (pln.d - dot(n̂, o)) / nd
    if t < zero(T)
        if inside(pln, o)
            # if the ray starts 'inside' the surface then we want to return a ray so intersection works
            return rayorigininterval(Infinity(T))
        else
            return EmptyInterval(T) # no ray plane intersection
        end
    end
    temp = Intersection(t, point(r, t), n̂, zero(T), zero(T), surface_id(pln), interface(pln))
    if nd < zero(T)
        return positivehalfspace(temp)
    else
        return rayorigininterval(temp)
    end
end

function BoundingBox(pln::Plane{T,3}) where {T<:Real}
    # TODO! this is far from ideal, we should try and do something better for intersection with non-axis-algined planes
    # valid for axis aligned planes, otherwise we have to assume an infinite bounding box
    if normal(pln) === SVector{3,T}(0, 0, 1)
        return BoundingBox(typemin(T), typemax(T), typemin(T), typemax(T), typemin(T), pln.pointonplane[3])
    elseif normal(pln) === SVector{3,T}(0, 0, -1)
        return BoundingBox(typemin(T), typemax(T), typemin(T), typemax(T), pln.pointonplane[3], typemax(T))
    elseif normal(pln) === SVector{3,T}(0, 1, 0)
        return BoundingBox(typemin(T), typemax(T), typemin(T), pln.pointonplane[2], typemin(T), typemax(T))
    elseif normal(pln) === SVector{3,T}(0, -1, 0)
        return BoundingBox(typemin(T), typemax(T), pln.pointonplane[2], typemax(T), typemin(T), typemax(T))
    elseif normal(pln) === SVector{3,T}(1, 0, 0)
        return BoundingBox(typemin(T), pln.pointonplane[1], typemin(T), typemax(T), typemin(T), typemax(T))
    elseif normal(pln) === SVector{3,T}(-1, 0, 0)
        return BoundingBox(pln.pointonplane[1], typemax(T), typemin(T), typemax(T), typemin(T), typemax(T))
    else
        return BoundingBox(typemin(T), typemax(T), typemin(T), typemax(T), typemin(T), typemax(T))
    end
end
