# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

"""
    Primitive{T<:Real}

`T` is the number type used to represent the primitive,  e.g., `Float64`.
Primitives are the basic elements which can be stored in bounding volume hierarchies and include surfaces and CSG objects

**Must** implement the following:
```julia
boundingbox(a::Primitive{T})::BoundingBox{T}
centroid(a::Primitive{T})::SVector{3,T}
```
"""
abstract type Primitive{T<:Real} end


"""
    Surface{T<:Real}

`T` is the number type used to represent the surface, e.g., `Float64`.
Basic `Surface`s are _not_ valid CSG objects, they function only in a stand-alone capacity.

**Must** implement the following:
```julia
surfaceintersection(surface::Surface{T}, ray::AbstractRay{T,3}) -> Union{EmptyInterval{T},Interval{T}}
normal(surface::Surface{T}) -> SVector{3,T}
interface(surface::Surface{T}) -> OpticalInterface{T}
makemesh(surface::Surface{T}) -> TriangleMesh{T}
```

In a conventional ray tracer the surface intersection function would only return the first surface the ray intersects. Because our ray tracer does CSG operations the surface intersection function intersects the ray with all leaf surfaces which are part of the CSG tree. 

Each leaf surface returns one or more 1D intervals along the ray. These intervals contain the part of the ray which is inside the surface. The intervals computed at the leaves are propagated upward through the CSG tree and the CSG operations of union, intersection, and difference are applied to generate new intervals which are themselves propagated upward.

The result is a union of 1D intervals, which may be disjoint, a single interval, or empty. The union of intervals represents the parts of the ray which are inside the CSG object.

Inside is well defined for halfspaces such as cylinders and spheres which divide space into two parts, but not for Bezier or NURBS patches which generally do not enclose a volume.  For surfaces which are not halfspaces the notion of inside is defined locally by computing the angle between the incoming ray and the normal of the surface at the point of intersection. All surfaces must be defined so that the normal points to the outside of the surface. 

A negative dot product between the incoming ray and the normal indicates the ray is coming from the outside of the surface and heading toward the inside. A positive dot product indicates the ray is coming from the inside of the surface and heading toward the outside.

Intervals are defined along the ray which is being intersected with the surface, so they are one dimensional. For example, assume we have a ray with origin o on the outside of a plane and an intersection with the plane at point int = o + td where t is a scalar and d is the unit direction of the ray. The inside interval will be (Intersection(t),Infinity). This interval begins at the intersection point on the plane and continues to positive infinity. The Intersection struct stores both the parametric value t and the 3D point of intersection to make various operations more efficient. But the interval operations only depend on the parametric value t.

If the origin o is on the inside of the plane then the inside interval will be (RayOrigin,Intersection(t)). Only the part of the ray from the ray origin to the intersection point is inside the plane. 

It is the programmer's responsibility to return Interval results from surfaceintersection that maintain these properties.

The following must be implemented only if the surface is being used as a detector
```julia
uv(surface::Surface{T}, p::SVector{3,T}) -> SVector{2,T}
uvtopix(surface::Surface{T}, uv::SVector{2,T}, imsize::Tuple{Int,Int}) -> Tuple{Int,Int}
onsurface(surface::Surface{T}, p::SVector{3,T}) -> Bool
```
"""
abstract type Surface{T<:Real} <: Primitive{T} end
export Surface, surfaceintersection, normal, interface

"""
    ParametricSurface{T,N} <: Surface{T}

`T` is the number type used to represent the surface, e.g., `Float64`.
`N` is the dimension of the space the surface is embedded in.
`ParametricSurface`s are valid CSG objects, in some cases (where analytic intersection isn't possible) they must be wrapped in an [`AcceleratedParametricSurface`](@ref) for use.

**Must** implement the following:
```julia
uv(surface::ParametricSurface{T,N}, p::SVector{N,T}) -> SVector{2,T}
uvrange(surface::ParametricSurface{T,N}) -> Tuple{Tuple{T,T},Tuple{T,T}}
point(surface::ParametricSurface{T,N}, u::T, v::T) -> SVector{N,T}
partials(surface::ParametricSurface{T,N}, u::T, v::T) -> Tuple{SVector{N,T}, SVector{N,T}}
normal(surface::ParametricSurface{T,N}, u::T, v::T) -> SVector{N,T}
inside(surface::ParametricSurface{T,N}, p: :SVector{N,T}) -> Bool
onsurface(surface::ParametricSurface{T,N}, p::SVector{N,T}) -> Bool
surfaceintersection(surface::ParametricSurface{T,N}, AbstractRay::Ray{T,N}) -> Union{EmptyInterval{T},Interval{T},DisjointUnion{T}}
interface(surface::ParametricSurface{T,N}) -> OpticalInterface{T}
```
"""
abstract type ParametricSurface{S<:Real,N} <: Surface{S} end
export ParametricSurface, point, partials, uvrange, inside, onsurface, uv

# all subclasses must implement one of these at least...
"""
    point(surf::ParametricSurface{T}, u::T, v::T) -> SVector{3,T}
    point(surf::ParametricSurface{T}, uv::SVector{2,T}) -> SVector{3,T}

Returns the 3D point on `surf` at the given uv coordinate.
"""
point(s::ParametricSurface{T}, u::T, v::T) where {T<:Real} = point(s, SVector{2,T}(u, v))
point(s::ParametricSurface{T}, uv::SVector{2,T}) where {T<:Real} = point(s, uv[1], uv[2])
"""
    normal(surf::ParametricSurface{T}, u::T, v::T) -> SVector{3,T}
    normal(surf::ParametricSurface{T}, uv::SVector{2,T}) -> SVector{3,T}

Returns the normal to `surf` at the given uv coordinate.
"""
normal(s::ParametricSurface{T}, u::T, v::T) where {T<:Real} = normal(s, SVector{2,T}(u, v))
normal(s::ParametricSurface{T}, uv::SVector{2,T}) where {T<:Real} = normal(s, uv[1], uv[2])
"""
    partials(surf::ParametricSurface{T}, u::T, v::T) -> (SVector{3,T}, SVector{3,T})
    partials(surf::ParametricSurface{T}, uv::SVector{2,T}) -> (SVector{3,T}, SVector{3,T})

Returns a tuple of the 3D partial derivatives of `surf` with respect to u and v at the given uv coordinate.
"""
partials(s::ParametricSurface{T}, u::T, v::T) where {T<:Real} = partials(s, SVector{2,T}(u, v))
partials(s::ParametricSurface{T}, uv::SVector{2,T}) where {T<:Real} = partials(s, uv[1], uv[2])
"""
    uv(surf::ParametricSurface{T}, p::SVector{3,T}) -> SVector{2,T}
    uv(surf::ParametricSurface{T}, x::T, y::T, z::T) -> SVector{2,T}

Returns the uv coordinate on `surf` of a point, `p`, in 3D space.
If `onsurface(surf, p)` is false then the behavior is undefined, it may return an inorrect uv, an invalid uv, NaN or crash.
"""
uv(s::ParametricSurface{T,3}, x::T, y::T, z::T) where {T<:Real} = uv(s, SVector{3,T}(x, y, z))
uv(s::ParametricSurface{T,3}, p::SVector{3,T}) where {T<:Real} = uv(s, p[1], p[2], p[3])
"""
    inside(surf::ParametricSurface{T}, p::SVector{3,T}) -> Bool
    inside(surf::ParametricSurface{T}, x::T, y::T, z::T) -> Bool

Tests whether a 3D point in world space is _inside_ `surf`.
"""
inside(s::ParametricSurface{T,3}, x::T, y::T, z::T) where {T<:Real} = inside(s, SVector{3,T}(x, y, z))
inside(s::ParametricSurface{T,3}, p::SVector{3,T}) where {T<:Real} = inside(s, p[1], p[2], p[3])
"""
    onsurface(surf::ParametricSurface{T}, p::SVector{3,T}) -> Bool
    onsurface(surf::ParametricSurface{T}, x::T, y::T, z::T) -> Bool

Tests whether a 3D point in world space is _on_ `surf`.
"""
onsurface(s::ParametricSurface{T,3}, x::T, y::T, z::T) where {T<:Real} = onsurface(s, SVector{3,T}(x, y, z))
onsurface(s::ParametricSurface{T,3}, p::SVector{3,T}) where {T<:Real} = onsurface(s, p[1], p[2], p[3])
"""
    uvrange(s::ParametricSurface)
    uvrange(::Type{S}) where {S<:ParametricSurface}

Returns a tuple of the form: `((umin, umax), (vmin, vmax))` specifying the limits of the parameterisation for this surface type.
Also implemented for some `Surface`s which are not `ParametricSurface`s (e.g. `Rectangle`).
"""
uvrange(::S) where {S<:ParametricSurface} = uvrange(S)

"""
    samplesurface(surf::ParametricSurface{T,N}, samplefunction::Function, numsamples::Int = 30)

Sample a parametric surface on an even `numsamples`×`numsamples` grid in UV space with provided function
"""
function samplesurface(surf::ParametricSurface{T,N}, samplefunction::Function, numsamples::Int = 30) where {T<:Real,N}
    urange, vrange = uvrange(surf)
    rt = typeof(samplefunction(surf, urange[1], vrange[1]))
    samples = Vector{rt}(undef, (numsamples + 1)^2)
    ustep = (urange[2] - urange[1]) / numsamples
    vstep = (vrange[2] - vrange[1]) / numsamples
    for ui in 0:numsamples
        for vi in 0:numsamples
            u = urange[1] + ui * ustep
            v = vrange[1] + vi * vstep
            samples[ui * (numsamples + 1) + vi + 1] = samplefunction(surf, u, v)
        end
    end
    return samples
end
export samplesurface

"""
    triangulate(surf::ParametricSurface{S,N}, quads_per_row::Int, extensionu::Bool = false, extensionv::Bool = false, radialu::Bool = false, radialv::Bool = false)

Create an array of triangles representing the parametric surface where vertices are sampled on an even grid in UV space.
The surface can be extended by 1% in u and v separately, and specifying either u or v as being radial - i.e. determining the radius on the surface e.g. rho for zernike - will result in that dimension being sampled using sqwrt so that area of triangles is uniform. The extension will also only apply to the maximum in this case.
"""
function triangulate(surf::ParametricSurface{T,N}, subdivisons::Int, extensionu::Bool = false, extensionv::Bool = false, radialu::Bool = false, radialv::Bool = false) where {T,N}
    triangles = newintrianglepool!(T)
    (umin, umax), (vmin, vmax) = uvrange(surf)

    if extensionu
        if !radialu
            umin -= TRIANGULATION_EXTENSION
        end
        umax += TRIANGULATION_EXTENSION
    end

    if extensionv
        if !radialv
            vmin -= TRIANGULATION_EXTENSION
        end
        vmax += TRIANGULATION_EXTENSION
    end

    # if we are using this for intersection then we need to expand the triangles very slightly to
    # avoid misses due to floating point precision
    expansion = extensionu || extensionv ? TRIANGULATION_EXPANSION : zero(T)

    for ui in 0:(subdivisons - 1)
        for vj in 0:(subdivisons - 1)
            # this can be inefficient as we evaluate each point twice
            # in practice we usually use extension in which case no point is the same anyway
            u1 = T(ui / subdivisons)
            u2 = T((ui + 1) / subdivisons)
            v1 = T(vj / subdivisons)
            v2 = T((vj + 1) / subdivisons)

            if radialu
                u1 = sqrt(max(u1, zero(T)))
                u2 = sqrt(max(u2, zero(T)))
            end
            if radialv
                v1 = sqrt(max(v1, zero(T)))
                v2 = sqrt(max(v2, zero(T)))
            end

            u1 = (umax - umin) * u1 + umin - expansion
            u2 = (umax - umin) * u2 + umin + expansion
            v1 = (vmax - vmin) * v1 + vmin - expansion
            v2 = (vmax - vmin) * v2 + vmin + expansion

            p1 = point(surf, u1, v1)
            p2 = point(surf, u2, v1)
            p3 = point(surf, u2, v2)
            p4 = point(surf, u1, v2)

            if validtri(p1, p2, p3)
                push!(triangles, Triangle(p1, p2, p3, SVector{2,T}(u1, v1), SVector{2,T}(u2, v1), SVector{2,T}(u2, v2)))
            end
            if validtri(p1, p3, p4)
                push!(triangles, Triangle(p1, p3, p4, SVector{2,T}(u1, v1), SVector{2,T}(u2, v2), SVector{2,T}(u1, v2)))
            end
        end
    end

    return triangles
end
export triangulate

"""
    makemesh(object, subdivisions::Int = 30) -> TriangleMesh

Creates a [`TriangleMesh`](@ref) from an object, either a [`ParametricSurface`](@ref), [`CSGTree`](@ref) or certain surfaces (e.g. `Circle`, `Rectangle`).
This is used for visualization purposes only.
"""
function makemesh(surface::ParametricSurface{S,N}, subdivisions::Int = 30)::TriangleMesh{S} where {S,N}
    m = TriangleMesh(triangulate(surface, subdivisions, false))
    emptytrianglepool!(S)
    return m
end
export makemesh
