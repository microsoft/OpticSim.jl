# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

struct ConvexPolygon{T,N,M} <: ParametricSurface{T,N}
    plane::Plane{T,N}
    origin::SVector{N,T}
    vertices::SMatrix{N,M,T}
end

# uv(surface::ParametricSurface{T,N}, p::SVector{N,T}) -> SVector{2,T}
# uvrange(surface::ParametricSurface{T,N}) -> Tuple{Tuple{T,T},Tuple{T,T}}
# point(surface::ParametricSurface{T,N}, u::T, v::T) -> SVector{N,T}
# partials(surface::ParametricSurface{T,N}, u::T, v::T) -> Tuple{SVector{N,T}, SVector{N,T}}
# normal(surface::ParametricSurface{T,N}, u::T, v::T) -> SVector{N,T}
# inside(surface::ParametricSurface{T,N}, p: :SVector{N,T}) -> Bool
# onsurface(surface::ParametricSurface{T,N}, p::SVector{N,T}) -> Bool
# surfaceintersection(surface::ParametricSurface{T,N}, AbstractRay::Ray{T,N}) -> Union{EmptyInterval{T},Interval{T},DisjointUnion{T}}
# interface(surface::ParametricSurface{T,N}) -> OpticalInterface{T}

