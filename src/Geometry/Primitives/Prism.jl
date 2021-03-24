# MIT License

# Copyright (c) Microsoft Corporation.

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE

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

