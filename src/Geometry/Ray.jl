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

abstract type AbstractRay{T<:Real,N} end
export AbstractRay

"""
    Ray{T,N} <: AbstractRay{T,N}

Purely geometric ray, defined as `origin + alpha * direction`.

```julia
Ray(origin::SVector{N,T}, direction::SVector{N,T})
```

Has the following accessor methods:
```julia
direction(ray::Ray{T,N}) -> SVector{N,T}
origin(ray::Ray{T,N}) -> SVector{N,T}
```
"""
struct Ray{T,N} <: AbstractRay{T,N}
    origin::SVector{N,T}
    direction::SVector{N,T}

    function Ray(origin::SVector{N,T}, direction::SVector{N,T}) where {T<:Real,N}
        return new{T,N}(origin, normalize(direction))
    end

    function Ray(origin::AbstractArray{T,1}, direction::AbstractArray{T,1}) where {T<:Real}
        @assert length(origin) == length(direction)
        N = length(origin)
        return new{T,N}(SVector{N,T}(origin), normalize(SVector{N,T}(direction)))
    end
end
export Ray

direction(ray::Ray{T,N}) where {T<:Real,N} = ray.direction
origin(ray::Ray{T,N}) where {T<:Real,N} = ray.origin
export origin, direction

function Base.print(io::IO, a::Ray{T,N}) where {T,N}
    println(io, "$(rpad("Origin:", 20)) $(origin(a))")
    println(io, "$(rpad("Direction:", 20)) $(direction(a))")
end

"""
    point(ray::AbstractRay{T,N}, alpha::T) -> SVector{T, N}

Returns a point on the ray at origin + alpha * direction. Alpha must be >= 0.
"""
function point(ray::AbstractRay{T,N}, alpha::T) where {N,T<:Real}
    @assert alpha >= zero(T) "Alpha must be nonnegative. alpha value: $alpha"
    return origin(ray) + alpha * direction(ray)
end

"""
    closestpointonray(r::Ray{T,N}, point::SVector{N,T}) -> SVector{T,N

Returns the point on the ray closest to point.
"""
function closestpointonray(r::AbstractRay{T,N}, point::SVector{N,T}) where {T,N}
    #find t that gives the smallest distance between the ray and the point
    #ray = p0 + tr̂, distance = ||line - point||
    o = origin(r)
    r̂ = direction(r)
    t = dot(point .- o, r̂)

    if t >= 0.0
        return o .+ t .* r̂
    else
        return o
    end
end

"""
    distance(r::Ray{T,N}, point::SVector{N,T}) -> Union{T,Nothing}

Returns distance to the position on the ray closest to point. If t < 0 returns nothing.
"""
function distance(r::AbstractRay{T,N}, point::SVector{N,T}) where {T,N}
    #find t that gives the smallest distance between the ray and the point
    #ray = p0 + tr̂, distance = ||line - point||
    # this has redundant code with closestpointonray but it is faster to do this way than to call closestpointonray and then comput norm(result .- point)

    o = origin(r)
    r̂ = direction(r)
    t = dot(point .- o, r̂)

    if t >= 0.0
        temp = o .+ t .* r̂
        return norm(temp .- point)
    else
        return nothing
    end
end

"""
    α(ray::AbstractRay{T,N}, point::SVector{N,T}) -> T

Computes the alpha corresponding to the closest position on the ray to point
"""
α(ray::AbstractRay{T,N}, point::SVector{N,T}) where {T<:Real,N} = dot(direction(ray), (point .- origin(ray)))
