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

"""
Each [`Interval`](@ref) consists of two `IntervalPoint`s, one of [`RayOrigin`](@ref), [`Intersection`](@ref) or [`Infinity`](@ref).
"""
abstract type IntervalPoint{T<:Real} end
Base.eltype(::IntervalPoint{T}) where {T<:Real} = T

abstract type FinitePoint{T} <: IntervalPoint{T} end
Base.eltype(::FinitePoint{T}) where {T<:Real} = T

"""
    Intersection{T,N} <: IntervalPoint{T}

Represents the point at which an [`Ray`](@ref) hits a [`Surface`](@ref).
This consists of the distance along the ray, the intersection point in world space, the normal in world space, the UV on the surface and the [`OpticalInterface`](@ref) hit.

Has the following accessor methods:
```julia
point(a::Intersection{T,N}) -> SVector{N,T}
normal(a::Intersection{T,N}) -> SVector{N,T}
uv(a::Intersection{T,N}) -> SVector{2,T}
u(a::Intersection{T,N}) -> T
v(a::Intersection{T,N}) -> T
α(a::Intersection{T,N}) -> T
interface(a::Intersection{T,N}) -> OpticalInterface{T}
flippednormal(a::Intersection{T,N}) -> Bool
```
"""
struct Intersection{T,N} <: FinitePoint{T}
    α::T # value of ray parameter at the point of intersection
    point::SVector{N,T}
    normal::SVector{N,T}
    u::T
    v::T
    interface::AllOpticalInterfaces{T} # returns a union of all OpticalInterface subtypes - can't have an abstract type here as it results in allocations
    flippednormal::Bool

    function Intersection(α::T, point::SVector{N,T}, normal::SVector{N,T}, u::T, v::T, interface::AllOpticalInterfaces{T}; flippednormal = false) where {T<:Real,N}
        new{T,N}(α, point, normalize(normal), u, v, interface, flippednormal)
    end

    function Intersection(α::T, point::AbstractVector{T}, normal::AbstractVector{T}, u::T, v::T, interface::AllOpticalInterfaces{T}; flippednormal = false) where {T<:Real}
        @assert length(point) == length(normal)
        N = length(point)
        new{T,N}(α, SVector{N,T}(point), SVector{N,T}(normal), u, v, interface, flippednormal)
    end
end
export Intersection

point(a::Intersection{T,N}) where {T<:Real,N} = a.point
normal(a::Intersection{T,N}) where {T<:Real,N} = a.normal
uv(a::Intersection{T,N}) where {T<:Real,N} = SVector{2,T}(a.u, a.v)
u(a::Intersection{T,N}) where {T<:Real,N} = a.u
v(a::Intersection{T,N}) where {T<:Real,N} = a.v
α(a::Intersection{T,N}) where {T<:Real,N} = a.α
interface(a::Intersection{T,N}) where {T<:Real,N} = a.interface
flippednormal(a::Intersection{T,N}) where {T<:Real,N} = a.flippednormal

function Base.print(io::IO, a::Intersection{T,N}) where {T<:Real,N}
    println(io, Intersection{T,N})
    println(io, "α \t$(α(a))")
    println(io, "Point \t$(point(a))")
    println(io, "Normal \t$(normal(a))")
    println(io, "normal is flipped? \t$(flippednormal(a))")
    println(io, "u \t$(u(a))")
    println(io, "v \t$(v(a))")
    println(io, "interface \t$(interface(a))")
end

"""
    reversenormal(a::Intersection{T,N})

Used by the CSG complement operator (i.e. [`csgdifference`](@ref)) to reverse the inside outside sense of the object.
"""
function reversenormal(a::Intersection{T,N}) where {T<:Real,N}
    return Intersection(α(a), point(a), -normal(a), u(a), v(a), interface(a), flippednormal = !flippednormal(a))
end

"""
    Infinity{T} <: IntervalPoint{T}

Point representing ∞ within an [`Interval`](@ref).

```julia
Infinity(T = Float64)
Infinity{T}()
```
"""
struct Infinity{T} <: IntervalPoint{T}
    Infinity(::Type{T} = Float64) where {T<:Real} = new{T}()
    Infinity{T}() where {T<:Real} = new{T}()
end
"""
    RayOrigin{T} <: IntervalPoint{T}

Point representing 0 within an [`Interval`](@ref), i.e. the start of the ray.

```julia
RayOrigin(T = Float64)
RayOrigin{T}()
```
"""
struct RayOrigin{T} <: FinitePoint{T}
    RayOrigin(::Type{T} = Float64) where {T<:Real} = new{T}()
    RayOrigin{T}() where {T<:Real} = new{T}()
end
export Infinity, RayOrigin

Base.:(==)(::RayOrigin{P}, ::RayOrigin{P}) where {P<:Real} = true
Base.:(==)(::Infinity{P}, ::Infinity{P}) where {P<:Real} = true
Base.:(==)(a::Intersection{P,N}, b::Intersection{P,N}) where {P<:Real,N} = α(a) == α(b)
Base.:(==)(::IntervalPoint{P}, ::IntervalPoint{P}) where {P<:Real} = false

Base.:(<=)(::RayOrigin{P}, ::RayOrigin{P}) where {P<:Real} = true
Base.:(<=)(::Infinity{P}, ::Infinity{P}) where {P<:Real} = true
Base.:(<=)(a::IntervalPoint{P}, b::Infinity{P}) where {P<:Real} = true
Base.:(<=)(a::Intersection{P,N}, b::Intersection{P,N}) where {P<:Real,N} = α(a) <= α(b)
Base.:(<=)(a::RayOrigin{P}, b::Infinity{P}) where {P<:Real} = true
Base.:(<=)(a::RayOrigin{P}, b::IntervalPoint{P}) where {P<:Real} = true
Base.:(<=)(::IntervalPoint{P}, ::IntervalPoint{P}) where {P<:Real} = false

Base.:(>=)(::RayOrigin{P}, ::RayOrigin{P}) where {P<:Real} = true
Base.:(>=)(::Infinity{P}, ::Infinity{P}) where {P<:Real} = true
Base.:(>=)(a::IntervalPoint{P}, b::Infinity{P}) where {P<:Real} = false
Base.:(>=)(a::Intersection{P,N}, b::Intersection{P,N}) where {P<:Real,N} = α(a) >= α(b)
Base.:(>=)(a::RayOrigin{P}, b::Infinity{P}) where {P<:Real} = false
Base.:(>=)(a::RayOrigin{P}, b::IntervalPoint{P}) where {P<:Real} = false
Base.:(>=)(::IntervalPoint{P}, ::IntervalPoint{P}) where {P<:Real} = false

Base.:(<)(::Infinity{P}, ::Infinity{P}) where {P<:Real} = false
Base.:(<)(::IntervalPoint{P}, ::Infinity{P}) where {P<:Real} = true
Base.:(<)(a::Intersection{P,N}, b::Intersection{P,N}) where {P<:Real,N} = α(a) < α(b)

Base.:(>)(::Infinity{P}, ::Infinity{P}) where {P<:Real} = false
Base.:(>)(::IntervalPoint{P}, ::Infinity{P}) where {P<:Real} = false
Base.:(>)(a::Intersection{P,N}, b::Intersection{P,N}) where {P<:Real,N} = α(a) > α(b)

Base.isless(a::IntervalPoint{P}, b::IntervalPoint{P}) where {P<:Real} = α(a) < α(b)

Base.eltype(::Intersection{T,N}) where {T<:Real,N} = T

"""
    isinfinity(a) -> Bool

Returns true if `a` is [`Infinity`](@ref). In performance critical contexts use `a isa Infinity{T}`.

"""
isinfinity(::Infinity) = true
isinfinity(::Any) = false

α(::Infinity{T}) where {T<:AbstractFloat} = typemax(T)
Base.eltype(::Infinity{T}) where {T<:Real} = T

"""
    israyorigin(a) -> Bool

Returns true if `a` is [`RayOrigin`](@ref). In performance critical contexts use `a isa RayOrigin{T}`.

"""
israyorigin(::Any) = false
israyorigin(::RayOrigin) = true

α(::RayOrigin{T}) where {T<:Real} = zero(T)
Base.eltype(::RayOrigin{T}) where {T<:Real} = T


"""
Apply a Transform to an Intersection object
"""
function Base.:*(a::Transform{T}, intsct::Intersection{T,3})::Intersection{T,3} where {T<:Real}
    u, v = uv(intsct)
    i = interface(intsct)
    if VERSION < v"1.6.0-DEV"
        # TODO REMOVE
        return @unionsplit OpticalInterface T i Intersection(α(intsct), a * point(intsct), Geometry.rotate(a, normal(intsct)), u, v, i, flippednormal = flippednormal(intsct))
    else
        return Intersection(α(intsct), a * point(intsct), Geometry.rotate(a, normal(intsct)), u, v, interface(intsct), flippednormal = flippednormal(intsct))
    end
end
