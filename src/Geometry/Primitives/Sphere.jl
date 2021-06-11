# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

"""
    Sphere{T,N} <: ParametricSurface{T,N}

Spherical surface centered at the origin.

```julia
Sphere(radius::T = 1.0; interface::NullOrFresnel{T} = nullinterface(T))
```
"""
struct Sphere{T,N} <: ParametricSurface{T,N}
    radius::T
    surface_id::UUID
    interface::NullOrFresnel{T}

    function Sphere(radius::T; interface::NullOrFresnel{T} = NullInterface(T)) where {T<:Real}
        @assert !isnan(radius)
        @assert radius > zero(T)
        return new{T,3}(radius, uuid1(Random.GLOBAL_RNG), interface)
    end
end
export Sphere

surface_id(a::Sphere{T,N}) where {T<:Real,N} = a.surface_id
interface(a::Sphere{T,N}) where {T<:Real,N} = a.interface
radius(a::Sphere{T}) where {T<:Real} = a.radius

uvrange(::Type{Sphere{T,N}}) where {T<:Real,N} = ((-T(π), T(π)), (zero(T), T(π)))

onsurface(sph::Sphere{T,3}, x::T, y::T, z::T) where {T<:Real} = samepoint(x^2 + y^2 + z^2, radius(sph)^2)

inside(sph::Sphere{T,3}, x::T, y::T, z::T) where {T<:Real} = x^2 + y^2 + z^2 - radius(sph)^2 < zero(T)

point(sph::Sphere{T,3}, ϕ::T, θ::T) where {T<:Real} = SVector{3,T}(radius(sph) * sin(θ) * cos(-ϕ), radius(sph) * sin(θ) * sin(-ϕ), radius(sph) * cos(θ))

normal(::Sphere{T,3}, ϕ::T, θ::T) where {T<:Real} = SVector{3,T}(sin(θ) * cos(-ϕ), sin(θ) * sin(-ϕ), cos(θ))

partials(sph::Sphere{T,3}, ϕ::T, θ::T) where {T<:Real} = SVector{3,T}(radius(sph) * sin(θ) * -sin(-ϕ), radius(sph) * sin(θ) * cos(-ϕ), 0.0), SVector{3,T}(radius(sph) * cos(θ) * cos(-ϕ), radius(sph) * cos(θ) * sin(-ϕ), radius(sph) * -sin(θ))
function uv(::Sphere{T,3}, x::T, y::T, z::T) where {T<:Real}
    # avoid divide by zero for ForwardDiff
    ϕ = -NaNsafeatan(y, x)
    if x == zero(T) && y == zero(T)
        θ = zero(T)
    else
        θ = NaNsafeatan(sqrt(x^2 + y^2), z)
    end
    return SVector{2,T}(ϕ, θ)
end

# Assumes the ray has been transformed into the canonical sphere coordinate frame which has the vertical axis passing through (0,0,0) and aligned with the z axis.
function surfaceintersection(sph::Sphere{T,N}, r::AbstractRay{T,N}) where {T<:Real,N}
    ox, oy, oz = origin(r)
    dx, dy, dz = direction(r)
    rad = radius(sph)

    a = dx^2 + dy^2 + dz^2
    b = 2 * (ox * dx + oy * dy + oz * dz)
    c = (ox^2 + oy^2 + oz^2) - rad^2

    temp = quadraticroots(a, b, c)

    if temp === nothing
        return EmptyInterval(T) # no intersection with sphere and ray not contained entirely in sphere
    end

    t1, t2 = temp

    if isapprox(t1, t2, rtol = 1e-12, atol = 2 * eps(T))
        return EmptyInterval(T) # single root which indicates a tangent ray sphere intersection
    end

    if t1 > zero(T)
        pt1 = point(r, t1)
    else
        pt1 = nothing
    end

    if t2 > zero(T)
        pt2 = point(r, t2)
    else
        pt2 = nothing
    end

    let int1 = nothing, int2 = nothing
        if pt1 !== nothing
            θ, ρ = uv(sph, pt1)
            int1 = Intersection(t1, pt1, pt1, θ, ρ, surface_id(sph), interface(sph))
        end

        if pt2 !== nothing
            θ, ρ = uv(sph, pt2)
            int2 = Intersection(t2, pt2, pt2, θ, ρ, surface_id(sph), interface(sph))
        end

        if int1 !== nothing && int2 !== nothing
            if t1 <= t2
                return Interval(int1, int2)
            else
                return Interval(int2, int1)
            end
        elseif int1 !== nothing
            return rayorigininterval(int1)
        elseif int2 !== nothing
            return rayorigininterval(int2)
        else
            return EmptyInterval(T)
        end
    end
end

BoundingBox(sph::Sphere{T,3}) where {T<:Real} = BoundingBox(-radius(sph), radius(sph), -radius(sph), radius(sph), -radius(sph), radius(sph))
