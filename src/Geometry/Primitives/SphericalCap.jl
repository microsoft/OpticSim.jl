# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

"""
    SphericalCap{T} <: ParametricSurface{T}

Spherical cap surface, creates a half-space which is essentially the subtraction of a sphere from an infinite plane.
Only the spherical cap itself is visualized, not the plane.
The positive normal side is outside the surface.

**Can be used as a detector in [`OpticalSystem`](@ref)s.**

```julia
SphericalCap(radius::T, ϕmax::T, [surfacenormal::SVector{3,T}, centrepoint::SVector{3,T}]; interface::NullOrFresnel{T} = nullinterface(T))
```

The minimal case returns a spherical cap centered at the origin with `surfacenormal = [0, 0, 1]`.
"""
struct SphericalCap{T} <: ParametricSurface{T,3}
    radius::T
    ϕmax::T
    zmax::T
    centrenormal::SVector{3,T}
    centrepoint::SVector{3,T}
    interface::NullOrFresnel{T}
    uvec::SVector{3,T}
    vvec::SVector{3,T}

    function SphericalCap(radius::T, ϕmax::T; interface::NullOrFresnel{T} = NullInterface(T)) where {T<:Real}
        @assert !isnan(radius)
        @assert radius > zero(T) && zero(T) < ϕmax < T(π)
        zmax = radius * (one(T) + sin(-(π / 2 - ϕmax)))
        new{T}(radius, ϕmax, zmax, SVector{3,T}(0.0, 0.0, 1.0), SVector{3,T}(0.0, 0.0, 0.0), interface, SVector{3,T}(1.0, 0.0, 0.0), SVector{3,T}(0.0, 1.0, 0.0))
    end

    function SphericalCap(radius::T, ϕmax::T, centrenormal::SVector{3,T}, centrepoint::SVector{3,T}; interface::NullOrFresnel{T} = NullInterface(T)) where {T<:Real}
        @assert radius > zero(T) && zero(T) < ϕmax < T(π)
        n̂ = normalize(centrenormal)
        rotationvec = SVector{3,T}(0.0, 1.0, 0.0)
        if abs(dot(rotationvec, n̂)) == one(T)
            rotationvec = SVector{3,T}(1.0, 0.0, 0.0)
        end
        uvec = normalize(cross(normalize(rotationvec), n̂))
        vvec = normalize(cross(n̂, uvec))
        zmax = radius * (one(T) + sin(-(π / 2 - ϕmax)))
        new{T}(radius, ϕmax, zmax, n̂, centrepoint, interface, uvec, vvec)
    end
end
export SphericalCap

Base.show(io::IO, a::SphericalCap{T}) where {T<:Real} = print(io, "SphericalCap{$T}($(a.radius), $(a.centrepoint), $(a.centrenormal), $(a.ϕmax), $(interface(a)))")

interface(a::SphericalCap{T}) where {T<:Real} = a.interface
radius(a::SphericalCap{T}) where {T<:Real} = a.radius
centroid(r::SphericalCap{T}) where {T<:Real} = r.centrepoint

uvrange(::Type{SphericalCap{T}}) where {T<:Real} = ((-T(π), T(π)), (zero(T), one(T))) # θ and ρ

function normal(r::SphericalCap{T}, θ::T, ρ::T) where {T<:Real}
    ϕ = π / 2 - ρ * r.ϕmax
    return -normalize(cos(-ϕ) * (cos(-θ) * r.uvec + sin(-θ) * r.vvec) + sin(-ϕ) * r.centrenormal)
end

function point(r::SphericalCap{T}, θ::T, ρ::T) where {T<:Real}
    ϕ = π / 2 - ρ * r.ϕmax
    return centroid(r) + radius(r) * (cos(-ϕ) * (cos(-θ) * r.uvec + sin(-θ) * r.vvec) + sin(-ϕ) * r.centrenormal + r.centrenormal)
end

uv(r::SphericalCap{T}, x::T, y::T, z::T) where {T<:Real} = uv(r, SVector{3,T}(x, y, z))
function uv(r::SphericalCap{T}, p::SVector{3,T}) where {T<:Real}
    prel = p - centroid(r) - radius(r) * r.centrenormal
    v = dot(prel, r.vvec) / radius(r)
    u = dot(prel, r.uvec) / radius(r)
    θ = -NaNsafeatan(v, u)
    z = clamp(dot(prel, -r.centrenormal) / radius(r), -one(T), one(T))
    ϕ = NaNsafeasin(z)
    ρ = (π / 2 - ϕ) / r.ϕmax
    return θ, ρ
end

function inside(r::SphericalCap{T}, p::SVector{3,T}) where {T<:Real}
    prel = (p - centroid(r) - radius(r) * r.centrenormal) / radius(r)
    v = dot(prel, r.vvec)
    u = dot(prel, r.uvec)
    if (u == zero(T) && v == zero(T)) || sqrt(u^2 + v^2) < sin(r.ϕmax)
        l = norm(prel)
        if l < one(T)
            return false
        end
    end
    return dot(r.centrenormal, p - centroid(r)) < r.zmax
end

function onsurface(r::SphericalCap{T}, p::SVector{3,T}) where {T<:Real}
    prel = (p - centroid(r) - radius(r) * r.centrenormal) / radius(r)
    v = dot(prel, r.vvec)
    u = dot(prel, r.uvec)
    if (u == zero(T) && v == zero(T)) || sqrt(u^2 + v^2) < sin(r.ϕmax)
        return norm(prel) == one(T) && dot(r.centrenormal, p - centroid(r)) < r.zmax
    else
        false
    end
end

function uvtopix(::SphericalCap{T}, uv::SVector{2,T}, imsize::Tuple{Int,Int}) where {T<:Real}
    θ, ρ = uv
    h, w = imsize
    u = (cos(θ) * ρ + one(T)) / 2
    v = (sin(θ) * ρ + one(T)) / 2
    pixu = Int(floor((w - 1) * u)) + 1
    pixv = Int(floor((h - 1) * v)) + 1
    return pixu, pixv
end

function surfaceintersection(sph::SphericalCap{T}, r::AbstractRay{T,3}) where {T<:Real}
    rad = radius(sph)
    orel = origin(r) - centroid(sph) - rad * sph.centrenormal

    ox, oy, oz = orel
    dx, dy, dz = direction(r)

    a = dx^2 + dy^2 + dz^2
    b = 2 * (ox * dx + oy * dy + oz * dz)
    c = (ox^2 + oy^2 + oz^2) - rad^2

    temp = quadraticroots(a, b, c)

    if temp === nothing
        if inside(sph, origin(r))
            return rayorigininterval(Infinity(T))
        else
            return EmptyInterval(T)
        end
    end

    t1, t2 = temp

    if isapprox(t1, t2, rtol = 1e-12, atol = 2 * eps(T))
        if inside(sph, origin(r))
            return rayorigininterval(Infinity(T))
        else
            return EmptyInterval(T)
        end
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
            if zero(T) <= ρ <= one(T)
                int1 = Intersection(t1, pt1, normal(sph, θ, ρ), θ, ρ, interface(sph))
            end
        end

        if pt2 !== nothing
            θ, ρ = uv(sph, pt2)
            if zero(T) <= ρ <= one(T)
                int2 = Intersection(t2, pt2, normal(sph, θ, ρ), θ, ρ, interface(sph))
            end
        end

        if int1 !== nothing && int2 !== nothing
            if t1 <= t2
                return Interval(int1, int2)
            else
                return Interval(int2, int1)
            end
        elseif int1 !== nothing
            if inside(sph, origin(r))
                return rayorigininterval(int1)
            else
                return positivehalfspace(int1)
            end
        elseif int2 !== nothing
            if inside(sph, origin(r))
                return rayorigininterval(int2)
            else
                return positivehalfspace(int2)
            end
        else
            if inside(sph, origin(r))
                return rayorigininterval(Infinity(T))
            else
                return EmptyInterval(T)
            end
        end
    end
end

BoundingBox(a::SphericalCap{T}) where {T<:Real} = BoundingBox(Plane(a.centrenormal, a.centrepoint + a.radius * a.centrenormal))
