# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

"""
    Cylinder{T,N} <: ParametricSurface{T,N}

Cylinder of infinite height centered at the origin, oriented along the z-axis.
`visheight` is used for visualization purposes only, **note that this does not fully represent the surface**.

```julia
Cylinder(radius::T, visheight::T = 2.0; interface::NullOrFresnel{T} = nullinterface(T))
```
"""
struct Cylinder{T,N} <: ParametricSurface{T,N}
    radius::T
    visheight::T # Only used for visualization purposes.
    interface::NullOrFresnel{T}

    function Cylinder(radius::T, visheight::T = T(2.0); interface::NullOrFresnel{T} = NullInterface(T)) where {T<:Real}
        @assert radius > zero(T) && visheight > zero(T)
        @assert !isnan(radius)
        new{T,3}(radius, visheight, interface)
    end
end
export Cylinder

Base.show(io::IO, a::Cylinder{T}) where {T<:Real} = print(io, "Cylinder{$T}($(a.radius), $(interface(a)))")

interface(a::Cylinder{T}) where {T<:Real} = a.interface
radius(a::Cylinder{T}) where {T<:Real} = a.radius

uvrange(::Type{Cylinder{T,N}}) where {T<:Real,N} = ((-T(π), T(π)), (-one(T), one(T)))

onsurface(cyl::Cylinder{T,3}, x::T, y::T, ::T) where {T<:Real} = samepoint(x^2 + y^2, radius(cyl)^2)

inside(cyl::Cylinder{T,3}, x::T, y::T, ::T) where {T<:Real} = x^2 + y^2 - radius(cyl)^2 < zero(T)

point(cyl::Cylinder{T,3}, u::T, v::T) where {T<:Real} = SVector{3,T}(radius(cyl) * cos(u), radius(cyl) * sin(u), v * cyl.visheight / 2)

normal(::Cylinder{T,3}, u::T, ::T) where {T<:Real} = SVector{3,T}(cos(u), sin(u), 0.0)

partials(::Cylinder{T,3}, u::T, ::T) where {T<:Real} = SVector{3,T}(-sin(u), cos(u), 0.0), SVector{3,T}(0.0, 0.0, 1.0)
uv(cyl::Cylinder{T,3}, x::T, y::T, z::T) where {T<:Real} = SVector{2,T}(atan(y, x), 2 * z / cyl.visheight)

# Assumes the ray has been transformed into the canonical cylinder coordinate frame which has the cylinder axis passing through (0,0,0) and aligned with the z axis.
function surfaceintersection(cyl::Cylinder{T,N}, r::AbstractRay{T,N}) where {T<:Real,N}
    # in the cylinder coordinate frame the cylinder axis is (0,0,1). project ray into plane perpendicular to cylinder axis by dropping z coordinate
    ox, oy, oz = origin(r)
    dx, dy, dz = direction(r)
    rad = radius(cyl)
    c = (ox^2 + oy^2) - rad^2

    if samepoint(one(T), abs(dz)) # ray is parallel to the cylinder axis
        if c > zero(T)
            # if strictly outside the cylinder then no intersection
            return EmptyInterval(T)
        else
            # ray is contained entirely in the cylinder
            return rayorigininterval(Infinity(T))
        end
    end

    a = dx^2 + dy^2
    b = 2 * (ox * dx + oy * dy)

    temp = quadraticroots(a, b, c)

    if temp === nothing
        return EmptyInterval(T) # no intersection with cylinder and ray not contained entirely in cylinder
    end

    t1, t2 = temp

    if isapprox(t1, t2, rtol = 1e-12, atol = 2 * eps(T))
        return EmptyInterval(T) # single root which indicates a tangent ray cylinder intersection
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
            u, v = uv(cyl, pt1)
            int1 = Intersection(t1, pt1, SVector{3,T}(pt1[1], pt1[2], 0.0), u, v, interface(cyl))
        end

        if pt2 !== nothing
            u, v = uv(cyl, pt2)
            int2 = Intersection(t2, pt2, SVector{3,T}(pt2[1], pt2[2], 0.0), u, v, interface(cyl))
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

BoundingBox(cyl::Cylinder{T,3}) where {T<:Real} = BoundingBox(-radius(cyl), radius(cyl), -radius(cyl), radius(cyl), typemin(T), typemax(T))
