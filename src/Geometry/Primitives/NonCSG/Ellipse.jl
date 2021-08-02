# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

"""
    Ellipse{T} <: Surface{T}

Elliptical surface, not a valid CSG object.
The rotation of the rectangle around its normal is defined by `rotationvec`.
`rotationvec×surfacenormal` is taken as the vector along the u axis.

**Can be used as a detector in [`AbstractOpticalSystem`](@ref)s.**

```julia
Ellipse(halfsizeu::T, halfsizev::T, [surfacenormal::SVector{3,T}, centrepoint::SVector{3,T}]; interface::NullOrFresnel{T} = nullinterface(T))
```

The minimal case returns an ellipse centered at the origin with `surfacenormal = [0, 0, 1]`.
"""
struct Ellipse{T} <: PlanarShapes{T}
    plane::Plane{T,3}
    halfsizeu::T
    halfsizev::T
    uvec::SVector{3,T}
    vvec::SVector{3,T}

    function Ellipse(halfsizeu::T, halfsizev::T; interface::NullOrFresnel{T} = NullInterface(T)) where {T<:Real}
        @assert halfsizeu > zero(T) && halfsizev > zero(T)
        new{T}(Plane(SVector{3,T}(0, 0, 1), SVector{3,T}(0, 0, 0), interface = interface), halfsizeu, halfsizev, SVector{3,T}(1.0, 0.0, 0.0), SVector{3,T}(0.0, 1.0, 0.0))
    end

    function Ellipse(halfsizeu::T, halfsizev::T, surfacenormal::SVector{3,T}, centrepoint::SVector{3,T}; interface::NullOrFresnel{T} = NullInterface(T), rotationvec::SVector{3,T} = SVector{3,T}(0.0, 1.0, 0.0)) where {T<:Real}
        @assert halfsizeu > zero(T) && halfsizev > zero(T)
        n̂ = normalize(surfacenormal)
        if abs(dot(rotationvec, n̂)) == one(T)
            rotationvec = SVector{3,T}(1.0, 0.0, 0.0)
        end
        uvec = normalize(cross(normalize(rotationvec), n̂))
        vvec = normalize(cross(n̂, uvec))
        new{T}(Plane(n̂, centrepoint, interface = interface), halfsizeu, halfsizev, uvec, vvec)
    end

    function Ellipse(plane::Plane{T,3}, halfsizeu::T, halfsizev::T, uvec::SVector{3,T}, vvec::SVector{3,T}) where {T<:Real}
        new{T}(plane, halfsizeu, halfsizev, uvec, vvec)
    end
end
export Ellipse

Base.show(io::IO, a::Ellipse{T}) where {T<:Real} = print(io, "Ellipse{$T}($(centroid(a)), $(normal(a)), $(a.halfsizeu), $(a.halfsizev), $(interface(a)))")



uvrange(::Type{Ellipse{T}}) where {T<:Real} = ((-T(π), T(π)), (zero(T), one(T))) # θ and ρ



point(r::Ellipse{T}, θ::T, ρ::T) where {T<:Real} = centroid(r) + ρ * (r.halfsizeu * cos(θ) * r.uvec + r.halfsizev * sin(θ) * r.vvec)
partials(r::Ellipse{T}, θ::T, ρ::T) where {T<:Real} = (ρ * (r.halfsizeu * -sin(θ) * r.uvec + r.halfsizev * cos(θ) * r.vvec), (r.halfsizeu * cos(θ) * r.uvec + r.halfsizev * sin(θ) * r.vvec))

uv(r::Ellipse{T}, x::T, y::T, z::T) where {T<:Real} = uv(r, SVector{3,T}(x, y, z))
function uv(r::Ellipse{T}, p::SVector{3,T}) where {T<:Real}
    v = dot(p - centroid(r), r.vvec)
    u = dot(p - centroid(r), r.uvec)
    θ = NaNsafeatan(v, u)
    rad = norm(r.halfsizeu * cos(θ) * r.uvec + r.halfsizev * sin(θ) * r.vvec)
    return SVector{2,T}(θ, norm(p - centroid(r)) / rad)
end

onsurface(a::Ellipse{T}, point::SVector{3,T}) where {T<:Real} = onsurface(a.plane, point) && zero(T) <= uv(a, point)[2] <= one(T)

function uvtopix(::Ellipse{T}, uv::SVector{2,T}, imsize::Tuple{Int,Int}) where {T<:Real}
    θ, ρ = uv
    h, w = imsize
    u = (cos(θ) * ρ + one(T)) / 2
    v = (sin(θ) * ρ + one(T)) / 2
    pixu = Int(floor((w - 1) * u)) + 1
    pixv = h - Int(floor((h - 1) * v))
    return pixu, pixv
end

centroid(r::Ellipse{T}) where {T<:Real} = r.plane.pointonplane

function surfaceintersection(ell::Ellipse{T}, r::AbstractRay{T,3}) where {T<:Real}
    interval = surfaceintersection(ell.plane, r)
    if interval isa EmptyInterval{T} || isinfiniteinterval(interval)
        return EmptyInterval(T) # no ray plane intersection or inside plane but no hit
    else
        intersect = halfspaceintersection(interval)
        p = point(intersect)
        θ, ρ = uv(ell, p)
        if ρ > one(T)
            return EmptyInterval(T) # no ray plane intersection
        else
            intuv = Intersection(α(intersect), p, normal(ell), θ, ρ, interface(ell))
            if dot(normal(ell), direction(r)) < zero(T)
                return positivehalfspace(intuv)
            else
                return rayorigininterval(intuv)
            end
        end
    end
end

vertices(e::Ellipse,subdivisions::Int = 10) = vertices3d(e,subdivisions)

function vertices3d(e::Ellipse{T},subdivisions::Int = 10) where{T}
    dθ = T(2π) / subdivisions
    centre = point(e, zero(T), zero(T))
    verts = MMatrix{3,subdivisions,T}(undef)
    for i in 0:(subdivisions - 1)
        θ1 = i * dθ - π
        verts[:,i+1] =  point(e, θ1, one(T))
    end
    return SMatrix{3,subdivisions,T}(verts)
end

function makemesh(c::Ellipse{T}, subdivisions::Int = 30) where {T<:Real}
    dθ = T(2π) / subdivisions
    centre = point(c, zero(T), zero(T))
    tris = Vector{Triangle{T}}(undef, subdivisions)
    for i in 0:(subdivisions - 1)
        θ1 = i * dθ - π
        θ2 = (i + 1) * dθ - π
        p1 = point(c, θ1, one(T))
        p2 = point(c, θ2, one(T))
        tris[i + 1] = Triangle(centre, p1, p2)
    end
    return TriangleMesh(tris)
end

###########################

function Circle(radius::T; interface::NullOrFresnel{T} = NullInterface(T)) where {T<:Real}
    return Ellipse(radius, radius; interface = interface)
end

"""
    Circle(radius, [surfacenormal, centrepoint]; interface = nullinterface(T))

Shortcut method to create a circle. The minimal case returns a circle centred at the origin with `normal = [0, 0, 1]`.
"""
function Circle(radius::T, surfacenormal::SVector{3,T}, centrepoint::SVector{3,T}; interface::NullOrFresnel{T} = NullInterface(T)) where {T<:Real}
    return Ellipse(radius, radius, surfacenormal, centrepoint; interface = interface, rotationvec = SVector{3,T}(0.0, 1.0, 0.0))
end

export Circle
