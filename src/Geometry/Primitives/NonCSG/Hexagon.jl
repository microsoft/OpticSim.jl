# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

"""
    Hexagon{T} <: Surface{T}

Hexagonal surface, not a valid CSG object.
The rotation of the hexagon around its normal is defined by `rotationvec`.
`rotationvec×surfacenormal` is taken as the vector along the u axis.

```julia
Hexagon(side_length::T, [surfacenormal::SVector{3,T}, centrepoint::SVector{3,T}]; rotationvec::SVector{3,T} = [0.0, 1.0, 0.0], interface::NullOrFresnel{T} = nullinterface(T))
```

The minimal case returns a rectangle centered at the origin with `surfacenormal = [0, 0, 1]`.
"""
struct Hexagon{T} <: PlanarShapes{T}
    plane::Plane{T,3}
    side_length::T
    uvec::SVector{3,T}
    vvec::SVector{3,T}

    function Hexagon(side_length::T; interface::NullOrFresnel{T} = NullInterface(T)) where {T<:Real}
        @assert side_length > zero(T)
        new{T}(Plane(zero(T), zero(T), one(T), zero(T), zero(T), zero(T), interface = interface), side_length, SVector{3,T}(1.0, 0.0, 0.0), SVector{3,T}(0.0, 1.0, 0.0))
    end

    function Hexagon(side_length::T, surfacenormal::AbstractArray{T,1}, centrepoint::AbstractArray{T,1}; interface::NullOrFresnel{T} = NullInterface(T), rotationvec::AbstractArray{T,1} = SVector{3,T}(0.0, 1.0, 0.0)) where {T<:Real}
        @assert length(surfacenormal) == 3 && length(centrepoint) == 3
        return Hexagon(side_length, SVector{3,T}(surfacenormal), SVector{3,T}(centrepoint), interface = interface, rotationvec = SVector{3,T}(rotationvec))
    end

    function Hexagon(side_length::T, surfacenormal::SVector{3,T}, centrepoint::SVector{3,T}; rotationvec::SVector{3,T} = SVector{3,T}(0.0, 1.0, 0.0), interface::NullOrFresnel{T} = NullInterface(T)) where {T<:Real}
        @assert side_length > zero(T)
        n̂ = normalize(surfacenormal)
        if abs(dot(rotationvec, n̂)) == one(T)
            rotationvec = SVector{3,T}(1.0, 0.0, 0.0)
        end
        uvec = normalize(cross(normalize(rotationvec), n̂))
        vvec = normalize(cross(n̂, uvec))
        new{T}(Plane(n̂, centrepoint, interface = interface), side_length, uvec, vvec)
    end
end
export Hexagon

Base.show(io::IO, hex::Hexagon{T}) where {T<:Real} = print(io, "Hexagon{$T}($(centroid(hex)), $(normal(hex)), $(hex.side_length), $(interface(hex)))")
centroid(hex::Hexagon{T}) where {T<:Real} = hex.plane.pointonplane
interface(hex::Hexagon{T}) where {T<:Real} = interface(hex.plane)
normal(hex::Hexagon{T}) where {T<:Real} = normal(hex.plane)

function surfaceintersection(hex::Hexagon{T}, r::AbstractRay{T,3}) where {T<:Real}
    interval = surfaceintersection(hex.plane, r)
    if interval isa EmptyInterval{T} || isinfiniteinterval(interval)
        return EmptyInterval(T) # no ray plane intersection or inside plane but no hit
    else
        intersect = halfspaceintersection(interval)
        p = point(intersect)
        q2u = abs(dot(p - centroid(hex), hex.uvec))
        q2v = abs(dot(p - centroid(hex), hex.vvec))
        h = hex.side_length * sqrt(3) / 2
        if q2u > h || q2v > hex.side_length
            return EmptyInterval(T)
        else
            if hex.side_length * h - hex.side_length / 2 * q2u - h * q2v >= 0
                if dot(normal(hex), direction(r)) < zero(T)
                    return positivehalfspace(intersect)
                else
                    return rayorigininterval(intersect)
                end
            else
                return EmptyInterval(T)
            end
        end
    end
end

function makemesh(hex::Hexagon{T}, ::Int = 0) where {T<:Real}
    uvec = hex.side_length * hex.uvec
    vvec = hex.side_length * hex.vvec
    c = centroid(hex)
    triangles = [Triangle(c + sin(0 * π / 3) * uvec + cos(0 * π / 3) * vvec, c, c + sin(1 * π / 3) * uvec + cos(1 * π / 3) * vvec), Triangle(c + sin(1 * π / 3) * uvec + cos(1 * π / 3) * vvec, c, c + sin(2 * π / 3) * uvec + cos(2 * π / 3) * vvec), Triangle(c + sin(2 * π / 3) * uvec + cos(2 * π / 3) * vvec, c, c + sin(3 * π / 3) * uvec + cos(3 * π / 3) * vvec), Triangle(c + sin(3 * π / 3) * uvec + cos(3 * π / 3) * vvec, c, c + sin(4 * π / 3) * uvec + cos(4 * π / 3) * vvec), Triangle(c + sin(4 * π / 3) * uvec + cos(4 * π / 3) * vvec, c, c + sin(5 * π / 3) * uvec + cos(5 * π / 3) * vvec), Triangle(c + sin(5 * π / 3) * uvec + cos(5 * π / 3) * vvec, c, c + sin(0 * π / 3) * uvec + cos(0 * π / 3) * vvec)]
    return TriangleMesh(triangles)
end
