# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

"""
    ParaxialLens{T} <: Surface{T}

`surfacenormal` is the **output** direction of the lens.
Paraxial lens cannot act as the interface between two materials, hence only a single outside material is specified, by default Air.

Create with the following functions
```julia
ParaxialLensEllipse(focaldistance, halfsizeu, halfsizev, surfacenormal, centrepoint; rotationvec = [0.0, 1.0, 0.0], outsidematerial = OpticSim.GlassCat.Air, decenteruv = (0.0, 0.0))
ParaxialLensRect(focaldistance, halfsizeu, halfsizev, surfacenormal, centrepoint; rotationvec = [0.0, 1.0, 0.0], outsidematerial = OpticSim.GlassCat.Air, decenteruv = (0.0, 0.0))
ParaxialLensHex(focaldistance, side_length, surfacenormal, centrepoint; rotationvec = [0.0, 1.0, 0.0], outsidematerial = OpticSim.GlassCat.Air, decenteruv = (0.0, 0.0))
ParaxialLensConvexPoly(focaldistance, local_frame, local_polygon_points, local_center_point; outsidematerial = OpticSim.GlassCat.Air)
```
"""
struct ParaxialLens{T} <: Surface{T}
    shape::Union{Rectangle{T},Ellipse{T},Hexagon{T},ConvexPolygon{T}}
    interface::ParaxialInterface{T}

    function ParaxialLens(shape::Union{Rectangle{T},Ellipse{T},Hexagon{T},ConvexPolygon{T}}, interface::ParaxialInterface{T}) where {T<:Real}
        new{T}(shape, interface)
    end
end

function ParaxialLensRect(focaldistance::T, halfsizeu::T, halfsizev::T, surfacenormal::AbstractArray{T,1}, centrepoint::AbstractArray{T,1}; rotationvec::AbstractArray{T,1} = SVector{3,T}(0.0, 1.0, 0.0), outsidematerial::OpticSim.GlassCat.AbstractGlass = OpticSim.GlassCat.Air, decenteruv::Tuple{T,T} = (zero(T), zero(T))) where {T<:Real}
    @assert length(surfacenormal) == 3 && length(centrepoint) == 3
    return ParaxialLensRect(focaldistance, halfsizeu, halfsizev, SVector{3,T}(surfacenormal), SVector{3,T}(centrepoint), rotationvec = SVector{3,T}(rotationvec), outsidematerial = outsidematerial, decenteruv = decenteruv)
end

function ParaxialLensRect(focaldistance::T, halfsizeu::T, halfsizev::T, surfacenormal::SVector{3,T}, centrepoint::SVector{3,T}; rotationvec::SVector{3,T} = SVector{3,T}(0.0, 1.0, 0.0), outsidematerial::OpticSim.GlassCat.AbstractGlass = OpticSim.GlassCat.Air, decenteruv::Tuple{T,T} = (zero(T), zero(T))) where {T<:Real}
    r = Rectangle(halfsizeu, halfsizev, surfacenormal, centrepoint)
    centrepoint = centrepoint + decenteruv[1] * r.uvec + decenteruv[2] * r.vvec
    return ParaxialLens(r, ParaxialInterface(focaldistance, centrepoint, outsidematerial))
end

function ParaxialLensHex(focaldistance::T, side_length::T, surfacenormal::AbstractArray{T,1}, centrepoint::AbstractArray{T,1}; rotationvec::AbstractArray{T,1} = SVector{3,T}(0.0, 1.0, 0.0), outsidematerial::OpticSim.GlassCat.AbstractGlass = OpticSim.GlassCat.Air, decenteruv::Tuple{T,T} = (zero(T), zero(T))) where {T<:Real}
    @assert length(surfacenormal) == 3 && length(centrepoint) == 3
    return ParaxialLensHex(focaldistance, side_length, SVector{3,T}(surfacenormal), SVector{3,T}(centrepoint), rotationvec = SVector{3,T}(rotationvec), outsidematerial = outsidematerial, decenteruv = decenteruv)
end

function ParaxialLensHex(focaldistance::T, side_length::T, surfacenormal::SVector{3,T}, centrepoint::SVector{3,T}; rotationvec::SVector{3,T} = SVector{3,T}(0.0, 1.0, 0.0), outsidematerial::OpticSim.GlassCat.AbstractGlass = OpticSim.GlassCat.Air, decenteruv::Tuple{T,T} = (zero(T), zero(T))) where {T<:Real}
    h = Hexagon(side_length, surfacenormal, centrepoint)
    centrepoint = centrepoint + decenteruv[1] * h.uvec + decenteruv[2] * h.vvec
    return ParaxialLens(h, ParaxialInterface(focaldistance, centrepoint, outsidematerial))
end

function ParaxialLensConvexPoly(focaldistance::T, local_frame::Transform{T}, local_polygon_points::Vector{SVector{2, T}}, local_center_point::SVector{2, T}; outsidematerial::OpticSim.GlassCat.AbstractGlass = OpticSim.GlassCat.Air) where {T<:Real}
    poly = ConvexPolygon(local_frame, local_polygon_points)
    centrepoint = SVector{3, T}(local2world(local_frame) * Vec3(local_center_point[1], local_center_point[2], zero(T)))
    return ParaxialLens(poly, ParaxialInterface(focaldistance, centrepoint, outsidematerial))
end

function ParaxialLensEllipse(focaldistance::T, halfsizeu::T, halfsizev::T, surfacenormal::AbstractArray{T,1}, centrepoint::AbstractArray{T,1}; rotationvec::AbstractArray{T,1} = SVector{3,T}(0.0, 1.0, 0.0), outsidematerial::OpticSim.GlassCat.AbstractGlass = OpticSim.GlassCat.Air, decenteruv::Tuple{T,T} = (zero(T), zero(T))) where {T<:Real}
    @assert length(surfacenormal) == 3 && length(centrepoint) == 3
    return ParaxialLensEllipse(focaldistance, halfsizeu, halfsizev, SVector{3,T}(surfacenormal), SVector{3,T}(centrepoint), rotationvec = SVector{3,T}(rotationvec), outsidematerial = outsidematerial, decenteruv = decenteruv)
end

function ParaxialLensEllipse(focaldistance::T, halfsizeu::T, halfsizev::T, surfacenormal::SVector{3,T}, centrepoint::SVector{3,T}; rotationvec::SVector{3,T} = SVector{3,T}(0.0, 1.0, 0.0), outsidematerial::OpticSim.GlassCat.AbstractGlass = OpticSim.GlassCat.Air, decenteruv::Tuple{T,T} = (zero(T), zero(T))) where {T<:Real}
    e = Ellipse(halfsizeu, halfsizev, surfacenormal, centrepoint)
    centrepoint = centrepoint + decenteruv[1] * e.uvec + decenteruv[2] * e.vvec
    return ParaxialLens(e, ParaxialInterface(focaldistance, centrepoint, outsidematerial))
end

export ParaxialLens, ParaxialLensRect, ParaxialLensEllipse, ParaxialLensHex, ParaxialLensConvexPoly

interface(r::ParaxialLens{T}) where {T<:Real} = r.interface
centroid(r::ParaxialLens{T}) where {T<:Real} = centroid(r.shape)
normal(r::ParaxialLens{T}) where {T<:Real} = normal(r.shape)
normal(r::ParaxialLens{T}, ::T, ::T) where {T<:Real} = normal(r)
point(r::ParaxialLens{T}, u::T, v::T) where {T<:Real} = point(r.shape, u, v)
uv(r::ParaxialLens{T}, x::T, y::T, z::T) where {T<:Real} = uv(r, SVector{3,T}(x, y, z))
uv(r::ParaxialLens{T}, p::SVector{3,T}) where {T<:Real} = uv(r.shape, p)

function surfaceintersection(l::ParaxialLens{T}, r::AbstractRay{T,3}) where {T<:Real}
    itvl = surfaceintersection(l.shape, r)
    if itvl isa EmptyInterval{T}
        return EmptyInterval(T)
    else
        intsct = halfspaceintersection(itvl)
        u, v = uv(intsct)
        intsct = Intersection(α(intsct), point(intsct), normal(intsct), u, v, interface(l), flippednormal = flippednormal(intsct))
        return positivehalfspace(intsct)
    end
end

makemesh(l::ParaxialLens{T}, ::Int = 0) where {T<:Real} = makemesh(l.shape)

function processintersection(opticalinterface::ParaxialInterface{T}, point::SVector{N,T}, normal::SVector{N,T}, incidentray::OpticalRay{T,N}, temperature::T, pressure::T, ::Bool, firstray::Bool = false) where {T<:Real,N}
    raypower = power(incidentray)
    m = outsidematerialid(opticalinterface) # necessarily both the same
    n = one(T)
    if !isair(m)
        mat = glassforid(m)::OpticSim.GlassCat.Glass
        n = index(mat, wavelength(incidentray), temperature = temperature, pressure = pressure)::T
    end
    geometricpathlength = norm(point - origin(incidentray)) + (firstray ? zero(T) : T(RAY_OFFSET))
    thisraypathlength = n * geometricpathlength
    raypathlength = pathlength(incidentray) + thisraypathlength
    # refraction calculated directly for paraxial lens
    raydirection = normalize((opticalinterface.centroid - point) + direction(incidentray) * opticalinterface.focallength / dot(normal, direction(incidentray)))
    if opticalinterface.focallength < 0
        # flip for negative focal lengths
        raydirection = -raydirection
    end
    return raydirection, raypower, raypathlength
end
