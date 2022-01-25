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
    shape::Union{Rectangle{T},Ellipse{T},Hexagon{T},ConvexPolygon{N, T} where {N}} 
    interface::ParaxialInterface{T}

    function ParaxialLens(shape::Union{Rectangle{T},Ellipse{T},Hexagon{T},ConvexPolygon{N, T} where {N}}, interface::ParaxialInterface{T}) where {T<:Real}
            new{T}(shape, interface)
    end
end

shape(a::ParaxialLens) = a.shape
export shape
opticalcenter(a::ParaxialLens) = opticalcenter(a.interface)
export opticalcenter
focallength(a::ParaxialLens) = focallength(a.interface)
export focallength
"""returns the 2 dimensional vertex points of the shape defining the lens aperture. These points lie in the plane of the shape"""
vertices(a::ParaxialLens) = vertices(a.shape)
"""All paraxial lens surfaces are planar"""
plane(a::ParaxialLens) = plane(a.shape)
struct VirtualPoint{T<:Real}
    center::SVector{3,T}
    direction::SVector{3,T}
    distance::T
 
    function VirtualPoint(center::AbstractVector{T},direction::AbstractVector{T},distance::T) where{T<:Real}
        @assert distance ≥ 0  #direction is encoded in the direction vector, but distance might have a sign. Need to discard it.
        new{T}(SVector{3,T}(center),SVector{3,T}(direction),distance)
    end

end

"""This will return (Inf,Inf,Inf) if the point is at infinity. In this case you probably should be using the direction of the VirtualPoint rather than its position"""
point(a::VirtualPoint) = a.center + a.direction*a.distance

distancefromplane(lens::ParaxialLens,point::AbstractVector) = distancefromplane(lens.shape,SVector{3}(point))
export distancefromplane

"""returns the virtual distance of the point from the lens plane. When |distance| == focallength then virtualdistance = ∞"""
function virtualdistance(focallength::T,distance::T) where{T<:Real}
    if abs(distance) == focallength
        return sign(distance)*T(Inf)
    else
        return distance*focallength/(focallength-abs(distance))
    end
end

"""computes the virtual point position corresponding to the input `point`, or returns nothing for points at infinity. `point` is specified in the lens coordinate frame"""
function virtualpoint(lens::ParaxialLens{T}, point::AbstractVector{T}) where{T}
    oc = opticalcenter(lens)
    point_oc = normalize(point - oc)
    distance = distancefromplane(lens,point)
    vdistance = virtualdistance(focallength(lens),distance)
    return VirtualPoint(oc,point_oc,abs(vdistance))
end
export virtualpoint

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

function ParaxialLensConvexPoly(focallength::T,convpoly::ConvexPolygon{N,T},local_center_point::SVector{2,T}; outsidematerial::OpticSim.GlassCat.AbstractGlass = OpticSim.GlassCat.Air) where {N, T<:Real}
    # the local frame information is stored in convpoly. The polygon points are 2D stored in the z = 0 local plane. The shape and vertices functions automatically apply local to world transformations so the points are represented in world space.
    centrepoint = SVector{3, T}(local2world(localframe(convpoly)) * Vec3(local_center_point[1], local_center_point[2], zero(T)))
    return ParaxialLens(convpoly, ParaxialInterface(focallength, centrepoint, outsidematerial))
end


function ParaxialLensConvexPoly(focaldistance::T, local_frame::Transform{T}, local_polygon_points::Vector{SVector{2, T}}, local_center_point::SVector{2, T}; outsidematerial::OpticSim.GlassCat.AbstractGlass = OpticSim.GlassCat.Air) where {N, T<:Real}
    # the local frame information is stored in convpoly. The polygon points are 2D stored in the z = 0 local plane. The shape and vertices functions automatically apply local to world transformations so the points are represented in world space.
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
point(r::ParaxialLens{T}, u::T, v::T) where {T<:Real} = point(r.shape, u, v)
uv(r::ParaxialLens{T}, x::T, y::T, z::T) where {T<:Real} = uv(r, SVector{3,T}(x, y, z))
uv(r::ParaxialLens{T}, p::SVector{3,T}) where {T<:Real} = uv(r.shape, p)

function surfaceintersection(l::ParaxialLens{T}, r::AbstractRay{T,3}) where {T<:Real}
    #this code seems unnecessary. If the ParaxialInterface is stored in the shape instead of the ParaxialLens then can do this: 
    #surfaceintersection(l::ParaxialLens....) = surfaceintersection(l.shape). Should be much faster than this code.
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
    #TODO there is an error in this formula for refraction. If the ray and the normal have a negative dot product the refracted ray will be the negative of the correct direction. Fix this.
    
    # refraction calculated directly for paraxial lens
    # raydirection = normalize((opticalinterface.centroid - point) + direction(incidentray) * opticalinterface.focallength / dot(normal, direction(incidentray)))
    #this works for both positive and negative focal lengths
    raydirection = normalize(sign(opticalinterface.focallength)*(opticalinterface.centroid - point) + direction(incidentray) * abs(opticalinterface.focallength / dot(normal, direction(incidentray))))
    
    # if opticalinterface.focallength < 0
    #     # flip for negative focal lengths
    #     raydirection = -raydirection
    # end
    return raydirection, raypower, raypathlength
end
