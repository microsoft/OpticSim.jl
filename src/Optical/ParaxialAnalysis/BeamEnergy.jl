# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.


# struct LensCoordinateFrame{T<:Real}
#     lens::ParaxialLens{T}
#     coordinateframe::Geometry.Transform{T}

#     function LensCoordinateFrame(lens::ParaxialLens{T}, pixelposition::AbstractVector{T}) where{T}
#         optcenter = opticalcenter(lens)
#         beamcenterray = optcenter-pixelposition
#         n̂ = normalize(beamcenterray) #unit normal to the plane defined by the ray passing from the pixel through the optical center of the lens.
#         #create local coordinate frame where n̂ is the z axis
#        xaxis,yaxis = Geometry.get_orthogonal_vectors(n̂)
#         Geometry.Transform()
#         return new{T}(lens,Geometry.Transform(xaxis,yaxis,n̂,optcenter)) #need to make sure this is a right handed system)
#     end
# end
# export LensCoordinateFrame

# """projects points in world space onto the x-y plane of the local coordinate frame. Returns array of 2D vectors"""
# function project(coordinateframe::Geometry.Transform{T},points::Vector{V}) where{T<:Real,V<:AbstractVector}
#     result = coordinateframe .* points
#     return map((x)-> SVector{2,T}(x[1],x[2]),result)
# end
# export project


# """project points along the direction of the plane normal onto plane"""
# function project(plane::Plane{T,N},points::SVector{M,SVector{3,T}}) where{T<:Real,M,N}
#     result = MVector{M,SVector{3,T}}(undef)

#     for i in 1:M
#         dist = -distancefromplane(plane,points[1])
#         result[i] = points[i] + dist*normal(plane)
#     end

#     return SVector{M,SVector{3,T}}(result)
# end

"""Projects pupilpoints onto the lens plane and converts them to two dimensional points represented in the local x,y coordinates of the lens coordinate frame"""
function project(lens::ParaxialLens{T},displaypoint::SVector{3,T},pupilpoints::SVector{N,SVector{3,T}}) where{T<:Real,N}
    #need local transform for lens. Not quite sure how to organize this yet
    projectedpoints = MVector{N,SVector{2,T}}(undef)
    locpupil = map(x->lens.transform*x,pupilpoints) #transform pupil vertices into local coordinate frame of lens
    virtpoint = virtualpoint(lens, displaypoint) #compute virtual point corresponding to physical display point
    for (i,ppoint) in pairs(pupilpoints)
        vec = ppoint-virtpoint
        vecdist = distancefromplane(lens,ppoint)
        virtdist = distancefromplane(lens,virtpoint)
        scale =  virtdist/vecdist
        planepoint = scale*vec + virtpoint
        projectedpoints[i] = SVector{2,T}(planepoint[1],planepoint[2]) #local lens coordinate frame has z axis aligned with the positive normal to the lens plane.
    end
    return SVector{N,SVector{2,T}}(projectedpoints) #may need to make this a Vector{SVector} to be compatible with LazySets VPolygon constructors. Or maybe an MVector.
end

function intersection(lens::ParaxialLens{T},projectedpoints::SVector{N,SVector{2,T}}) where{T<:Real,N}
    lpoly = LazySets.VPolygon(vertices(lens))
    ppoly = LazySets.VPolygon(projectedpoints)
    return lpoly ∩ ppoly
end

"""Returns a number between 0 and 1 representing the ratio of the lens radiance to the pupil radiance. Assume lᵣ is the radiance transmitted through the lens from the display point. Some of this radiance, pᵣ, passes through the pupil. The beamenergy is the ratio pᵣ/lᵣ."""
function beamenergy(lens::ParaxialLens{T},displaypoint::AbstractVector{T},pupilpoints::AbstractVector{SVector{3,T}}) where{T<:Real}
    projectedpoints = project(lens,displaypoint,pupilpoints)
    virtpoint = virtualpoint(lens,displaypoint)
    beamlens = SphericalPolygon(vertices(lens),virtpoint,T(1)) #assumes lens vertices are represented in the local lens coordinate frame
    
    intsct = projectedpoints ∩ vertices(lens)
    if isempty(intsct)
        return T(0)
    else
        beamintsct = SphericalPolygon(intsct,virtpoint,T(1))
        return area(beamintsct)/area(beamlens)
    end
end

