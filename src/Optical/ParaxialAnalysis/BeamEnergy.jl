# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

#compute beam from pixel and paraxial lens
#project polygon onto plane perpendicular to ray through optical closestintersection
#compute beam energy passing through another polygonal aperture

struct LensCoordinateFrame{T<:Real}
    lens::ParaxialLens{T}
    coordinateframe::Geometry.Transform{T}

    function LensCoordinateFrame(lens::ParaxialLens{T}, pixelposition::AbstractVector{T}) where{T}
        optcenter = opticalcenter(lens)
        beamcenterray = optcenter-pixelposition
        n̂ = normalize(beamcenterray) #unit normal to the plane defined by the ray passing from the pixel through the optical center of the lens.
        #create local coordinate frame where n̂ is the z axis
       xaxis,yaxis = Geometry.get_orthogonal_vectors(n̂)
        Geometry.Transform()
        return new{T}(lens,Geometry.Transform(xaxis,yaxis,n̂,optcenter))
    end
end
export LensCoordinateFrame

"""projects points in world space onto the x-y plane of the local coordinate frame. Returns array of 2D vectors"""
function project(coordinateframe::Geometry.Transform{T},points::Vector{V}) where{T<:Real,V<:AbstractVector}
    result = coordinateframe .* points
    return map((x)-> SVector{2,T}(x[1],x[2]),result)
end
export project

function virtualimagedistance(lens::ParaxialLens, pixelposition::AbstractVector)
    fl = focallength(lens)
    normal = normal(a)

end

function beamenergy()
