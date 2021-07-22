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
        return new{T}(lens,Geometry.Transform(xaxis,yaxis,n̂,optcenter)) #need to make sure this is a right handed system)
    end
end
export LensCoordinateFrame

"""projects points in world space onto the x-y plane of the local coordinate frame. Returns array of 2D vectors"""
function project(coordinateframe::Geometry.Transform{T},points::Vector{V}) where{T<:Real,V<:AbstractVector}
    result = coordinateframe .* points
    return map((x)-> SVector{2,T}(x[1],x[2]),result)
end
export project

function beam(lens::ParaxialLens,point::AbstractVector)
    
end

struct SphericalTriangle{T<:Real}
    ptvectors::SVector{3,SVector{3,T}}
    spherecenter::SVector{3,T}
    radius::T

    function SphericalTriangle(points::SVector{3,SVector{3,T}},spherecenter::SVector{3,T},radius::T) where{T<:Real}
        vecs = MVector{3,SVector{3,T}}(undef)

        for (i,point) in pairs(points)
            vecs[i] = normalize(point-spherecenter) #vector from the center of the sphere of radius 1 to a point on the sphere
        end
        return new{T}(SVector{3,SVector{3,T}}(vecs),spherecenter,radius)
    end
end
export SphericalTriangle

SphericalTriangle(points::Vector{Vector{T}},spherecenter::Vector{T},radius::T) where{T<:Real} = SphericalTriangle(SVector{3,SVector{3,T}}(points),SVector{3,T}(spherecenter),radius)

function area(tri::SphericalTriangle{T}) where{T<:Real}
    sum = T(0)
    sum += acos(tri.ptvectors[1]⋅tri.ptvectors[2])
    sum += acos(tri.ptvectors[2]⋅tri.ptvectors[3])
    sum += acos(tri.ptvectors[3]⋅tri.ptvectors[1])
    return (sum-π) * tri.radius^2
end


struct SphericalPolygon{T<:Real,N}
    ptvectors::SVector{N,SVector{3,T}}
    spherecenter::SVector{3,T}
    radius::T

    function SphericalPolygon(points::Vector{Vector{T}},spherecenter::Vector{T},radius::T) where{T<:Real}
        N = length(points)
        for (i,point) in pairs(points)
            println(point)
            points[i] = normalize(point-spherecenter)
        end
        return new{T,N}(SVector{N,SVector{3,T}}(points),SVector{3,T}(spherecenter),radius)
    end
end


"""Breaks the convex spherical polygon into spherical triangles and computes the sum of the angles of all the triangles. Because the edges from the centroid to the poly vertices are shared by 2 spherical triangles these angles are multiplied by 2. The sum of all the angles around the centroid is 2π. Have to subtract π for each of the N triangles so total polygon area is 2π -Nπ + 2∑(angle from centroid vector to vertex vector)."""
function area(poly::SphericalPolygon{T,N}) where{T<:Real,N}
    accum = T(0)
 
    centroid = normalize(sum(poly.ptvectors)) #point somewhere in the middle of the convex spherical polygon, assuming the polygon isn't a hemisphere or more.
    for i in 1:N
        accum += acos(poly.ptvectors[i]⋅centroid)
    end
    accum *= 2
    return (accum + 2π - N*π)*poly.radius^2
end

typetest(poly::SphericalPolygon{T,N}) where{T<:Real,N} = T(0)
export typetest

testdatapoly() = SphericalPolygon([[0.0,1.0,0.0],[1.0,0.0,0.0],[0.0,0.0,1.0]],[0.0,0.0,0.0],1.0)
export testdatapoly

function testarea()
    tri = SphericalTriangle([[0.0,1.0,0.0],[1.0,0.0,0.0],[0.0,0.0,1.0]],[0.0,0.0,0.0],1.0)
    poly =  testdatapoly()

    return area(tri), area(poly)
end
export testarea

function beamenergy()
end
