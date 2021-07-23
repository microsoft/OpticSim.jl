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


"""project points along the direction of the plane normal onto plane"""
function project(plane::Plane{T,N},points::SVector{M,SVector{3,T}}) where{T<:Real,M,N}
    result = MVector{M,SVector{3,T}}(undef)

    for i in 1:M
        dist = -distancefromplane(plane,points[1])
        result[i] = points[i] + dist*normal(plane)
    end

    return SVector{M,SVector{3,T}}(result)
end


"""Vector guaranteed to have unit length. Identical to SVector otherwise"""
struct UnitVector{N,T<:Real}
    vec::SVector{N,T}

    function UnitVector(a::AbstractVector{T}) where{T<:Real}
        N = length(a)
        return new{N,T}(normalize(SVector{N,T}(a)))
    end
end

vector(a::UnitVector) = a.vec

SVector(a::UnitVector) = a.vec

Base.:+(a::UnitVector,b::UnitVector) = vector(a) + vector(b)
Base.:+(a::SVector{3},b::UnitVector{3}) = a + vector(b)

"""returns the spherical angle formed by the cone with centervector at its center with neighbor1,neighbor2 the edges"""
function sphericalangle(centervector::UnitVector{3,T}, neighbor1::UnitVector{3,T}, neighbor2::UnitVector{3,T}) where{T<:Real}
    vec1 = normalize(neighbor1.vec - (neighbor1.vec⋅centervector.vec)*centervector.vec) #subtract off the component of the neighbor vectors from the center vector. This leaves only the component orthogonal to center vector
    vec2 = normalize(neighbor2.vec - (neighbor2.vec⋅centervector.vec)*centervector.vec)
    return acos(vec1⋅vec2)
end

struct SphericalTriangle{T<:Real}
    ptvectors::SVector{3,UnitVector{3,T}}
    spherecenter::SVector{3,T}
    radius::T

    function SphericalTriangle(points::SVector{3,SVector{3,T}},spherecenter::SVector{3,T},radius::T) where{T<:Real}
        vecs = MVector{3,UnitVector{3,T}}(undef)

        for (i,point) in pairs(points)
            vecs[i] = UnitVector(point-spherecenter) #vector from the center of the sphere of radius 1 to a point on the sphere
        end
        return new{T}(SVector{3,UnitVector{3,T}}(vecs),spherecenter,radius)
    end

    """Constructor which takes vectors which are guaranteed to be unitlength"""
    SphericalTriangle(vectors::SVector{3,UnitVector{3,T}},spherecenter::SVector{3,T},radius::T) where{T<:Real} = new{T}(vectors,spherecenter,radius)
        
end
export SphericalTriangle

SphericalTriangle(points::Vector{Vector{T}},spherecenter::Vector{T},radius::T) where{T<:Real} = SphericalTriangle(SVector{3,SVector{3,T}}(points),SVector{3,T}(spherecenter),radius)

function area(tri::SphericalTriangle{T}) where{T<:Real}
    vec1,vec2,vec3 = tri.ptvectors[1],tri.ptvectors[2],tri.ptvectors[3]
    sum = T(0)
    sum += sphericalangle(vec1,vec3,vec2)
    sum += sphericalangle(vec2,vec1,vec3)
    sum += sphericalangle(vec3,vec2,vec1)
    return (sum-π) * tri.radius^2
end

struct SphericalPolygon{T<:Real,N}
    ptvectors::SVector{N,UnitVector{3,T}}
    spherecenter::SVector{3,T}
    radius::T

    function SphericalPolygon(points::SVector{N,SVector{3,T}},spherecenter::SVector{3,T},radius::T) where{T<:Real,N}
        normpoints = MVector{N,UnitVector{3,T}}(undef)

        for (i,point) in pairs(points)
            normpoints[i] = UnitVector(point-spherecenter)
        end
        return new{T,N}(SVector{N,UnitVector{3,T}}(normpoints),SVector{3,T}(spherecenter),radius)
    end
end


"""Conceptually breaks the convex spherical polygon into spherical triangles and computes the sum of the angles of all the triangles. The sum of all the angles around the centroid is 2π. Have to subtract π for each of the N triangles. Rather than compute the angles of triangles formed by taking edges from the centroid to each vertex, can instead just compute the internal angle of neighboring edges. Total polygon area is 2π -Nπ + ∑(interior angles)."""
function area(poly::SphericalPolygon{T,N}) where{T<:Real,N}
    accum = T(0)
    ptvecs = poly.ptvectors

    for i in 2:N-1  
        accum += sphericalangle(ptvecs[i],ptvecs[i-1],ptvecs[i+1])
    end
    # finish up first and last interior angles which have different indexing because of wraparound
    accum += sphericalangle(ptvecs[N],ptvecs[N-1],ptvecs[1])
    accum += sphericalangle(ptvecs[1],ptvecs[N],ptvecs[2])

    accum = (2π -N*π + accum)*poly.radius^2
    return accum
end

function circlepoly(nsides; offset = [0.0,0.0,1.0])
    step = 2π/nsides
    result = MVector{nsides,SVector{3,Float64}}(undef)
    for i in 0:nsides-1
        y,x = sincos(step*i)
        result[i+1] = [x,y,0.0] .+ offset
    end
    return SVector{nsides,SVector{3,Float64}}(result)
end
export circlepoly

testtri() = SphericalTriangle(SVector{3,SVector{3,Float64}}(SVector{3,Float64}(0.0,1.0,0.),SVector{3,Float64}(1.0,0.0,0.0),SVector{3,Float64}(0.0,0.0,1.0)),SVector{3,Float64}(0.0,0.0,0.0),1.0)
export testtri

testdatapoly() = SphericalPolygon(SVector{3,SVector{3,Float64}}(SVector{3,Float64}(0.0,1.0,0.),SVector{3,Float64}(1.0,0.0,0.0),SVector{3,Float64}(0.0,0.0,1.0)),SVector{3,Float64}(0.0,0.0,0.0),1.0)
export testdatapoly

foursidedpoly() = SphericalPolygon(SVector{4,SVector{3,Float64}}(SVector{3,Float64}(0.0,1.0,0.0),SVector{3,Float64}(1.0,1.0,-.9),SVector{3,Float64}(1.0,0.0,0.0),SVector{3,Float64}(0.0,0.0,1.0)),SVector{3,Float64}(0.0,0.0,0.0),1.0)

sphericalcircle(nsides = 10) = SphericalPolygon(circlepoly(nsides),SVector(0.0,0.0,0.0),1.0)
export sphericalcircle

function testarea()
    tri = SphericalTriangle(SVector{3,SVector{3,Float64}}(SVector{3,Float64}(0.0,1.0,0.),SVector{3,Float64}(1.0,0.0,0.0),SVector{3,Float64}(0.0,0.0,1.0)),SVector{3,Float64}(0.0,0.0,0.0),1.0)
    poly =  testdatapoly()

    return area(tri), area(poly), area(foursidedpoly())
end
export testarea

