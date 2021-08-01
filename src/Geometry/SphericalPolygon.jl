"""returns the spherical angle formed by the cone with centervector at its center with neighbor1,neighbor2 the edges"""
function sphericalangle(neighbor1::SVector{3,T}, centervector::SVector{3,T}, neighbor2::SVector{3,T}) where{T<:Real}
    vec1 = normalize(neighbor1 - (neighbor1⋅centervector)*centervector) #subtract off the component of the neighbor vectors from the center vector. This leaves only the component orthogonal to center vector
    vec2 = normalize(neighbor2 - (neighbor2⋅centervector)*centervector)
    return acos(vec1⋅vec2)
end

struct SphericalTriangle{T<:Real}
    ptvectors::SMatrix{3,3,T}
    spherecenter::SVector{3,T}
    radius::T

    """Constructor which uses SMatrix to store vectors. For compatibility with LazySets this is more efficient"""
    function SphericalTriangle(vectors::SMatrix{3,3,T},spherecenter::SVector{3,T},radius::T) where{T} 
        temp = MMatrix{3,3,T}(undef)
        for col in 1:3
            temp[:,col] = normalize(vectors[:,col]-spherecenter)
        end

        new{T}(temp,spherecenter,radius)
    end
end
export SphericalTriangle

vector(a::SphericalTriangle,i::Int) = a.ptvectors[:,i]

function area(tri::SphericalTriangle{T}) where{T<:Real}
    vec1,vec2,vec3 = vector(tri,1),vector(tri,2),vector(tri,3)
    sum = T(0)
    sum += sphericalangle(vec1,vec3,vec2)
    sum += sphericalangle(vec2,vec1,vec3)
    sum += sphericalangle(vec3,vec2,vec1)
    return (sum-π) * tri.radius^2
end

struct SphericalPolygon{N,T<:Real}
    ptvectors::SMatrix{3,N,T}
    spherecenter::SVector{3,T}
    radius::T

    function SphericalPolygon(points::SMatrix{3,N,T},spherecenter::SVector{3,T},radius::T) where{T<:Real,N}
        normpoints = MMatrix{3,N,T}(undef)

        for i in 1:N
            normpoints[:,i] = normalize(points[:,i] - spherecenter)
        end
        return new{N,T}(SMatrix{3,N,T}(normpoints),SVector{3,T}(spherecenter),radius)
    end
end
export SphericalPolygon


"""Conceptually breaks the convex spherical polygon into spherical triangles and computes the sum of the angles of all the triangles. The sum of all the angles around the centroid is 2π. Have to subtract π for each of the N triangles. Rather than compute the angles of triangles formed by taking edges from the centroid to each vertex, can instead just compute the internal angle of neighboring edges. Total polygon area is 2π -Nπ + ∑(interior angles)."""
function area(poly::SphericalPolygon{N,T}) where{T<:Real,N}
    return area(poly.ptvectors,poly.radius)
end

function area(ptvecs::SMatrix{3,N,T},radius::T) where{N,T<:Real}
    accum = T(0)

    for i in 2:N-1 
        v1 = SVector{3,T}(ptvecs[1,i-1],ptvecs[2,i-1],ptvecs[3,i-1])      
        v2 = SVector{3,T}(ptvecs[1,i],ptvecs[2,i],ptvecs[3,i])      
        v3 = SVector{3,T}(ptvecs[1,i+1],ptvecs[2,i+1],ptvecs[2,i+1])      
        accum += sphericalangle(v1,v2,v3)
    end
    # finish up first and last interior angles which have different indexing because of wraparound
    v1 = SVector{3,T}(ptvecs[1,N-1],ptvecs[2,N-1],ptvecs[3,N-1])      
    v2 = SVector{3,T}(ptvecs[1,N],ptvecs[2,N],ptvecs[3,N])      
    v3 = SVector{3,T}(ptvecs[1,1],ptvecs[2,1],ptvecs[3,1])      

    accum += sphericalangle(v1,v2,v3)
    
    v1 = SVector{3,T}(ptvecs[1,2],ptvecs[2,2],ptvecs[3,2])      
    v2 = SVector{3,T}(ptvecs[1,1],ptvecs[2,1],ptvecs[3,1])      
    v3 = SVector{3,T}(ptvecs[1,N],ptvecs[2,N],ptvecs[3,N])      
   
    accum += sphericalangle(v1,v2,v3)

    accum = (2π -N*π + accum)*radius^2
    return accum
end
    
function circlepoly(nsides; offset = [0.0,0.0,1.0])
    step = 2π/nsides
    result = MMatrix{3,nsides,Float64}(undef)
    for i in 0:nsides-1
        y,x = sincos(step*i)
        result[:,i+1] = [x,y,0.0] .+ offset
    end
    return SMatrix{3,nsides,Float64}(result)
end
export circlepoly


"""creates a circular polygon that subtends a half angle of θ"""
function sphericalcircle(θ, nsides = 10)
    temp = MMatrix{3,nsides,Float64}(undef)
    for i in 0:1:(nsides-1)
        ϕ = i*2π/nsides
        temp[1,i+1] = sin(θ)*cos(ϕ)
        temp[2,i+1] = cos(θ)
        temp[3,i+1] = sin(θ)*sin(ϕ)
    end
    return SphericalPolygon(SMatrix{3,nsides,Float64}(temp),SVector(0.0,0.0,0.0),1.0)
end
export sphericalcircle


oneeigthsphere() = SphericalTriangle(SMatrix{3,3,Float64}(
    0.0,1.0,0.0,
    1.0,0.0,0.0,
    0.0,0.0,1.0),
    SVector(0.0,0.0,0.0),
    1.0)
export oneeigthsphere

onesixteenthphere() = SphericalTriangle(SMatrix{3,3,Float64}(
    0.0,1.0,0.0,
    1.0,1.0,0.0,
    0.0,0.0,1.0),
    SVector(0.0,0.0,0.0),
    1.0)
export onesixteenthphere

threesidedpoly() = SphericalPolygon(SMatrix{3,3,Float64}(
    0.0,1.0,0.0,
    1.0,0.0,0.0,
    0.0,0.0,1.0),
    SVector(0.0,0.0,0.0),
    1.0)
export threesidedpoly


foursidedpoly() = SphericalPolygon(SMatrix{3,4,Float64}(
    0.0,1.0,0.0,
    1.0,1.0,-.9,
    1.0,0.0,0.0,
    0.0,0.0,1.0),
    SVector(0.0,0.0,0.0),
    1.0)
export foursidedpoly

