# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.


#display lenslet systems: emitter, paraxial lenslet, paraxial eye lens, retina detector.
#start with eyebox plane and hexagonal tiling. Project onto display surface, generate lenslets automatically.

using OpticSim.Geometry:Transform,world2local
using OpticSim:plane_from_points,surfaceintersection,closestintersection,Ray,Plane,ConvexPolygon,Sphere
using OpticSim.Repeat:tilevertices,HexBasis1,tilesinside

"""compute the mean of the columns of `a`. If `a` is an `SMatrix` this is very fast and will not allocate."""
centroid(a::AbstractMatrix) = sum(eachcol(a))/size(a)[2] #works except for the case of zero dimensional matrix.
export centroid

function project(point::AbstractVector{T},projectionvector::AbstractVector{T},surface::OpticSim.Surface{T}) where{T}
    ray = Ray(point, projectionvector)
    intsct = surfaceintersection(surface, ray)
    pointintsct = closestintersection(intsct,false)
    if pointintsct === nothing #point didn't project to a point on the surface
        return nothing
    else
        return OpticSim.point(pointintsct)
    end
end
export project

"""project the vertices of a polygon represented by `vertices` onto `surface` using the point as the origin and `projectionvector` as the projection direction. Return nothing if any of the projected points do not intersect the surface. The projected vertices are not guaranteed to be coplanar."""
function project(vertices::AbstractMatrix{T}, projectionvector::AbstractVector{T}, surface::OpticSim.Surface{T}) where {T}
    result = similar(vertices)

    for i in 1:size(vertices)[2]
        origin = vertices[:, i]
        pt = project(origin,projectionvector,surface)

        if pt === nothing #one of the polygon points didn't project to a point on the surface so reject the polygon
            return nothing
        else
            result[:, i] = pt
        end
    end
    return result
end

project(vertices::AbstractMatrix{T}, projectionvector::AbstractVector{T}) where{T} = project(SMatrix{size(vertices)...,T}(vertices),projectionvector)

"""Finds the best fit plane to `vertices` then projects `vertices` onto this plane by transforming from the global to the local coordinate frame. The projected points are represented in the local coordinate frame of the plane."""
function projectonbestfitplane(vertices::AbstractMatrix{T}) where{T}
    @assert size(vertices)[1] == 3 "projection only works for 3D points"

    center, _, localrotation  = plane_from_points(vertices) 
    toworld = Transform(localrotation,center) #compute local to world transformation
    tolocal = world2local(toworld)
    result = tolocal * vertices 
    return vcat(result[1:2,:],[0 for _ in 1:size(vertices)[2]]'),toworld,tolocal #project onto best fit plane by setting z coordinate to zero.
end
export projectonbestfitplane

function testproject()
    normal = SVector(0.0, 0, 1)
    hex = HexBasis1()
    # verts = SMatrix{3,6}(vcat(tilevertices((0, 0), hex), [0.0 0 0 0 0 0]))
    verts =vcat(tilevertices((0, 0), hex), [0.0 0 0 0 0 0])
    surf = Plane(normal, SVector(0.0, 0, 10))

    project(verts, normal, surf)
end
export testproject

function testprojectonplane()
    verts = tilevertices(HexBasis1())
    verts = vcat(verts,[0 for _ in 1:6]')
    projectonbestfitplane(SMatrix{size(verts)...}(verts))
end
export testprojectonplane

"""projects convex polygon, represented by `vertices`, onto `surface` along vector `normal`. Assumes original polygon is convex and that the projection will be convex. No guarantee that this will be true bur for smoothly curved surfaces that are not varying too quickly relative to the size of the polygon it should be true."""
function planarpoly(vertices,normal,surface)
    projectedpoints = project(vertices,normal,surface)
    planarpoints,toworld,tolocal = projectonbestfitplane(vertices)
    return ConvexPolygon(toworld,planarpoints[1:2,:])
end

spherepoint(radius,θ,ϕ) = radius .* SVector(cos(θ)sin(ϕ),sin(θ),cos(θ)cos(ϕ))
export spherepoint

"""Computes points on the edges of the spherical rectangle defined by the range of θ,ϕ. This is used to determine lattice boundaries on the eyebox surface."""
function spherepoints(radius, θmin,θmax,ϕmin,ϕmax)
    

    θedges =  [spherepoint(radius,ϕ,θ) for θ in θmin:.5:θmax, ϕ in (ϕmin,ϕmax)]
    ϕedges =  [spherepoint(radius,ϕ,θ) for ϕ in ϕmin:.5:ϕmax, θ in (θmin,θmax)]
    allpoints = vcat(reshape(θedges,reduce(*,size(θedges))),reshape(ϕedges,reduce(*,size(ϕedges))))
    # allpoints = vcat(reshape(θedges,reduce(*,size(θedges))))
    
    reshape(reinterpret(Float64,allpoints),3,length(allpoints)) #return points as 3xn matrix with points as columns
end
export spherepoints

"""given a total fov in θ  and ϕ compute sample points on the edges of the spherical rectangle."""
spherepoints(radius,θ,ϕ) = spherepoints(radius,-θ/2,θ/2,-ϕ/2,ϕ/2)

function testspherepoints()
    pts = spherepoints(1.0,-.2,-.2,1.0,1.1)
    surf = Plane(0.0,0.0,-1.0,0.0,0.0,0.0)
    dir = [0.0,0.0,-1.0]
    project(pts,dir,surf)
end
export testspherepoints

function bounds(pts::AbstractMatrix{T}) where{T} 
    println(pts)
    return [extrema(row) for row in eachrow(pts)]
end
export bounds

function eyeboxbounds(eyebox::OpticSim.Plane,dir::AbstractVector, radius,fovθ,fovϕ) 
    pts = spherepoints(radius,fovθ,fovϕ)
    display(pts)
    projectedpts = project(pts,dir,eyebox)
    return bounds(projectedpts)
end
export eyeboxbounds

function boxtiles(bbox,lattice)
    tiles = tilesinside(bbox[1][1],bbox[2][1],bbox[1][2],bbox[2][2],lattice)
end
export boxtiles

eyeboxtiles(eyebox,dir,radius,fovθ,fovϕ,lattice) = boxtiles(eyeboxbounds(eyebox,dir,radius,fovθ,fovϕ),lattice)
export eyeboxtiles