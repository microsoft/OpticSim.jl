# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.


#display lenslet systems: emitter, paraxial lenslet, paraxial eye lens, retina detector.
#start with eyebox plane and hexagonal tiling. Project onto display surface, generate lenslets automatically.

using OpticSim.Geometry:Transform,world2local
using OpticSim:plane_from_points,surfaceintersection,closestintersection,Ray,Plane,ConvexPolygon,Sphere 
using OpticSim
using OpticSim.Repeat:tilevertices,HexBasis1,tilesinside
using Unitful:upreferred
using Unitful.DefaultSymbols:°

# using OpticSim.Vis:drawcells

"""compute the mean of the columns of `a`. If `a` is an `SMatrix` this is very fast and will not allocate."""
columncentroid(a::AbstractMatrix) = sum(eachcol(a))/size(a)[2] #works except for the case of zero dimensional matrix.
export columncentroid

function project(point::AbstractVector{T},projectionvector::AbstractVector{T},surface::Union{OpticSim.LeafNode{T,N}, S}) where {T,N,S<:ParametricSurface{T,N}}
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
function project(vertices::R, projectionvector::AbstractVector{T}, surface::Union{OpticSim.LeafNode{T,N}, S}) where {T,N,S<:ParametricSurface{T,N},R<:AbstractMatrix{T}}
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


"""Finds the best fit plane to `vertices` then projects `vertices` onto this plane by transforming from the global to the local coordinate frame. The projected points are represented in the local coordinate frame of the plane."""
function projectonbestfitplane(vertices::AbstractMatrix{T}) where{T}
    @assert size(vertices)[1] == 3 "projection only works for 3D points"

    center, _, localrotation  = plane_from_points(vertices) 
    toworld = Transform(localrotation,center) #compute local to world transformation
    tolocal = world2local(toworld)
    numpts = size(vertices)[2]
    result = tolocal * SMatrix{3,numpts}(vertices) 

    return result,toworld,tolocal #project onto best fit plane by setting z coordinate to zero.
end
export projectonbestfitplane

function testproject()
    normal = SVector(0.0, 0, 1)
    hex = HexBasis1()
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

"""projects convex polygon, represented by `vertices`, onto `surface` along vector `normal`. Assumes original polygon is convex and that the projection will be convex. No guarantee that this will be true but for smoothly curved surfaces that are not varying too quickly relative to the size of the polygon it should be true."""
function planarpoly(projectedpoints::AbstractMatrix{T}) where{T}
    planarpoints,toworld,_ = projectonbestfitplane(projectedpoints)
    numpts = size(planarpoints)[2]
    temp = SMatrix{2,numpts}(planarpoints[1:2,:])
    vecofpts = collect(reinterpret(reshape,SVector{2,T},temp))
    #this code is inefficient because of different data structures used by convex polygon and most other code for representing lists of vertices.
     return ConvexPolygon(toworld,vecofpts)
end

spherepoint(radius,θ,ϕ) = radius .* SVector(cos(θ)sin(ϕ),sin(θ),cos(θ)cos(ϕ))
export spherepoint

"""Computes points on the edges of the spherical rectangle defined by the range of θ,ϕ. This is used to determine lattice boundaries on the eyebox surface."""
function spherepoints(radius, θmin,θmax,ϕmin,ϕmax)
    θedges =  [spherepoint(radius,ϕ,θ) for θ in θmin:.5:θmax, ϕ in (ϕmin,ϕmax)]
    ϕedges =  [spherepoint(radius,ϕ,θ) for ϕ in ϕmin:.5:ϕmax, θ in (θmin,θmax)]
    allpoints = vcat(reshape(θedges,reduce(*,size(θedges))),reshape(ϕedges,reduce(*,size(ϕedges))))
    
    reshape(reinterpret(Float64,allpoints),3,length(allpoints)) #return points as 3xn matrix with points as columns
end
export spherepoints

"""given a total fov in θ  and ϕ compute sample points on the edges of the spherical rectangle."""
spherepoints(radius,θ,ϕ) = spherepoints(radius,-θ/2,θ/2,-ϕ/2,ϕ/2)

function testprojection()
    pts = spherepoints(1.0,-.2,-.2,1.0,1.1)
    surf = Plane(0.0,0.0,-1.0,0.0,0.0,0.0)
    dir = [0.0,0.0,-1.0]
    project(pts,dir,surf)
end
export testprojection

function bounds(pts::AbstractMatrix{T}) where{T} 
    return [extrema(row) for row in eachrow(pts)]
end
export bounds

"""projects points on the spher onto the eyebox along the -z axis direction."""
function eyeboxbounds(eyebox::OpticSim.Plane,dir::AbstractVector, radius,fovθ,fovϕ) 
    pts = spherepoints(radius,fovθ,fovϕ)
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

function spherepolygon(vertices,projectiondirection,sph::OpticSim.LeafNode)
    @assert typeof(vertices) <: SMatrix
    projectedpts = project(vertices,projectiondirection,sph)
    return planarpoly(projectedpts) #polygon on best fit plane to projectedpts
end

#TODO need to figure out what to use as the normal (eventually this will need to take into account the part of the eyebox the lenslet should cover) 
#TODO write test for planarpoly and generation of Paraxial lens system using planar hexagon as lenslet.
""" 
# Generate hexagonal polygons on a spherical display surface. 

Each hexagonal polygon corresponds to one lenslet. A hexagonal lattice is first generated in the eyebox plane and then the vertices of these hexagonal tiles are projected onto the spherical display surface. The projected points are not necessarily planar so they are projected onto a best fit plane. For systems with fields of view less that 90° the error is small.

`dir` is the 3d direction vector to project points from the eyebox plane to the spherical display surface. 

`radius` is the radius of the spherical display surface.

`fovθ,fovϕ` correspond to the field of view of the display as seen from the center of the eyebox plane.

`lattice` is the hexagonal lattice to tile the sphere with, HexBasis1 or HexBasis3, which are rotated versions of each other."""
function spherepolygons(eyebox::Plane{T,N},eyerelief,sphereradius,dir,fovθ,fovϕ,lattice) where{T,N}
    # if fovθ, fovϕ are in degrees convert to radians. If they are unitless then the assumption is that they represent radians
    eyeboxz = eyebox.pointonplane[3] 
    vertexofsphere = eyeboxz + ustrip(mm,eyerelief - sphereradius) #we don't use mm when creating shapes because Transform doesn't work properly with unitful values. Add the units back on here.
    sph = OpticSim.LeafNode(Sphere(ustrip(mm,sphereradius)),Geometry.translation(0.0,0.0,-vertexofsphere)) #Sphere objects are always centered at the origin so have to make 
    θ = upreferred(fovθ) #converts to radians if in degrees
    ϕ = upreferred(fovϕ) #converts to radians if in degrees
    tiles = eyeboxtiles(eyebox,dir,ustrip(mm,sphereradius),θ,ϕ,lattice)
    shapes = Vector{ConvexPolygon{6,T}}(undef,0)
    coordinates = Vector{Tuple{Int64,Int64}}(undef,0)
    for coords in eachcol(tiles)
        twodverts = tilevertices(coords,lattice)
        numverts = size(twodverts)[2]
        temp = MMatrix{1,numverts,T}(undef)
        fill!(temp,eyeboxz)
        zerorow = SMatrix{1,numverts,T}(temp)
        vertices = SMatrix{N,numverts,T}(vcat(twodverts,zerorow))
        @assert typeof(vertices) <: SMatrix
        push!(shapes,spherepolygon(vertices,-dir,sph))
        push!(coordinates,Tuple{T,T}(coords))
    end
    shapes,coordinates
end

function spherelenslets(eyeboxplane::Plane{T,N},eyerelief,focallength,dir,sphereradius,fovθ,fovϕ,lattice) where{T,N}
    lenspolys,tilecoords = spherepolygons(eyeboxplane,eyerelief, sphereradius, dir,fovθ,fovϕ,lattice)
    result = Vector{ParaxialLens{T}}(undef,length(lenspolys))
    empty!(result)
    for poly in lenspolys
        lenslet = ParaxialLensConvexPoly(focallength,poly,SVector{2,T}(T.((0,0))))
        push!(result,lenslet)
    end
    return result,tilecoords
end
export spherelenslets


# function testeyeboxtiles()
#     tiles = eyeboxtiles(Plane(0.0,0.0,1.0,0.0,0.0,12.0),[0.0,0.0,-1.0],30,deg2rad(55),deg2rad(45),HexBasis1())
#     drawcells(HexBasis1(),10,tiles)
#     tiles = eyeboxtiles(Plane(0.0,0.0,1.0,0.0,0.0,12.0),[0.0,0.0,-1.0],30,deg2rad(25),deg2rad(45),HexBasis1())
#     drawcells(HexBasis1(),10,tiles)
# end
# export testeyeboxtiles