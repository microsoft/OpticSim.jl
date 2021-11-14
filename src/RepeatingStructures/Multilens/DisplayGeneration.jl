#display lenslet systems: emitter, paraxial lenslet, paraxial eye lens, retina detector.
#start with eyebox plane and hexagonal tiling. Project onto display surface, generate lenslets automatically.

using OpticSim.Geometry:Transform,world2local
using OpticSim:plane_from_points,surfaceintersection,closestintersection,Ray,Plane
using OpticSim.Repeat:tilevertices,HexBasis1

"""compute the mean of the columns of `a`. If `a` is an `SMatrix` this is very fast and will not allocate."""
centroid(a::AbstractMatrix) = sum(eachcol(a))/size(a)[2] #works except for the case of zero dimensional matrix.
export centroid

"""project the vertices of a polygon represented by `vertices` onto `surface` using the point as the origin and `projectionvector` as the projection direction. Return nothing if any of the pronected points do not intersect the surface. The projected vertices are not guaranteed to be coplanar."""
function project(vertices::SMatrix{3,N,T}, projectionvector::AbstractVector{T}, surface::OpticSim.Surface{T}) where {N,T}
    result = MMatrix{3,N,T,3*N}(undef)

    for i in 1:size(vertices)[2]
        origin = vertices[:, i]
        ray = Ray(origin, projectionvector)
        intsct = surfaceintersection(surface, ray)
        pointintsct = closestintersection(intsct,false)
        if pointintsct === nothing #one of the polygon points didn't project to a point on the surface so reject the polygon
            return nothing
        else
            result[:, i] = OpticSim.point(pointintsct)
        end
    end
    return SMatrix{3,N,T,3*N}(result)
end

"""Finds the best fit plane to `vertices` then projects `vertices` onto this plane by transforming from the global to the local coordinate frame. The projected points are represented in the local coordinate frame of the plane."""
function projectonplane(vertices::AbstractMatrix{T}) where{T}
    @assert size(vertices)[1] == 3 "projection only works for 3D points"

    center, normal, localrotation  = plane_from_points(vertices) 
    toworld = Transform(localrotation,center) #compute local to world transformation
    tolocal = world2local(toworld)
    result = tolocal * vertices 
    return vcat(result[1:2,:],[0 for _ in 1:size(vertices)[2]]'),toworld,tolocal #project onto plane by setting z coordinate to zero. Return local coordinate frame so other functions can convert local points to world points.
end
export projectonplane

function testproject()
    normal = SVector(0.0, 0, 1)
    hex = HexBasis1()
    verts = SMatrix{3,6}(vcat(tilevertices((0, 0), hex), [0.0 0 0 0 0 0]))
    surf = Plane(normal, SVector(0.0, 0, 10))

    project(verts, normal, surf)
end
export testproject

function testprojectonplane()
    verts = tilevertices(HexBasis1())
    verts = vcat(verts,[0 for _ in 1:6]')
    projectonplane(SMatrix{size(verts)...}(verts))
end
export testprojectonplane

