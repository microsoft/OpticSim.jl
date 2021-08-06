# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

import Statistics

const _abs_err_orientation_2d = 2*eps(Float64)

"""
    ConvexPolygon{N, T<:Real} <: PlanarShape{T} 

General Convex Polygon surface, not a valid CSG object.
The rotation of the polygon around its normal is defined by `rotationvec`.
`rotationvec×surfacenormal` is taken as the vector along the u axis.

```julia
ConvexPolygon(local_frame::Transform{T}, local_polygon_points::Vector{SVector{2, T}}, interface::NullOrFresnel{T} = nullinterface(T))
```

The local frame defines the plane (spans by the right and up vectors) with the plane normal given by the forward vector.
the local_polygon_points are given with respect to the local frame and are 2D points.
NOTE: This class uses static vectors to hold the points which will lead to more efficient performance, but should not be used with polygons with more than 20-30 points.
"""
struct ConvexPolygon{N, T<:Real}  <: PlanarShape{T} 
    plane::Plane{T,3}
    local_frame::Transform{T}
    local_points::SMatrix{2,N,T}
    # for efficency
    _local_frame_inv::Transform{T}                                  # cache the inverse matrix to avoid computing it for every intersection test
    _local_lines::Vector{SVector{3, SVector{2, T}}}                 # defines the edge points + a third point representing the slopes in order to save some calculationsduring ray checking
    _length::Int64                                                  # cache the length of the polygon 

    function ConvexPolygon(
            local_frame::Transform{T},
            local_polygon_points::Vector{SVector{2, T}}, 
            interface::NullOrFresnel{T} = NullInterface(T)
        ) where {T<:Real}

        # need at least 3 points to define apolygon
        @assert length(local_polygon_points) > 3         

        local_center = Statistics.mean(local_polygon_points)
        world_center = local2world(local_frame) * Vec3(local_center[1], local_center[2], zero(T))
        
        # poly_points = [GeometricalPredicates.Point2D(p[1], p[2]) for p in local_polygon_points]
        # poly = GeometricalPredicates.Polygon2D(poly_points...)
        pts = local_polygon_points
        local_lines = SVector(
            [SVector(
                pts[i],                                             # A
                pts[mod(i, length(pts))+1],                         # B
                                                                    # bx = Bx - Ax,  by = By - Ay
                SVector(pts[mod(i, length(pts))+1][1] - pts[i][1], pts[mod(i, length(pts))+1][2] - pts[i][2]))
                for i in 1:length(pts)]...
        )

        plane = Plane(forward(local_frame), world_center, interface = interface)
        N = length(local_polygon_points)
        temp = MMatrix{2,N,T}(undef)
        for (i,pt) in pairs(local_polygon_points)
            temp[:,i] = pt
        end

        new{N, T}(plane, local_frame, SMatrix{2,N,T}(temp), inv(local_frame), local_lines, length(local_lines))
    end
end
export ConvexPolygon

centroid(poly::ConvexPolygon) = poly.plane.pointonplane

#function barrier to make vertices allocate less and be faster.
function to3d(pts::SMatrix{2,N,T,L}) where{N,L,T}
    temp = MMatrix{3,N,T}(undef)
    for row in 1:2
        for col in 1:N
            temp[row,col] = pts[row,col]
        end
    end

    for col in 1:N 
        temp[3,col] = T(0)
    end

    return SMatrix{3,N,T}(temp)
return temp
end

#this function allocates. Don't know why, it shouldn't but it does.
function vertices(poly::ConvexPolygon{N,T}) where{N,T<:Real}
   return poly.local_frame * to3d(poly.local_points)
end


function surfaceintersection(poly::ConvexPolygon{N,T}, r::AbstractRay{T,3}) where {N,T<:Real}
    interval = surfaceintersection(poly.plane, r)
    if interval isa EmptyInterval{T} || isinfiniteinterval(interval)
        return EmptyInterval(T) # no ray plane intersection or inside plane but no hit
    else
        intersect = halfspaceintersection(interval)
        p = point(intersect)

        local_p = poly._local_frame_inv * p
        
        @inline function orientation(l::SVector{3, SVector{2, T}})::Int8 where {T<:Real}
            cx = local_p[1] - l[1][1]                         # getx(p) - getx(geta(l))
            cy = local_p[2] - l[1][2]                         # gety(p) - gety(geta(l))
            _pr2 = -(l[3][1]*cy) + l[3][2]*cx           # _pr2 = -l._bx*cy + l._by*cx
            if _pr2 < -_abs_err_orientation_2d
                1
            elseif _pr2 > _abs_err_orientation_2d
                -1
            else
                0   #  can implement something more accurate in the future to minimize numerical errors
            end
        end

        @inline function inpolygon()
            side = orientation(poly._local_lines[1])
            for i = 2:poly._length
                orientation(poly._local_lines[i]) == side || return false   
            end
            return true
        end

        if !inpolygon()
            return EmptyInterval(T)
        else
            if dot(normal(poly), direction(r)) < zero(T)
                return positivehalfspace(intersect)
            else
                return rayorigininterval(intersect)
            end
        end
    end
end


"""
    makemesh(poly::ConvexPolygon{N, T}, ::Int = 0) where {N, T<:Real} -> TriangleMesh

Create a triangle mesh that can be rendered by iterating on the polygon's edges and for each edge use the centroid as the third vertex of the triangle.
"""
function makemesh(poly::ConvexPolygon{N,T}, ::Int = 0) where {N,T<:Real}
    c = centroid(poly)

    l2w = local2world(poly.local_frame)
    len = Size(poly.local_points)[2] 

    triangles = []
    for i in 1:len
        p1 = poly.local_points[:,i]
        p2 = poly.local_points[:,mod(i,len) + 1]

        tri = Triangle(
            Vector(l2w * Vec3(p2[1], p2[2], zero(T))), 
            Vector(c), 
            Vector(l2w * Vec3(p1[1], p1[2], zero(T))))
        push!(triangles, tri)
    end

    triangles = Vector{Triangle{T}}(triangles)

    return TriangleMesh(triangles)
end
