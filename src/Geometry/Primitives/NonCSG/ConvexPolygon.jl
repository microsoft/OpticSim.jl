# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

import Statistics
import GeometricalPredicates

"""
    ConvexPolygon{T} <: Surface{T}

General Convex Polygon surface, not a valid CSG object.
The rotation of the polygon around its normal is defined by `rotationvec`.
`rotationvec×surfacenormal` is taken as the vector along the u axis.

```julia
ConvexPolygon(local_frame::Transform{T}, local_polygon_points::Vector{Vector{T}}, interface::NullOrFresnel{T} = nullinterface(T))
```

The local frame defines the plane (spans by the right and up vectors) with the plane normal given by the forward vector.
the local_polygon_points are given with respect to the local frame and are 2D points.
"""
struct ConvexPolygon{T} <: Surface{T}
    plane::Plane{T,3}
    local_frame::Transform{T}
    local_points::Vector{Vector{T}}

    # for efficency
    _poly2d::GeometricalPredicates.Polygon2D{GeometricalPredicates.Point2D}

    function ConvexPolygon(
            local_frame::Transform{T},
            local_polygon_points::Vector{Vector{T}}, 
            interface::NullOrFresnel{T} = NullInterface(T)
        ) where {T<:Real}

        @assert length(local_polygon_points) > 3         # need at least 3 points to define apolygon
        @assert size(local_polygon_points[1])[1] == 2    # check that first point is a 2D point

        local_center = Statistics.mean(local_polygon_points)
        world_center = local2world(local_frame) * Vec3(local_center[1], local_center[2], zero(T))

        poly_points = [GeometricalPredicates.Point2D(p[1], p[2]) for p in local_polygon_points]
        poly = GeometricalPredicates.Polygon2D(poly_points...)

        plane = Plane(forward(local_frame), world_center, interface = interface)
        new{T}(plane, local_frame, local_polygon_points, poly)
    end
end
export ConvexPolygon

# Base.show(io::IO, poly::ConvexPolygon{T}) where {T<:Real} = print(io, "ConvexPolygon{$T}($(centroid(hex)), $(normal(hex)), $(hex.side_length), $(interface(hex)))")
centroid(poly::ConvexPolygon{T}) where {T<:Real} = poly.plane.pointonplane
interface(poly::ConvexPolygon{T}) where {T<:Real} = interface(poly.plane)
normal(poly::ConvexPolygon{T}) where {T<:Real} = normal(poly.plane)
normal(poly::ConvexPolygon{T}, ::T, ::T) where {T<:Real} = normal(poly)


function surfaceintersection(poly::ConvexPolygon{T}, r::AbstractRay{T,3}) where {T<:Real}
    interval = surfaceintersection(poly.plane, r)
    if interval isa EmptyInterval{T} || isinfiniteinterval(interval)
        return EmptyInterval(T) # no ray plane intersection or inside plane but no hit
    else
        intersect = halfspaceintersection(interval)
        p = point(intersect)

        local_p = world2local(poly.local_frame) * p
        @assert abs(local_p[3]) < 0.00001   # need to find a general epsilon - check that the point lies on the plain
        in_poly = GeometricalPredicates.inpolygon(poly._poly2d, GeometricalPredicates.Point(local_p[1], local_p[2]))

        if !in_poly
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


function makemesh(poly::ConvexPolygon{T}, ::Int = 0) where {T<:Real}
    c = centroid(poly)

    l2w = local2world(poly.local_frame)
    l = length(poly.local_points)

    triangles = []
    for i in 1:l
        p1 = poly.local_points[i]
        p2 = poly.local_points[mod(i,l) + 1]

        tri = Triangle(
            Vector(l2w * Vec3(p2[1], p2[2], zero(T))), 
            Vector(c), 
            Vector(l2w * Vec3(p1[1], p1[2], zero(T))))
        push!(triangles, tri)
    end

    triangles = Vector{Triangle{T}}(triangles)

    return TriangleMesh(triangles)
end
