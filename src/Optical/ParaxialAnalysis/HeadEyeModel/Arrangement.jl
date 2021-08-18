# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

abstract type ArrangementConvexShape{N, T<:Real} end


struct Shape{N, T} <: ArrangementConvexShape{N, T} 
    _center::SVector{N, T}
    _points::Vector{SVector{N, T}}
    _coordinates::Tuple{Int,Int}

    function Shape(center::SVector{N, T}, points::AbstractVector{SVector{N, T}}, coordinates::Tuple{Int,Int}) where {N, T<:Real}
        return new{N, T}(center, points, coordinates)
    end

    """convenience function to make Shape constructor more compatible with functions in Repeat module, which use SMatrix{N1,N2} to represent vertices rather than Vector{SVector}"""
    function Shape(center::SVector{N, T}, points::SMatrix{N,N2,T}, coordinates::Tuple{Int,Int}) where {N, N2, T<:Real}
        return new{N,T}(center, Vector{SVector{N,T}}([points[:,i] for i in 1:N2]),coordinates)
    end
end

center(h::Shape) = h._center
points(h::Shape) = h._points
coordinates(h::Shape) = h._coordinates

function get_shapes(basis::Basis{2,T},cells::AbstractMatrix{Int}, radius::T) where{T}
    basic_tile = Repeat.tilevertices(basis) * radius
    cols = Base.size(cells)[2]
    res = Vector{Shape{2, T}}(undef, 0)

get_shapes(type::Type{Any};) = @error "Unknown Type [$type] - available types are :hexagon, :rectangle"

function get_shapes(::Type{Hexagon}; resolution::Tuple{Int,Int}=(2,2), radius=1.5)::Vector{Shape{2, Float64}}
    cells = Repeat.hexcellsinbox(resolution[1],resolution[2]) 
    hexbasis = Repeat.HexBasis1()
    basic_tile = Repeat.tilevertices(hexbasis) * radius

    res = Vector{Shape{2, Float64}}(undef, 0)
    for c in cells
        center = SVector(hexbasis[c[1], c[2]]) * radius
        points = [(SVector(p...) + center) for p in eachrow(basic_tile)]
        center = center
        push!(res, Shape(center, points, (cell[1],cell[2])))
    end
    return res
end

get_shapes(type::Type{Any};) = @error "Unknown Type [$type] - available types are Hexagon, Rectangle"

function get_shapes(basis::HexBasis1{2,T};resolution::Tuple{Int,Int}=(2,2), size=1.5) where{T}
    println(typeof(Repeat.hexcellsinbox(resolution[1],resolution[2])))
    return get_shapes(basis,Repeat.hexcellsinbox(resolution[1],resolution[2]),size)
end

function get_shapes(basis::RectangularBasis{2,T};resolution::Tuple{Int,Int}=(2,2), size=1.5) where{T}
    temp = collect(Iterators.product(-resolution[1]:resolution[1], -resolution[2]:resolution[2]))
    cells = reshape(temp,length(temp))
    result = Matrix{Int64}(undef,2,length(temp))
    
    for i in 1:length(cells)
        result[1,i] = cells[i][1]
        result[2,i] = cells[i][2]
    end

    get_shapes(basis,result,size)
end

function project(shapes::Vector{Shape{2, T}}, csg::OpticSim.CSGTree{T})::Vector{Shape{3, T}} where {T<:Real}
    ray_origin_distance = -10.0
    projected_shapes = Vector{Shape{3, Float64}}(undef, 0)
    for shape in shapes
        ok = true               # will turn to false if one of the hexagon vertces is not hitting the CSGTree and we reject the hexagon
        projected_shape_points = Vector{SVector{3, Float64}}(undef, 0)
        for pt in points(shape)
            ray = Ray(Vec3(pt[1], pt[2], ray_origin_distance), Vec3(0.0, 0.0, 1.0))
    
            interval = surfaceintersection(csg, ray)
    
            if interval isa EmptyInterval || OpticSim.isinfiniteinterval(interval)
                ok = false;
                @warn "Hexagon Rejected - Coordinates $(coordinates(shape))"
                break;
            end
    
            p = point(OpticSim.lower(interval))
            push!(projected_shape_points, p)
        end
        if (ok)
            pt = center(shape)
            ray = Ray(Vec3(pt[1], pt[2], ray_origin_distance), Vec3(0.0, 0.0, 1.0))
            interval = surfaceintersection(csg, ray)
            p = point(OpticSim.lower(interval))
            projected_shape = Shape(p, projected_shape_points, coordinates(shape))
            push!(projected_shapes, projected_shape)
        end
    end
    return projected_shapes
end


function build_paraxial_lens(shape::Shape{3, T}; local_center_point = SVector(0.0, 0.0), parent_transform::Transform = identitytransform(), scale=1.0, focal_length=10.0)::ParaxialLens{T} where {T<:Real}
    pts = HeadEye.points(shape)
    cen = HeadEye.center(shape)

    # estimate a best fitting plane (least-squares wise)
    # switching to SMatrix format - should consider using SMatrix in the basic shapes instead of a vector of points
    pts_mat = reshape(SVector(reduce(vcat, pts)...), StaticArrays.Size(length(pts[1]),length(pts)))
    fitted_center, fitted_normal = HeadEye.plane_from_points2(pts_mat)    

    # create a local frame for the fitted plane
    local_frame = Transform(fitted_center, fitted_normal)

    # project the center on the estimated plane
    local_center = world2local(local_frame) * cen
    projected_center = local2world(local_frame) * Vec3(local_center[1], local_center[2], 0.0)

    # convert the local frame so the origin will align with the projected center - not a must but more elegant
    local_frame = Transform(projected_center, fitted_normal)

    # project the rest of the points and get their local coordinates in order to create the convex polygon
    local_points = Vector{SVector{2, Float64}}(undef, 0)
    for p in pts
        local_p = world2local(local_frame) * p * scale
        push!(local_points, SVector(local_p[1], local_p[2]))
    end

    # polygon = ConvexPolygon(local_frame, local_points)
    # push!(planar_polygons, polygon)

    lens = ParaxialLensConvexPoly(focal_length, parent_transform * local_frame, local_points, local_center_point)
    return lens
end

function build_emitter(lens::ParaxialLens{T}; parent_transform::Transform = identitytransform(), distance=0.05, size=(0.01, 0.01))::Emitters.Sources.AbstractSource{T} where {T<:Real}
    local_frame = Transform(Vec3(OpticSim.centroid(lens)), Vec3(normal(lens)))
    S = Emitters.Spectrum.Uniform()
    P = Emitters.AngularPower.Lambertian()
    O = Emitters.Origins.RectGrid(size[1], size[2], 2, 2)
    D = Emitters.Directions.Constant()
    local Tr = parent_transform * Transform(local2world(local_frame) * (unitZ3() * -distance), forward(local_frame) )
    source = Emitters.Sources.Source(Tr, S, O, D, P)    
    return source
end