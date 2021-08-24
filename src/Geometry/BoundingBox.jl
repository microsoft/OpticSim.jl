# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

"""
    BoundingBox{T<:Real}

Axis-aligned three-dimensional bounding box.

```julia
BoundingBox(xmin::T, xmax::T, ymin::T, ymax::T, zmin::T, zmax::T)
BoundingBox(s::Surface{T})
BoundingBox(s::ParametricSurface{T,3}, transform::Transform{T} = identitytransform(T))
BoundingBox(c::CSGTree{T})
BoundingBox(tri::Triangle{T})
BoundingBox(triangles::AbstractVector{Triangle{T}})
BoundingBox(points::AbstractArray{SVector{3,T}})
BoundingBox(la::LensAssembly{T})
```
"""
struct BoundingBox{T<:Real}
    xmin::T
    ymin::T
    zmin::T
    xmax::T
    ymax::T
    zmax::T

    function BoundingBox(xmin::T, xmax::T, ymin::T, ymax::T, zmin::T, zmax::T) where {T<:Real}
        # if anything is NaN then we want to fall back to the infinite bounding box
        # if the input bounds are infinite then we also want to fall back to the infinite bounding box - this is necessary because the signs of the Infinity can sometimes get messed up
        # isinf(-Inf) == true so this works for any infinity
        xmin = isnan(xmin) || isinf(xmin) ? typemin(T) : xmin
        ymin = isnan(ymin) || isinf(ymin) ? typemin(T) : ymin
        zmin = isnan(zmin) || isinf(zmin) ? typemin(T) : zmin
        xmax = isnan(xmax) || isinf(xmax) ? typemax(T) : xmax
        ymax = isnan(ymax) || isinf(ymax) ? typemax(T) : ymax
        zmax = isnan(zmax) || isinf(zmax) ? typemax(T) : zmax
        if xmin <= xmax && ymin <= ymax && zmin <= zmax
            return new{T}(xmin, ymin, zmin, xmax, ymax, zmax)
        else
            throw(ErrorException("Invalid bounding box"))
        end
    end
end
export BoundingBox

BoundingBox(object::Object) = BoundingBox(object.object)

function BoundingBox(s::ParametricSurface{T,3}, transform::Transform{T} = identitytransform(T)) where {T<:Real}
    # get the bounding box of a transformed bounding box
    bbox = BoundingBox(s)
    if transform == identitytransform(T)
        return bbox
    else
        p1 = transform * SVector(bbox.xmin, bbox.ymin, bbox.zmin)
        p2 = transform * SVector(bbox.xmin, bbox.ymax, bbox.zmin)
        p3 = transform * SVector(bbox.xmin, bbox.ymax, bbox.zmax)
        p4 = transform * SVector(bbox.xmin, bbox.ymin, bbox.zmax)
        p5 = transform * SVector(bbox.xmax, bbox.ymin, bbox.zmin)
        p6 = transform * SVector(bbox.xmax, bbox.ymax, bbox.zmin)
        p7 = transform * SVector(bbox.xmax, bbox.ymax, bbox.zmax)
        p8 = transform * SVector(bbox.xmax, bbox.ymin, bbox.zmax)
        return BoundingBox(SVector(p1, p2, p3, p4, p5, p6, p7, p8))
    end
end

function BoundingBox(tri::Triangle{T}) where {T<:Real}
    big = fill(typemin(T), SVector{3,T})
    small = fill(typemax(T), SVector{3,T})

    for vert in vertices(tri)
        small = min.(small, vert)
        big = max.(big, vert)
    end

    return BoundingBox(small[1], big[1], small[2], big[2], small[3], big[3])
end

function BoundingBox(triangles::AbstractVector{Triangle{T}}) where {T<:Real}
    big = fill(typemin(T), SVector{3,T})
    small = fill(typemax(T), SVector{3,T})

    for tri in triangles
        for vert in vertices(tri)
            small = min.(small, vert)
            big = max.(big, vert)
        end
    end

    return BoundingBox(small[1], big[1], small[2], big[2], small[3], big[3])
end

function BoundingBox(points::AbstractArray{SVector{3,T}}) where {T<:Real}
    xmax = maximum((x) -> x[1], points)
    ymax = maximum((x) -> x[2], points)
    zmax = maximum((x) -> x[3], points)
    xmin = minimum((x) -> x[1], points)
    ymin = minimum((x) -> x[2], points)
    zmin = minimum((x) -> x[3], points)
    return BoundingBox(xmin, xmax, ymin, ymax, zmin, zmax)
end

function area(a::BoundingBox{T}) where {T<:Real}
    dx = a.xmax - a.xmin
    dy = a.ymax - a.ymin
    dz = a.zmax - a.zmin
    return 2 * (dx * dy + dy * dz + dz * dx)
end

union(::Nothing, b::BoundingBox{T}) where {T<:Real} = b
union(a::BoundingBox{T}, ::Nothing) where {T<:Real} = a

function union(a::BoundingBox{T}, b::BoundingBox{T}) where {T<:Real}
    return BoundingBox(min(a.xmin, b.xmin), max(a.xmax, b.xmax), min(a.ymin, b.ymin), max(a.ymax, b.ymax), min(a.zmin, b.zmin), max(a.zmax, b.zmax))
end

function intersection(a::BoundingBox{T}, b::BoundingBox{T}) where {T<:Real}
    if a.xmax < b.xmin || a.ymax < b.ymin || a.zmax < b.zmin || b.xmax < a.xmin || b.ymax < a.ymin || b.zmax < b.zmin
        @info a, b
        return nothing
    else
        return BoundingBox(max(a.xmin, b.xmin), min(a.xmax, b.xmax), max(a.ymin, b.ymin), min(a.ymax, b.ymax), max(a.zmin, b.zmin), min(a.zmax, b.zmax))
    end
end

"""
    doesintersect(bbox::BoundingBox{T}, r::AbstractRay{T,3}) -> Bool

Tests whether `r` intersects an axis-aligned [`BoundingBox`](@ref), `bbox`.
"""
function doesintersect(a::BoundingBox{T}, r::AbstractRay{T,3}) where {T<:Real}
    if inside(a, origin(r)) || onsurface(a, origin(r))
        return true
    end

    d = direction(r)
    o = origin(r)

    # Infs and zeros get us into all kinds of problems with AutoDiff here..
    # work arounds are possible in all cases, the code just gets messy :(
    # really we are doing this:
    #     dfx = one(T) / d[1]
    #     dfy = one(T) / d[2]
    #     dfz = one(T) / d[3]
    #     t1 = (a.xmin - o[1]) * dfx
    #     t2 = (a.xmax - o[1]) * dfx
    #     t3 = (a.ymin - o[2]) * dfy
    #     t4 = (a.ymax - o[2]) * dfy
    #     t5 = (a.zmin - o[3]) * dfz
    #     t6 = (a.zmax - o[3]) * dfz
    # tmin = max(max(min(t1, t2), min(t3, t4)), min(t5, t6))
    # tmax = min(min(max(t1, t2), max(t3, t4)), max(t5, t6))

    tmin = typemin(T)
    tmax = typemax(T)

    if d[1] != zero(T)
        # don't really understand why the gradient fails here is xmin or xmax is Inf, but this fixes it
        tx1 = isinf(a.xmin) ? sign(d[1]) * a.xmin : (a.xmin - o[1]) / d[1]
        tx2 = isinf(a.xmax) ? sign(d[1]) * a.xmax : (a.xmax - o[1]) / d[1]
    else
        # avoiding divide by zero (when d[1] == 0) to preserve gradients
        tx1 = sign(a.xmin - o[1]) * typemax(T)
        tx2 = sign(a.xmax - o[1]) * typemax(T)
    end
    # avoid min/max on Infs to preserve gradients
    if tx1 > tx2
        if tmin < tx2
            tmin = tx2
        end
        if tmax > tx1
            tmax = tx1
        end
    else
        if tmin < tx1
            tmin = tx1
        end
        if tmax > tx2
            tmax = tx2
        end
    end

    # below uses the same NaN avoidance techniques...
    if d[2] != zero(T)
        ty1 = isinf(a.ymin) ? sign(d[2]) * a.ymin : (a.ymin - o[2]) / d[2]
        ty2 = isinf(a.ymax) ? sign(d[2]) * a.ymax : (a.ymax - o[2]) / d[2]
    else
        ty1 = sign(a.ymin - o[2]) * typemax(T)
        ty2 = sign(a.ymax - o[2]) * typemax(T)
    end
    if ty1 > ty2
        if tmin < ty2
            tmin = ty2
        end
        if tmax > ty1
            tmax = ty1
        end
    else
        if tmin < ty1
            tmin = ty1
        end
        if tmax > ty2
            tmax = ty2
        end
    end

    if d[3] != zero(T)
        tz1 = isinf(a.zmin) ? sign(d[3]) * a.zmin : (a.zmin - o[3]) / d[3]
        tz2 = isinf(a.zmax) ? sign(d[3]) * a.zmax : (a.zmax - o[3]) / d[3]
    else
        tz1 = sign(a.zmin - o[3]) * typemax(T)
        tz2 = sign(a.zmax - o[3]) * typemax(T)
    end
    if tz1 > tz2
        if tmin < tz2
            tmin = tz2
        end
        if tmax > tz1
            tmax = tz1
        end
    else
        if tmin < tz1
            tmin = tz1
        end
        if tmax > tz2
            tmax = tz2
        end
    end

    return tmax >= tmin && tmax > zero(T)
end
export doesintersect

function inside(a::BoundingBox{T}, p::SVector{3,T}) where {T<:Real}
    return a.xmin < p[1] < a.xmax && a.ymin < p[2] < a.ymax && a.zmin < p[3] < a.zmax
end

function onsurface(a::BoundingBox{T}, p::SVector{3,T}) where {T<:Real}
    return ((p[1] === a.xmin || p[1] === a.xmax) && a.ymin < p[2] < a.ymax && a.zmin < p[3] < a.zmax) || (a.xmin < p[1] < a.xmax && (p[2] === a.ymin || p[2] === a.ymax) && a.zmin < p[3] < a.zmax) || (a.xmin < p[1] < a.xmax && a.ymin < p[2] < a.ymax && (p[3] === a.zmin || p[3] === a.zmax))
end

"""
    surfaceintersection(bbox::BoundingBox{T}, r::AbstractRay{T,3}) -> Union{EmptyInterval{T},Interval{T}}

Calculates the intersection of `r` with an axis-aligned [`BoundingBox`](@ref), `bbox`.

Returns an [`EmptyInterval`](@ref) if there is no intersection or an [`Interval`](@ref) if there is one or two intersections.
Note that the uv of the returned intersection is always **0**.
"""
function surfaceintersection(a::BoundingBox{T}, r::AbstractRay{T,3}) where {T<:Real}
    d = direction(r)
    o = origin(r)

    # Infs and zeros get us into all kinds of problems with AutoDiff here..
    # work arounds are possible in all cases, the code just gets messy :(
    # really we are doing this:
    #     dfx = one(T) / d[1]
    #     dfy = one(T) / d[2]
    #     dfz = one(T) / d[3]
    #     t1 = (a.xmin - o[1]) * dfx
    #     t2 = (a.xmax - o[1]) * dfx
    #     t3 = (a.ymin - o[2]) * dfy
    #     t4 = (a.ymax - o[2]) * dfy
    #     t5 = (a.zmin - o[3]) * dfz
    #     t6 = (a.zmax - o[3]) * dfz
    # tmin = max(max(min(t1, t2), min(t3, t4)), min(t5, t6))
    # tmax = min(min(max(t1, t2), max(t3, t4)), max(t5, t6))

    tmin = typemin(T)
    tmax = typemax(T)
    upper_normal = nothing
    lower_normal = nothing

    if d[1] != zero(T)
        # don't really understand why the gradient fails here is xmin or xmax is Inf, but this fixes it
        tx1 = isinf(a.xmin) ? sign(d[1]) * a.xmin : (a.xmin - o[1]) / d[1]
        tx2 = isinf(a.xmax) ? sign(d[1]) * a.xmax : (a.xmax - o[1]) / d[1]
    else
        # avoiding divide by zero (when d[1] == 0) to preserve gradients
        tx1 = sign(a.xmin - o[1]) * typemax(T)
        tx2 = sign(a.xmax - o[1]) * typemax(T)
    end
    # avoid min/max on Infs to preserve gradients
    if tx1 > tx2
        if tmin < tx2
            lower_normal = SVector{3,T}(1, 0, 0)
            tmin = tx2
        end
        if tmax > tx1
            upper_normal = SVector{3,T}(-1, 0, 0)
            tmax = tx1
        end
    else
        if tmin < tx1
            lower_normal = SVector{3,T}(-1, 0, 0)
            tmin = tx1
        end
        if tmax > tx2
            upper_normal = SVector{3,T}(1, 0, 0)
            tmax = tx2
        end
    end

    # below uses the same NaN avoidance techniques...
    if d[2] != zero(T)
        ty1 = isinf(a.ymin) ? sign(d[2]) * a.ymin : (a.ymin - o[2]) / d[2]
        ty2 = isinf(a.ymax) ? sign(d[2]) * a.ymax : (a.ymax - o[2]) / d[2]
    else
        ty1 = sign(a.ymin - o[2]) * typemax(T)
        ty2 = sign(a.ymax - o[2]) * typemax(T)
    end
    if ty1 > ty2
        if tmin < ty2
            lower_normal = SVector{3,T}(0, 1, 0)
            tmin = ty2
        end
        if tmax > ty1
            upper_normal = SVector{3,T}(0, -1, 0)
            tmax = ty1
        end
    else
        if tmin < ty1
            lower_normal = SVector{3,T}(0, -1, 0)
            tmin = ty1
        end
        if tmax > ty2
            upper_normal = SVector{3,T}(0, 1, 0)
            tmax = ty2
        end
    end

    if d[3] != zero(T)
        tz1 = isinf(a.zmin) ? sign(d[3]) * a.zmin : (a.zmin - o[3]) / d[3]
        tz2 = isinf(a.zmax) ? sign(d[3]) * a.zmax : (a.zmax - o[3]) / d[3]
    else
        tz1 = sign(a.zmin - o[3]) * typemax(T)
        tz2 = sign(a.zmax - o[3]) * typemax(T)
    end
    if tz1 > tz2
        if tmin < tz2
            lower_normal = SVector{3,T}(0, 0, 1)
            tmin = tz2
        end
        if tmax > tz1
            upper_normal = SVector{3,T}(0, 0, -1)
            tmax = tz1
        end
    else
        if tmin < tz1
            lower_normal = SVector{3,T}(0, 0, -1)
            tmin = tz1
        end
        if tmax > tz2
            upper_normal = SVector{3,T}(0, 0, 1)
            tmax = tz2
        end
    end

    if !(tmax >= tmin && tmax > zero(T))
        return EmptyInterval(T)
    else
        if tmin <= zero(T)
            if tmax == typemax(T)
                return rayorigininterval(Infinity(T))
            else
                return rayorigininterval(Intersection(tmax, point(r, tmax), upper_normal, zero(T), zero(T), NullInterface(T)))
            end
        else
            lower = Intersection(tmin, point(r, tmin), lower_normal, zero(T), zero(T), NullInterface(T))
            if tmax == typemax(T)
                return positivehalfspace(lower)
            else
                return Interval(lower, Intersection(tmax, point(r, tmax), upper_normal, zero(T), zero(T), NullInterface(T)))
            end
        end
    end
end
