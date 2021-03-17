# MIT License

# Copyright (c) Microsoft Corporation.

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE

"""
    AcceleratedParametricSurface{T,N,S} <: ParametricSurface{T,N}

Wrapper class for [`ParametricSurface`](@ref)s where analytical intersection isn't feasible (e.g. [`ZernikeSurface`](@ref), [`ChebyshevSurface`](@ref)).
The surface is instead triangulated and an iterative (newton raphson) process carried out to determine precise ray intersection points.
`S` is the type of the ParametricSurface being wrapped.

```julia
AcceleratedParametricSurface(surf::ParametricSurface{T,N}, numsamples::Int = 17; interface::NullOrFresnel{T} = nullinterface(T))
```
"""
struct AcceleratedParametricSurface{T,N,S<:ParametricSurface{T,N}} <: ParametricSurface{T,N}
    surface::S
    triangles::Vector{Triangle{T}}
    sidelengths::Vector{T}
    triangles_bbox::BoundingBox{T}
    interface::NullOrFresnel{T}

    function AcceleratedParametricSurface(surf::S, numsamples::Int = 17; interface::NullOrFresnel{T} = NullInterface(T)) where {T<:Real,N,S<:ParametricSurface{T,N}}
        a = AcceleratedParametricSurface(surf, triangulate(surf, numsamples, true, true), interface = interface)
        emptytrianglepool!(T)
        return a
    end

    function AcceleratedParametricSurface(surf::S, triangles::Vector{Triangle{T}}; interface::NullOrFresnel{T} = NullInterface(T)) where {T<:Real,N,S<:ParametricSurface{T,N}}
        sidelengths = Vector{T}(undef, length(triangles))
        @inbounds @simd for i in 1:length(triangles)
            t = triangles[i]
            sidelengths[i] = max(norm(t.BA), norm(t.CA))
        end
        return new{T,N,S}(surf, copy(triangles), sidelengths, BoundingBox(triangles), interface)
    end
end
export AcceleratedParametricSurface


uvrange(::Type{AcceleratedParametricSurface{T,N,S}}) where {T<:Real,N,S<:ParametricSurface{T,N}} = uvrange(S)
uvrange(surf::AcceleratedParametricSurface{T,N,S}) where {T<:Real,N,S<:ParametricSurface{T,N}} = uvrange(surf.surface)
point(surf::AcceleratedParametricSurface{T,N,S}, u::T, v::T) where {T<:Real,N,S<:ParametricSurface{T,N}} = point(surf.surface, u, v)
partials(surf::AcceleratedParametricSurface{T,N,S}, u::T, v::T) where {T<:Real,N,S<:ParametricSurface{T,N}} = partials(surf.surface, u, v)
normal(surf::AcceleratedParametricSurface{T,N,S}, u::T, v::T) where {T<:Real,N,S<:ParametricSurface{T,N}} = normal(surf.surface, u, v)
inside(surf::AcceleratedParametricSurface{T,N,S}, x::T, y::T, z::T) where {T<:Real,N,S<:ParametricSurface{T,N}} = inside(surf.surface, x, y, z)
onsurface(surf::AcceleratedParametricSurface{T,N,S}, x::T, y::T, z::T) where {T<:Real,N,S<:ParametricSurface{T,N}} = onsurface(surf.surface, x, y, z)
"""
    interface(surf::Surface{T}) -> OpticalInterface{T}

Return the [`OpticalInterface`](@ref) associated with `surf`.
"""
interface(a::AcceleratedParametricSurface{T}) where {T<:Real} = a.interface
uv(surf::AcceleratedParametricSurface{T,3}, x::T, y::T, z::T) where {T<:Real} = uv(surf.surface, x, y, z)

function BoundingBox(surf::AcceleratedParametricSurface{T,3}) where {T<:Real}
    # the stored bbox is a tight one around the triangulated surface only, in reality the half-space created is an infinite prism capped by the triangulated surface
    # this bounding box is not accurate to the actual surface, but it is accurate to the triangulated surface
    # if a ray hits the accurate bounding box but not the triangulated bounding box then it would miss the triangulated surface
    # anyway (erroneously) and not return an intersection - so we don't lose anything by making this approximation here
    tbox = surf.triangles_bbox
    return BoundingBox(tbox.xmin, tbox.xmax, tbox.ymin, tbox.ymax, typemin(T), tbox.zmax)
end

# function makemesh(surface::AcceleratedParametricSurface{S,N,T}, ::Int) where {S<:Real,N,T<:ParametricSurface{S,N}}
#     # If we've already triangulated the surface then use this directly as it will be the most
#     # accurate representation of the surface being intersected with
#     return TriangleMesh(surface.triangles)
# end

# "computes all the intersections of the ray with the surface. Each intersection is represented as an Interval that encloses the interior of the object. If the ray is entering the surface, dot(ray,normal) < 0, then the interval is (intersection,infinity) where α is the parametric value of the line surface intersection. If the ray is leaving the surface, dot(ray,normal) > 0, then the interval is (rayorigin, α). If there are two intersections with the surface then returns an Interval(intersection,intersection). If there are multiple segments enclosing the object then will return a list that will contain some combination of the above types of Interval."
# fallback method for any accelerated surface which doesn't have an overridden method

"""
    surfaceintersection(surf::Surface{T}, r::AbstractRay{T}) where {T}

Calculates the intersection of `r` with a surface of any type, `surf`.
Note that some surfaces cannot be intersected analytically so must be wrapped in an [`AcceleratedParametricSurface`](@ref) in order to be intersected.

Returns an [`EmptyInterval`](@ref) if there is no [`Intersection`](@ref), an [`Interval`](@ref) if there is one or two intersections and a [`DisjointUnion`](@ref) if there are more than two intersections.
"""
function surfaceintersection(surf::AcceleratedParametricSurface{T,N}, r::AbstractRay{T,N}) where {T,N}
    if doesintersect(surf.triangles_bbox, r)
        i = triangulatedintersection(surf, r)
        if i isa EmptyInterval{T}
            if inside(surf, origin(r))
                return rayorigininterval(Infinity(T))
            else
                return EmptyInterval(T)
            end
        else
            return i
        end
    else
        if inside(surf, origin(r))
            return rayorigininterval(Infinity(T))
        else
            return EmptyInterval(T)
        end
    end
end

"""
    triangulatedintersection(surf::AcceleratedParametricSurface{T,N,S}, r::AbstractRay{T,N})

Intersection of a ray, `r`, with a triangulated surface, `surf`, no concept of inside so never returns a [`RayOrigin`](@ref) [`Interval`](@ref).
"""
function triangulatedintersection(surf::AcceleratedParametricSurface{T,N,S}, r::AbstractRay{T,N}) where {T<:Real,N,S}
    result = newinintervalpool!(T)
    offset_ray = Ray(origin(r) - ACCEL_SURF_RAY_OFFSET * direction(r), direction(r))
    for (i, tri) in enumerate(surf.triangles)
        # test if the ray is within the maximum side length from one of the points, if not then skip it
        w = vertex(tri, 1) - origin(r)
        dist = norm(w - dot(w, direction(r)) * direction(r))
        if dist > surf.sidelengths[i]
            continue
        end
        # find the triangle intersection - it's quite possible that for rays starting near the surface we will we miss
        # an intersection because of imprecise triangulation, so introduce an offset
        intv = surfaceintersection(tri, offset_ray)
        if !(intv isa EmptyInterval{T})
            ints = halfspaceintersection(intv)
            if α(ints) <= ACCEL_SURF_RAY_OFFSET
                # skip any intersection behind the true start point (with some wiggle room)
                continue
            end
            let dup = false
                # check if a point with the same alpha had already been included, if so ignore this one
                # doing this check now avoids executing newton many times unnecessarily
                for eintv in result
                    eints = halfspaceintersection(eintv)
                    if samepoint(α(eints), α(ints) + ACCEL_SURF_RAY_OFFSET)
                        dup = true
                        break
                    end
                end
                if !dup
                    # evaluate the precise point of the intersection with the surface
                    intuv = newton(surf.surface, r, uv(ints))
                    # check that the evaluated uv point lies on the (correct side of the) ray
                    if intuv !== nothing
                        u, v = intuv
                        surfpt = point(surf, u, v)
                        t = α(r, surfpt) # this is the 'true' α again (i.e. not offset)
                        if t > zero(T)
                            # check again that this point hasn't already been included (α might be slightly different, i.e. more accurate, now)
                            for eintv in result
                                eints = halfspaceintersection(eintv)
                                if samepoint(α(eints), t)
                                    dup = true
                                    break
                                end
                            end
                            if !dup && samepoint(point(r, t), surfpt)
                                # intersection is good so calc the normal and add to list
                                surfnormal = normal(surf, u, v)
                                intsct = Intersection(t, surfpt, surfnormal, u, v, interface(surf))
                                if dot(direction(r), surfnormal) < zero(T)
                                    intvl = positivehalfspace(intsct)
                                else
                                    intvl = rayorigininterval(intsct)
                                end
                                push!(result, intvl)
                            end
                        end
                    end
                end
            end
        end
    end

    nr = length(result)
    if nr === 0
        return EmptyInterval(T)
    elseif nr === 1
        return result[1]
    end

    # sort!(result, lt = (a, b) -> α(halfspaceintersection(a)) < α(halfspaceintersection(b))) #this will sort the intervals by α
    # alloc free sort of result
    @inbounds for i in 2:nr
        value = result[i]
        j = i - 1
        while j > 0 && α(halfspaceintersection(result[j])) > α(halfspaceintersection(value))
            result[j + 1] = result[j]
            j = j - 1
        end
        result[j + 1] = value
    end

    let start = 1
        temp = newinintervalpool!(T)

        if lower(result[1]) isa RayOrigin{T}
            # Starting "inside" the surface, i.e., n̂⋅d̂, > 0.
            # Can't truly define inside for Bezier patches because they don't define a half-space.
            start = 2
            push!(temp, result[1])
        end

        @inbounds for i in start:2:(nr - 1)
            intvlintsct = intervalintersection(result[i], result[i + 1])
            if !(intvlintsct isa EmptyInterval{T})
                push!(temp, intvlintsct)
            else
                # this can happen in some edge cases, seems like the best thing to do is juse ignore it...
                @warn "Funny behavior in triangulated intersection" maxlog = 1
                # throw(ErrorException("Triangulated intersection error - this should never happen."))
            end
        end

        if mod(nr - (start - 1), 2) === 1 #one unpaired interval
            push!(temp, last(result))
        end

        return DisjointUnion(temp)
    end
end

"""
    jacobian(surf::ParametricSurface{T,N}, u::T, v::T, P1::SVector{M,T}, P2::SVector{M,T})

Computes Jacobian of `f(t,u,v) = ( dot(P1,[surf(u,v),1],P2,[surf(u,v),1]) )`.
`P1`, `P2` are orthogonal planes that pass through the ray.
`J = [ ∂f1/∂u ∂f1/∂v ; ∂f2/∂u ∂f2/∂v]`
"""
function jacobian(surf::ParametricSurface{T,N}, u::T, v::T, P1::SVector{M,T}, P2::SVector{M,T}) where {T<:Real,N,M}
    (du, dv) = partials(surf, u, v)
    return SMatrix{2,2,T}(dot(view(P1, 1:(M - 1)), du), dot(view(P2, 1:(M - 1)), du), dot(view(P1, 1:(M - 1)), dv), dot(view(P2, 1:(M - 1)), dv))
end

"""
    newton(surf::ParametricSurface{T,N}, r::AbstractRay{T,N}, startingpoint::SVector{2,T})

Newton iteration to find the precise intersection of a parametric surface with a ray given a starting point (in uv space) on the surface.
"""
function newton(surf::ParametricSurface{T,N}, r::AbstractRay{T,N}, startingpoint::SVector{2,T}) where {T,N}
    tolerance = 1e-12
    maxiterations = 6

    urange, vrange = uvrange(surf)

    dx, dy, dz = direction(r)
    if abs(dx) > abs(dy) && abs(dx) > abs(dz)
        N1 = SVector{3,T}(dy, -dx, zero(T))
    else
        N1 = SVector{3,T}(zero(T), dz, -dy)
    end

    N2 = cross(N1, direction(r))
    P1 = SVector{4,T}(N1[1], N1[2], N1[3], -dot(N1, origin(r)))
    P2 = SVector{4,T}(N2[1], N2[2], N2[3], -dot(N2, origin(r)))

    @inline function f(uv::SVector{2,T})::SVector{2,T}
        psurf = point(surf, uv[1], uv[2])
        res1 = zero(T)
        res2 = zero(T)
        @inbounds for i in 1:3
            res1 += P1[i] * psurf[i]
            res2 += P2[i] * psurf[i]
        end
        res1 += P1[4]
        res2 += P2[4]
        # return [dot(P1, vcat(psurf, 1)), dot(P2, vcat(psurf, 1))]
        return SVector{2,T}(res1, res2)
    end

    xilast = startingpoint
    xi = zeros(SVector{2,T})
    error = typemax(T)
    i = 0
    finalpass = false

    while i < maxiterations
        if any(isnan.(xilast))
            # if something is NaN then just return the starting point and hope that the surface point matches
            xilast = startingpoint
            error = zero(T)
            break
        end

        j = jacobian(surf, xilast[1], xilast[2], P1, P2) # xi has values u,v

        if det(j) == zero(T)
            # if the jacobian is singular then break at the current iteration and return whatever we have
            error = zero(T)
            break
        end

        lstf = f(xilast)
        xi = xilast - (j \ lstf)
        error = sum(abs.(xi - xilast))

        if !finalpass && error < tolerance
            i = maxiterations - 3 # iterate a few more times after the xi and xi-1 are the same value, at least to 64 bit precision. Seems to give more accurate solutions. Perhaps because of higher precision for intermediate results?
            finalpass = true
        end

        i += 1

        xilast = xi
    end

    if error > tolerance
        return nothing
    else
        if xi[1] < urange[1] || xi[1] > urange[2] || xi[2] < vrange[1] || xi[2] > vrange[2]
            return nothing
        else
            return xilast
        end
    end
end

# Only show the underlying Bezier surface control points, not the approximating triangle mesh which is generally huge
Base.show(io::IO, a::AcceleratedParametricSurface{S,N,T}) where {S<:Real,N,T<:ParametricSurface{S,N}} = println(io, "Accelerated$(string(a.surface))")
