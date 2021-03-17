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
    BSplineCurve{P,S,N,M} <: Spline{P,S,N,M}

`N` is the spatial dimension of the curve.
`M` is the curve order, i.e., the highest power of the parameterizing variable, `u`.
All curve segments are assumed to be of the same order.

```julia
BSplineCurve{P,S,N,M}(knots::KnotVector{S}, controlpoints::AbstractArray{MVector{N,S},1})
```
"""
struct BSplineCurve{P,S,N,M} <: Spline{P,S,N,M}
    knotvector::KnotVector{S}
    controlpolygon::Array{MVector{N,S},1} # may have to allow different orders on different curve segments

    function BSplineCurve{P,S,N,M}(knots::KnotVector{S}, controlpoints::AbstractArray{MVector{N,S},1}) where {P,S,N,M}
        bsplineinvariant(numknots(knots), length(controlpoints), M)
        newpoints = [MVector{N,S}(point) for point in controlpoints]
        new{P,S,N,M}(knots, newpoints)
    end
end
export BSplineCurve

numknots(curve::BSplineCurve) = numknots(curve.knotvector)

function findspan(curve::BSplineCurve{T,S,N,M}, u) where {T,S,N,M}
    findspan(curve.knotvector, M, u)
end

function numspans(curve::BSplineCurve{T,S,N,M}) where {T,S,N,M}
    return numknots(curve) - M - 1 # number of spans = numknots - curveorder - 1
end

spatialdimension(::Spline{T,S,N,M}) where {T,S,N,M} = N
curveorder(::Spline{T,S,N,M}) where {T,S,N,M} = M

function point(curve::BSplineCurve{P,S,N,M}, u::T)::SVector{N,S} where {T<:Real,P,S,N,M}
    # returns the raw point type of the curve - for Homogeneous curve types this will be a homogenous point
    span = findspan(curve, u)
    bases = basisfunctions(curve.knotvector, u, M)
    c = zeros(MVector{N,S})
    for i in 1:(M + 1)
        pt = curve.controlpolygon[span - (M + 1) + i]
        c .= c .+ bases[i] .* pt
    end
    return SVector{N,S}(c)
end

function euclideanpoint(curve::BSplineCurve{Euclidean,S,N,M}, ::T)::SVector where {T<:Real,S,N,M}
    # might seem weird to have this function for Euclidean curves. But in user code can use euclidean point without worrying about whether curve is homogeneous and still get correct results.
    return point(curve)
end

function euclideanpoint(curve::BSplineCurve{Rational,S,N,M}, u::T)::SVector where {T<:Real,S,N,M}
    return toeuclidean(point(curve, u))
end

"""
    BSplineSurface{P,S,N,M} <: SplineSurface{P,S,N,M}

Curve order is the same in the u and v direction and fixed over all spans.
u and v knot vectors are allowed to be different - _may change this to make them both the same_.

Control points in the u direction correspond to columns, with the lowest value of u corresponding to row 1.
Control points in the v direction correspond to rows, with the lowest value of v corresponding to col 1.

!!! danger
    This surface does not create a valid half-space, requires updates to function correctly.

```julia
BSplineSurface{P,S,N,M}(knots::KnotVector{S}, controlpoints::AbstractArray{<:AbstractArray{S,1},2})
BSplineSurface{P,S,N,M}(uknots::KnotVector{S}, vknots::KnotVector{S}, controlpoints::AbstractArray{<:AbstractArray{S,1},2})
```
"""
struct BSplineSurface{P,S,N,M} <: SplineSurface{P,S,N,M}
    uknotvector::KnotVector{S}
    vknotvector::KnotVector{S}
    controlpolygon::Array{MVector{N,S},2} # may have to allow different orders on different curve segments

    function BSplineSurface{P,S,N,M}(knots::KnotVector{S}, controlpoints::AbstractArray{<:AbstractArray{S,1},2}) where {P,S,N,M}
        return BSplineSurface(knots, knots, controlpoints)
    end

    function BSplineSurface{P,S,N,M}(uknots::KnotVector{S}, vknots::KnotVector{S}, controlpoints::AbstractArray{<:AbstractArray{S,1},2}) where {P,S,N,M}
        # m == n + M + 1
        # relation between number of knots = m + 1, number of control points = n + 1, and curve order M
        n_u, n_v = size(controlpoints)
        bsplineinvariant(numknots(uknots), n_u, M)
        bsplineinvariant(numknots(vknots), n_v, M)
        newpoints = [MVector{N,S}(point) for point in controlpoints]
        new{P,S,N,M}(uknots, vknots, newpoints)
    end
end
export BSplineSurface

numknots(surf::BSplineSurface) = (numknots(surf.uknotvector), numknots(surf.vknotvector))

uvrange(surface::BSplineSurface) = (urange(surface.uknotvector), urange(surface.vknotvector))

function findspan(surface::BSplineSurface{T,S,N,M}, u, v) where {T,S,N,M}
    return (findspan(surface.uknotvector, M, u), findspan(surface.vknotvector, M, v))
end

function insertknot(knots::AbstractArray{S,1}, knotindex::Int, controlpoints::Array{MVector{N,S},1}, curveorder::Int) where {N,S}
    knotvalue = knots[knotindex]
    M = curveorder
    K = [zeros(MVector{N,S}) for i in 1:(M + 1)]
    newknots = vcat(knots[1:knotindex], [knotvalue], knots[(knotindex + 1):end])
    newcontrolpoints = vcat(controlpoints[1:(knotindex - M)], K, controlpoints[(knotindex + 1):end])

    for i in (knotindex - M + 1):knotindex
        αᵢ = (knotvalue - knots[i]) / (knots[i + M] - knots[i])
        @. newcontrolpoints[i] = αᵢ * controlpoints[i] + (1 - αᵢ) * controlpoints[i - 1]
    end
    for i in (knotindex + 1):length(newcontrolpoints)
        newcontrolpoints[i] .= controlpoints[i - 1]
    end

    return (newknots, newcontrolpoints)
end

function insertknot(curve::BSplineCurve{T,S,N,M}, knotindex::Int) where {T,S,N,M}
    newknots, newcontrolpoints = insertknot(curve.knotvector.knots, knotindex, curve.controlpolygon, M)
    return BSplineCurve{T,S,N,M}(KnotVector{S}(newknots), newcontrolpoints)
end

function insertknots(curve::BSplineCurve{T,S,N,M}) where {T,S,N,M}
    newknots, newcontrolpoints = insertknots(curve.knotvector.knots, curve.controlpolygon, M)
    return BSplineCurve{T,S,N,M}(newknots, newcontrolpoints)
end

function insertknots(surf::BSplineSurface{T,S,N,M}) where {T,S,N,M}
    # TODO avoid using '
    newvknots, temppoints = expandcontrolpoints(surf.vknotvector.knots, surf.controlpolygon, M)
    newuknots, temppoints = expandcontrolpoints(surf.uknotvector.knots, copy(temppoints'), M)
    # TODO this adjoint is a pain for messing up all the types, ideally we wouldn't need to deepcopy
    # should have an arg for expandcontrolpoints that controls direction or something
    return BSplineSurface{T,S,N,M}(KnotVector{S}(newuknots), KnotVector{S}(newvknots), copy(temppoints'))
end

function expandcontrolpoints(vknots::AbstractArray{S,1}, controlpolygon::Array{MVector{N,S},2}, curveorder::Int) where {N,S}
    upts, vpts = size(controlpolygon)
    temp = Array{Array{MVector{N,S},1},1}(undef, 0) # TODO static array
    let newknots = vknots
        for i in 1:upts
            newknots, newcontrolpoints = insertknots(vknots, view(controlpolygon, i, :), curveorder)
            # not the most efficient since new knot arrays are being computed for reach row when all rows of knots are the same. Optimize later if necessary.
            # println("temp $temp newpts $newcontrolpoints")
            push!(temp, newcontrolpoints)
        end
        newcontrolpolygon = hcat(temp[:]...) # TODO more efficient way
        # make new 2D control point array
        return newknots, newcontrolpolygon
    end
end

function insertknots(knots::AbstractArray{S,1}, controlpoints::AbstractArray{MVector{N,S},1}, curveorder::Int) where {N,S}
    insertionknots = knotstoinsert(knots, curveorder)
    newknots = copy(knots)
    newcontrolpoints = Array{MVector{N,S},1}(controlpoints)
    offset = 0
    for insertion in insertionknots
        for i in 1:insertion[2]
            newknots, newcontrolpoints = insertknot(newknots, insertion[1] + offset, newcontrolpoints, curveorder)
        end
        offset += insertion[2]  # inserting knots changes the index where the next set of knots needs to be inserted.
    end
    return (newknots, newcontrolpoints)
end

function knotstoinsert(curve::Spline{T,S,N,M}) where {T,S,N,M}
    return knotstoinsert(curve.knotvector)
end

function knotstoinsert(knots::AbstractArray{S,1}, curveorder::Int) where {S}
    numknots = length(knots)
    index = curveorder + 2
    knotcounts = Array{Tuple{Int64,Int64},1}(undef, 0)

    while index < numknots - curveorder
        stop = index + 1
        while knots[index] == knots[stop]
            stop += 1       # stop will have the index after the last repeated knot in this sequence
        end
        numtoinsert = (curveorder) - (stop - index)
        if numtoinsert != 0
            push!(knotcounts, (index, numtoinsert))
        end
        # println("index $index stop $stop")
        index = stop
    end
    return knotcounts
end

function tobeziersegments(curve::BSplineCurve{T,S,N,M}) where {T,S,N,M}
    # returns an array of arrays of groups of M+1 Bezier control points, where M is the curve order.
    return tobeziersegments(curve.knotvector.knots, curve.controlpolygon, M)
end

function tobeziersegments(knots::AbstractArray{S,1}, controlpoints::Array{MVector{N,S},1}, curveorder::Int) where {N,S}
    # returns an array of arrays of groups of M+1 Bezier control points, where M is the curve order. Inefficient because it creates a new curve for every knot insertion. Optimize when necessary
    M = curveorder
    _, newcontrolpoints = insertknots(knots, controlpoints, curveorder)
    # each knot is now inserted M times, except for the first and last, which are repeated M+1 times.
    numsegments = length(knots) - 2 * (M + 1) + 1
    # println(newcurve.controlpolygon)
    # println("numsegments $numsegments")
    segments = Array{Array{MVector{N,S},1}}(undef, numsegments)

    for i in 1:numsegments
        start = (i - 1) * M + 1
        segments[i] = newcontrolpoints[start:(start + M)]
    end

    return segments
end

function tobeziersegments(surf::BSplineSurface{P,S,N,M}) where {P,S,N,M}
    newsurf = insertknots(surf)
    newcontrolpoints = newsurf.controlpolygon
    # each knot is now inserted M times, except for the first and last, which are repeated M+1 times.
    uknots, vknots = numknots(newsurf)
    numusegments = (uknots - 2 * (M + 1)) ÷ M + 1
    numvsegments = (vknots - 2 * (M + 1)) ÷ M + 1
    # println(newcurve.controlpolygon)
    # println("numsegments $numsegments")
    segments = Array{BezierSurface{P,S,N,M}}(undef, numusegments, numvsegments)
    # println("usges $numusegments vseg $numvsegments")
    for i in 1:numusegments
        for j in 1:numvsegments
            startu = (i - 1) * M + 1
            startv = (j - 1) * M + 1
            segments[i, j] = BezierSurface{P,S,N,M}(newcontrolpoints[startu:(startu + M), startv:(startv + M)])
        end
    end

    return segments
end

function point(surface::BSplineSurface{P,S,N,M}, u::T, v::T)::SVector{N,S} where {T<:Real,P,S,N,M}
    spanu, spanv = findspan(surface, u, v)

    ubases = basisfunctions(surface.uknotvector, u, M)
    vbases = basisfunctions(surface.vknotvector, v, M)
    sum = zeros(MVector{N,S})
    c = zeros(MVector{N,S})

    for j in 1:(M + 1)
        for i in 1:(M + 1)
            pt = surface.controlpolygon[spanu - (M + 1) + i, spanv - (M + 1) + j]
            c .= c .+ ubases[i] .* pt
        end
        sum .= sum .+ vbases[j] .* c
        c .= 0
    end

    return SVector{N,S}(sum)
end

function bsplineinvariant(numknots::Int, numcontrolpoints::Int, curveorder::Int)
    m = numknots - 1
    n = numcontrolpoints - 1
    @assert m == n + curveorder + 1 "This invariant should hold: m = n + curveorder + 1. Actual values: m=$m n=$n curveorder=$curveorder"
end

# Note TO SELF: need code to extract curves in the u and v directions from BSpline surfaces to work with the conversion to Bezier form.
# "converts a NURBS curve into Bezier curve segments, represented as Bezier control points. For an m degree NURBS with n knots this will return n - (m+1) + 1 Bezier segments, each with m+1 control points."
# function toBezier(curve::BSplineCurve{T,S,N,M}) where {T,S,N,M}
#     #this code follows that in the NURBS book closely. Use OffsetArrays to index by zero to avoid the nightmare of converting to indexing by 1.
#     U = curve.knots
#     Pw = curve.controlpolygon

#     Qw = OffsetArray{Array{OffestArray{Array{SVector{N,S},1}},1}}(undef,0)
#     temppoints = Array{SVector{N,S},1}(undef,0)

#     m = M + numknots(curve) + 1
#     a = M
#     b = M+1
#     nb = 0

#     for i in 0:M
#         push!(temppoints,Pw[i])
#     end
#     Qw[nb] = temppoints

#     while b < m
#         i = b
#         while b < m && U[b+1] == U[b]
#             b += 1;
#         end
#         mult = b-i+1
#         if mult < M
#             numerator = U[b] - U[a] #numerator of alpha
#             #compute and store alphas

#             for j in p:-1:mult+1
#                 alphas[j-mult-1] = numerator/(U[a+j] - U[a])
#             end

#             r = p - mult #insert knot r times
#             for j in 1:r
#                 save = r-j
#                 s = mult+j # This many new points
#                 for k in p:-1:s
#                     alpha = alphas[k-s]
#                     Qw[nb][k] = alpha*Qw[nb][k] + (1-alpha)*Qw[nb][k-1]
#                 end

#                 if b < m #control point of next segment
#                     Qw[nb+1][save] = Qw[nb][p]
#                 end
#             end
#         end

#         nb += 1 #Bezier segment completed
#         if b < m  #initialize for next segment
#             for i in p-mult:pairs
#                 Qw[nb][i] = Pw[b-p+1]
#             end
#             a = b
#             b = b+1
#         end
#     end
# end
