# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

"""
    BezierCurve{P,S,N,M} <: Spline{P,S,N,M}

`N` is the dimension of the curve, `M` is the curve order

```julia
BezierCurve{P,S,N,M}(controlpoints::AbstractArray{<:AbstractArray{S,1},1})
```
"""
struct BezierCurve{P,S,N,M} <: Spline{P,S,N,M} # this has a funny name because of name conflict with Plots.
    controlpolygon::Array{SVector{N,S},1}

    function BezierCurve{P,S,N,M}(controlpoints::AbstractArray{<:AbstractArray{S,1},1}) where {P<:CurveType,S<:Number,N,M}
        @assert length(controlpoints) == M + 1 && length(controlpoints[1]) == N
        return new{P,S,N,M}([SVector{N,S}(point) for point in controlpoints])
    end
end
export BezierCurve

# m is the order of the moving lines, n is the order of the curve, d is the spatial dimension of the curve
_arraysize(m, n, d) = (2m + n + 2, (m + 1) * (d + 1))
# computes the proper row index so control point Pj gets multiplied by correct non-redundant basis function.
# The m term is to offset all the Pj values below the first m+1 rows to allow for the correct summation of the final g(d+1,k) terms.
_rowindex(k, j, movinglineorder) = k + j + movinglineorder + 2

function movinglinearray(curve::BezierCurve{P,S,N,M}) where {P,S,N,M}
    # curveorder = M, spatialdimension = N
    movinglineorder = M - 1 # minimum order of moving line polynomial necessary to represent curve of order M.
    result = zeros(S, _arraysize(movinglineorder, M, N))
    ctrlpts = curve.controlpolygon

    # upper right corner will contain just Bc(k,m) coefficients that get multiplied by the moving line vector g[movinglineorder*(d+1) + k] term
    start = movinglineorder * (N + 1) + 1

    for index in start:(start + movinglineorder)
        k = index - start

        println("$k $movinglineorder $index")

        result[k + 1, index] = Bc(k, movinglineorder)
    end

    for spatialindex in 0:(N - 1)
        colblock = spatialindex * (movinglineorder + 1) + 1
        for j in 0:M
            for k in 0:movinglineorder
                row = _rowindex(k, j, movinglineorder)
                result[row, colblock + k] = Bc(j, M) * Bc(k, movinglineorder) * ctrlpts[j + 1][spatialindex + 1]
            end
        end
    end
    return result
end

###################################################################################################################################

"""
    BezierSurface{P,S,N,M} <: SplineSurface{P,S,N,M}

Bezier surface defined by grid of control points.

!!! danger
    This surface does not create a valid half-space, requires updates to function correctly.

```julia
BezierSurface{P,S,N,M}(controlpoints::AbstractArray{<:AbstractArray{S,1},2})
```
"""
struct BezierSurface{P,S,N,M} <: SplineSurface{P,S,N,M}
    controlpolygon::Array{SVector{N,S},2}

    function BezierSurface{P,S,N,M}(controlpoints::AbstractArray{<:AbstractArray{S,1},2}) where {P<:CurveType,S<:Real,N,M}
        @assert (M + 1, M + 1) == size(controlpoints) && length(controlpoints[1, 1]) == N
        return new([SVector{N,S}(point) for point in controlpoints])
    end

    # function BezierSurface{P,S,N,M}(controlpoints::Array{SVector{N,S},2}) where {P<:CurveType,S<:Number,N,M}
    #     @assert (M + 1, M + 1) == size(controlpoints)
    #     newpoints = [copy(point) for point in controlpoints]
    #     # temp1 = SVector{M + 1,SVector{N,S}}(undef, M + 1)
    #     # temp2 = SVector{M + 1,SVector{N,S}}(undef, M + 1)
    #     return new(newpoints)
    # end
end
export BezierSurface

# """this function transforms the control points of the Bezier patch before ray tracing. Probably a better idea to use transforms in CSG operators,
# although this could be a little less computation at run time."""
# function Base.:*(a::Transform{T}, surf::BezierSurface{P,T,N,M}) where {P,T,N,M}
#     result = deepcopy(surf)
#     for index in CartesianIndices(surf.controlpolygon)
#         result.controlpolygon[index] = a * result.controlpolygon[index]
#     end
#     return result
# end

deepcopy(a::BezierSurface{P,S,N,M}) where {P,S,N,M} = BezierSurface{P,S,N,M}(a.controlpolygon)
uvrange(::Type{BezierSurface{P,S}}) where {P,S<:Real} = ((zero(S), one(S)), (zero(S), one(S)))
uvrange(::BezierSurface{P,S}) where {P,S<:Real} = uvrange(BezierSurface{P,S})
onsurface(::BezierSurface{P,S,3,M}, ::SVector{3,S}) where {P,S<:Real,M} = false

# DECASTELJAU
@inline function decasteljau(controlpoints::AbstractArray{SVector{N,S},1}, ::Val{2}, u::T) where {T<:Real,S,N}
    w = 1 - u
    return (controlpoints[1] * w + controlpoints[2] * u) * w + (controlpoints[2] * w + controlpoints[3] * u) * u
end

@inline function decasteljau(controlpoints::AbstractArray{SVector{N,S},1}, ::Val{3}, u::T) where {T<:Real,S,N}
    w = 1 - u
    np1 = controlpoints[1] * w + controlpoints[2] * u
    np2 = controlpoints[2] * w + controlpoints[3] * u
    np3 = controlpoints[3] * w + controlpoints[4] * u
    np1 = np1 * w + np2 * u
    np2 = np2 * w + np3 * u
    np1 = np1 * w + np2 * u
    return np1
end

@inline function decasteljau(controlpoints::AbstractArray{SVector{N,S},1}, ::Val{4}, u::T) where {T<:Real,S,N}
    w = 1 - u
    np1 = controlpoints[1] * w + controlpoints[2] * u
    np2 = controlpoints[2] * w + controlpoints[3] * u
    np3 = controlpoints[3] * w + controlpoints[4] * u
    np4 = controlpoints[4] * w + controlpoints[5] * u
    np1 = np1 * w + np2 * u
    np2 = np2 * w + np3 * u
    np3 = np3 * w + np4 * u
    np1 = np1 * w + np2 * u
    np2 = np2 * w + np3 * u
    np1 = np1 * w + np2 * u
    return np1
end

@inline function decasteljau(controlpoints::AbstractArray{SVector{N,S},1}, ::Val{M}, u::T) where {T<:Real,S,N,M}
    # general solution for decasteljau for arbitrary curve order
    newpoints = MVector{M + 1,SVector{N,T}}(controlpoints)
    # do this to force newpoints to have the type of u rather than controlpoints. When u is a dual number this makes the array the correct type.
    @inbounds for i in 1:M
        @inbounds for k in 1:(M + 1 - i)
            newpoints[k] = newpoints[k] * (1 - u) + newpoints[k + 1] * u
            # this kills forwarddiff and zygote. ForwardDiff doesn't work because newpoints has the type of controlpoints, which is going to be Float64 most times,
            # when it needs to be of type dual number. Zygote doesn't support mutable arrays. Arrrrgggghhh! Useless!
        end
    end
    return newpoints[1]
end

# POINT
function point(surf::BezierSurface{P,S,N,3}, u::T, v::T) where {T<:Real,P,S,N}
    q1 = decasteljau(view(surf.controlpolygon, :, 1), Val(3), u)
    q2 = decasteljau(view(surf.controlpolygon, :, 2), Val(3), u)
    q3 = decasteljau(view(surf.controlpolygon, :, 3), Val(3), u)
    q4 = decasteljau(view(surf.controlpolygon, :, 4), Val(3), u)
    return decasteljau(SVector{4,SVector{N,S}}(q1, q2, q3, q4), Val(3), v)
end

function point(surf::BezierSurface{P,S,N,4}, u::T, v::T) where {T<:Real,P,S,N}
    q1 = decasteljau(view(surf.controlpolygon, :, 1), Val(4), u)
    q2 = decasteljau(view(surf.controlpolygon, :, 2), Val(4), u)
    q3 = decasteljau(view(surf.controlpolygon, :, 3), Val(4), u)
    q4 = decasteljau(view(surf.controlpolygon, :, 4), Val(4), u)
    q5 = decasteljau(view(surf.controlpolygon, :, 5), Val(4), u)
    return decasteljau(SVector{5,SVector{N,S}}(q1, q2, q3, q4, q5), Val(4), v)
end

function point(surf::BezierSurface{P,S,N,M}, u::T, v::T) where {T<:Real,P,S,N,M}
    qi = zeros(MVector{M + 1,SVector{N,S}})
    @inbounds for i in 1:(M + 1)
        qi[i] = decasteljau(view(surf.controlpolygon, :, i), Val(M), u)
    end
    return decasteljau(qi, Val(M), v)
end

# DERIVATIVE
@inline function derivative(controlpoints::SVector{4,SVector{N,S}}, u::T) where {T<:Real,N,S}
    d1 = controlpoints[2] .- controlpoints[1]
    d2 = controlpoints[3] .- controlpoints[2]
    d3 = controlpoints[4] .- controlpoints[3]
    return 3 * decasteljau(SVector{3,SVector{N,S}}(d1, d2, d3), Val(2), u)
end

@inline function derivative(controlpoints::SVector{5,SVector{N,S}}, u::T) where {T<:Real,N,S}
    d1 = controlpoints[2] .- controlpoints[1]
    d2 = controlpoints[3] .- controlpoints[2]
    d3 = controlpoints[4] .- controlpoints[3]
    d4 = controlpoints[5] .- controlpoints[4]
    return 4 * decasteljau(SVector{4,SVector{N,S}}(d1, d2, d3, d4), Val(3), u)
end

@inline function derivative(controlpoints::MVector{Q,SVector{N,S}}, u::T) where {T<:Real,N,S,Q}
    derivpoints = zeros(MVector{Q - 1,SVector{N,T}})
    @inbounds for i in 1:(Q - 1)
        derivpoints[i] = controlpoints[i + 1] .- controlpoints[i]
    end
    return (Q - 1) * decasteljau(derivpoints, Val(Q - 2), u)
end

# PARTIALS
function partials(surf::BezierSurface{Euclidean,S,N,3}, u::T, v::T) where {T<:Real,S,N}
    q1 = decasteljau(view(surf.controlpolygon, :, 1), Val(3), u)
    q2 = decasteljau(view(surf.controlpolygon, :, 2), Val(3), u)
    q3 = decasteljau(view(surf.controlpolygon, :, 3), Val(3), u)
    q4 = decasteljau(view(surf.controlpolygon, :, 4), Val(3), u)
    partialv = derivative(SVector{4,SVector{N,S}}(q1, q2, q3, q4), v)
    q1 = decasteljau(view(surf.controlpolygon, 1, :), Val(3), v)
    q2 = decasteljau(view(surf.controlpolygon, 2, :), Val(3), v)
    q3 = decasteljau(view(surf.controlpolygon, 3, :), Val(3), v)
    q4 = decasteljau(view(surf.controlpolygon, 4, :), Val(3), v)
    partialu = derivative(SVector{4,SVector{N,S}}(q1, q2, q3, q4), u)
    return (partialu, partialv)
end

function partials(surf::BezierSurface{Euclidean,S,N,4}, u::T, v::T) where {T<:Real,S,N}
    q1 = decasteljau(view(surf.controlpolygon, :, 1), Val(4), u)
    q2 = decasteljau(view(surf.controlpolygon, :, 2), Val(4), u)
    q3 = decasteljau(view(surf.controlpolygon, :, 3), Val(4), u)
    q4 = decasteljau(view(surf.controlpolygon, :, 4), Val(4), u)
    q5 = decasteljau(view(surf.controlpolygon, :, 5), Val(4), u)
    partialv = derivative(SVector{5,SVector{N,S}}(q1, q2, q3, q4, q5), v)
    q1 = decasteljau(view(surf.controlpolygon, 1, :), Val(4), v)
    q2 = decasteljau(view(surf.controlpolygon, 2, :), Val(4), v)
    q3 = decasteljau(view(surf.controlpolygon, 3, :), Val(4), v)
    q4 = decasteljau(view(surf.controlpolygon, 4, :), Val(4), v)
    q5 = decasteljau(view(surf.controlpolygon, 5, :), Val(4), v)
    partialu = derivative(SVector{5,SVector{N,S}}(q1, q2, q3, q4, q5), u)
    return (partialu, partialv)
end

function partials(surf::BezierSurface{Euclidean,S,N,M}, u::T, v::T) where {T<:Real,S,N,M}
    qi = zeros(MVector{M + 1,SVector{N,T}})
    # do this to force newpoints to have the type of u rather than controlpoints. When u is a dual number this makes the array the correct type.
    @inbounds for i in 1:(M + 1)
        qi[i] = decasteljau(view(surf.controlpolygon, :, i), Val(M), u)
    end
    partialv = derivative(qi, v)
    @inbounds for i in 1:(M + 1)
        qi[i] = decasteljau(view(surf.controlpolygon, i, :), Val(M), v)
    end
    partialu = derivative(qi, u)
    return (partialu, partialv)
end

# "computes partials and point simultaneously for Euclidean surface"
# function pointandpartials(surf::BezierSurface{Euclidean,S,N,M}, u::T, v::T) where {T<:Real,S,N,M}
#     qi = zeros(MVector{M + 1,SVector{N,T}})
#     # do this to force newpoints to have the type of u rather than controlpoints. When u is a dual number this makes the array the correct type.
#     @inbounds for i in 1:(M + 1)
#         qi[i] = decasteljau(view(surf.controlpolygon, :, i), M, u)
#     end
#     pt = decasteljau(qi, M, v)
#     partialv = derivative(qi, v)
#     @inbounds for i in 1:(M + 1)
#         qi[i] = decasteljau(view(surf.controlpolygon, i, :), M, v)
#     end
#     partialu = derivative(qi, u)
#     return (pt, partialu, partialv)
# end

function partials(surf::BezierSurface{Rational,S,N,M}, u::T, v::T) where {T<:Real,S,N,M}
    # computes Euclidean space partials for Rational surface. Not consistent with point, which has euclideanpoint function, but not sure why user would ever want Rational partials.
    pt = point(surf, u, v)
    du, dv = partials(surf, u, v)
    rdu = rationalcorrection(pt, du)
    rdv = rationalcorrection(pt, dv)
    return (rdu, rdv)
end

function normal(surf::BezierSurface, u::T, v::T) where {T<:Real}
    du, dv = partials(surf, u, v)
    return normalize(cross(du, dv))
end

function rationalcorrection(pt::SVector{N,S}, dpt::SVector{N,S}) where {N,S}
    # have form A/B where A = (x,y,z) and B = w component of vector. Want partial(A/B,u). Get (B*partial(A,u) - A*partial(B,u))/B*B
    A = pt[1:(N - 1)]
    B = pt[N]
    dB = dpt[N]
    dA = dpt[1:(N - 1)]
    return (B .* dA - A .* dB) / B^2
end

# """attempt to generate inline code so Zygote and Forwarddiff could differentatiate properly. Didn't work."""
# function inlinedecasteljau(curveorder, N)
#     a = "(u,$(reduce(*, ("newpoints$(k)_$j," for k in 1:curveorder,j in 1:N)))) -> let
#     "
#     for i in 1:curveorder
#         for k in 1:(curveorder + 1 - i)
#             for j in 1:N
#                 a *= "newpoints$(k)_$j = newpoints$(k)_$j * (1-u) + newpoints$(k + 1)_$j * u
#                  "
#             end
#         end
#     end
#     a *= "return ($(reduce(*, ("newpoints1_$j," for j in 1:N))))
#     end"

#     return eval(Meta.parse(a))
# end

# "Creates a differentiable function to evaluate a surface point. Zygote has trouble with mutable arrays so this creates a function
# which unrolls all loops and converts array indexing to variables. Stopped working on this because it became so diffult to make it work."
# function pointinline(surf::BezierSurface{P,S,N,M}, u::T, v::T) where {T<:Real,P,S,N,M}
#     controlpoints = surf.controlpolygon
#     qi = zeros(MVector{N * (M + 1),MVector{N,S},1})
#     pointinl = inlinedecasteljau(M, N)
#     @inbounds for i in 1:(M + 1)
#         start = (i - 1) * N + 1
#         @inbounds for j in 1:N
#             temp = pointinl(u, controlpoints[:, i]...)
#             qi[start + j - 1] = temp[j]
#         end
#     end
#     return pointinl(v, qi...)
# end

function point(curve::BezierCurve{P,S,N,M}, u::T) where {T<:Real,P,S,N,M}
    # Will return homogeneous curve point if curve is homogeneous and euclidean point if curve is Euclidean.
    return decasteljau(curve.controlpolygon, M, u)
end

euclideanpoint(curve::BezierCurve{Euclidean,S,N,M}, u::T) where {T<:Real,S,N,M} = point(curve, u)

function euclideanpoint(curve::BezierCurve{Rational,S,N,M}, u::T) where {T<:Real,S,N,M}
    temp = point(curve, u)
    return toeuclidean(temp)
end

#########################################################################################################################

function inside(surf::AcceleratedParametricSurface{S,3,BezierSurface{P,S,3,M}}, p::SVector{3,S}) where {P,S<:Real,M}
    # approximate the normal to the bezier surface as a whole by taking the average of the normals of the triangles
    # defined by the extreme points of the surface
    A = point(surf, 0.0, 0.0)
    B = point(surf, 0.0, 1.0)
    C = point(surf, 1.0, 0.0)
    D = point(surf, 1.0, 1.0)
    n1 = cross(C - A, B - A)
    n2 = cross(B - D, C - D)
    approx_normal = (n1 + n2) / (norm(n1) + norm(n2))
    r = Ray(p, approx_normal)
    @inbounds for i in 1:length(surf.triangles)
        intv = surfaceintersection(surf.triangles[i], r)
        if !(intv isa EmptyInterval)
            return true
        end
    end
    return false
end

function BoundingBox(surf::BezierSurface{P,S,N,M}) where {P,S,N,M}
    big = fill(typemin(S), MVector{N,S})
    small = fill(typemax(S), MVector{N,S})

    for pt in surf.controlpolygon
        small .= min.(small, pt)
        big .= max.(big, pt)
    end

    return BoundingBox(small[1], big[1], small[2], big[2], small[3], big[3])
end
