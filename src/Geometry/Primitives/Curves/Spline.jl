# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

"""
Either `Rational` or `Euclidean`, used for [`Spline`](@ref)s and [`SplineSurface`](@ref)s.
"""
abstract type CurveType end
abstract type Rational <: CurveType end
abstract type Euclidean <: CurveType end
export CurveType

"""
    Spline{P<:CurveType,S<:Number,N,M}

`M` is the curve order, i.e., the highest power of the parameterizing variable, u.
`P` determines the [`CurveType`](@ref).

All Spline types must implement:
```julia
point(curve,u)
```
and have field `controlpolygon`
"""
abstract type Spline{P<:CurveType,S<:Number,N,M} end

"""
    SplineSurface{P,S,N,M} <: ParametricSurface{S,N}

Curve order, `M`, is the same in the u and v direction and fixed over all spans.
`P` determines the [`CurveType`](@ref).
"""
abstract type SplineSurface{P<:CurveType,S,N,M} <: ParametricSurface{S,N} end
export Spline, SplineSurface

function toeuclidean(point)
    return SVector{length(point) - 1}([point[i] / point[end] for i in 1:(lastindex(point) - 1)])
end

function euclideancontrolpoints(curve::Spline{Rational,S,N,M}) where {S,N,M}
    return toeuclidean.(curve.controlpolygon)
end

function euclideancontrolpoints(curve::Spline{Euclidean,S,N,M}) where {S,N,M}
    return curve.controlpolygon
end

#Bernstein polynomial functions
#B(i,n,u) is broken into Bc,Bu because these two functions will be needed when constructing the matrics to compute moving lines or moving planes
Bc(i, n) = factorial(n) ÷ (factorial(i) * factorial(n - i)) #coefficient of the i,n Bernstein basis polynomial
Bu(i, n, u) = u^i * (1 - u)^(n - i) #polynomial part of Bernstein basis polynomial
B(i, n, u) = Bc(i, n) * Bu(i, n, u)
