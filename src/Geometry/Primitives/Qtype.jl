"""
Module to enclose QType polynomial specific functionality. For reference see:
- [_Robust, efficient computational methods for axially symmetric optical aspheres_ - G. W. Forbes, 2010](https://www.osapublishing.org/viewmedia.cfm?uri=oe-18-19-19700&seq=0)
- [_Characterizing the shape of freeform optics_ - G. W. Forbes, 2012](https://www.osapublishing.org/viewmedia.cfm?uri=oe-20-3-2483&seq=0)

"""
module QType
using StaticArrays
using Plots
using ..Opticks: QTYPE_PRECOMP

function F(m::Int, n::Int)::Float64
    @assert m > 0
    if n === 0
        return Float64((m^2 * factorial2(big(2m - 3))) / (2^(m + 1) * factorial(big(m - 1))))
    elseif n > 0 && m === 1
        return Float64((4(n - 1)^2 * n^2 + 1) / (8 * (2n - 1)^2) + 11 * Int(n === 1) / 32)
    elseif n > 0 && m > 1
        χ = m + n - 2
        return Float64((2n * χ * (3 - 5m + 4n * χ) + m^2 * (3 - m + 4n * χ)) / ((m + 2n - 3) * (m + 2n - 2) * (m + 2n - 1) * (2n - 1)) * γ(m, n))
    else
        throw(ErrorException("Invalid n and m"))
    end
end

γ(a::Int, b::Int) = (factorial(big(b)) * factorial2(big(2a + 2b - 3))) / (2^(a + 1) * factorial(big(a + b - 3)) * factorial2(big(2b - 1)))

function G(m::Int, n::Int)::Float64
    @assert m > 0
    if n === 0
        return Float64(factorial2(big(2m - 1)) / (2^(m + 1) * factorial(big(m - 1))))
    elseif n > 0 && m === 1
        return Float64(-((2 * n^2 - 1) * (n^2 - 1)) / (8 * (4 * n^2 - 1)) - Int(n === 1) / 24)
    elseif n > 0 && m > 0
        return Float64(-((2n * (m + n - 1) - m) * (n + 1) * (2m + 2n - 1)) / ((m + 2n - 2) * (m + 2n - 1) * (m + 2n) * (2n + 1)) * γ(m, n))
    else
        throw(ErrorException("Invalid n and m"))
    end
end

function factorial2(n::I)::I where {I<:Signed}
    @assert n <= 21 || n isa BigInt # not sure what number is the limit but this should do it...
    if n <= 0
        return I(1)
    else
        return n * factorial2(n - 2)
    end
end

function f(m::Int, n::Int, force::Bool = false)::Float64
    @assert m > 0
    if m <= QTYPE_PRECOMP && n < QTYPE_PRECOMP && !force
        return @inbounds PRECOMP_f[m, n + 1]
    else
        if n === 0
            return sqrt(F(m, 0))
        else
            return sqrt(F(m, n) - g(m, n - 1, force)^2)
        end
    end
end

function g(m::Int, n::Int, force::Bool = false)::Float64
    @assert m > 0
    if m <= QTYPE_PRECOMP && n < QTYPE_PRECOMP && !force
        return @inbounds PRECOMP_g[m, n + 1]
    else
        return G(m, n) / f(m, n, force)
    end
end

@inline function A(m::Int, n::Int, force::Bool = false)::Float64
    @assert m > 0
    if m <= QTYPE_PRECOMP && n < QTYPE_PRECOMP && !force
        return @inbounds PRECOMP_A[m, n + 1]
    else
        if m === 1 && n === 0
            return 2.0
        elseif m === 1 && n === 1
            return -4.0 / 3.0
        elseif n === 0 && m > 1
            return 2m - 1.0
        else
            return (2n - 1) * (m + 2n - 2) * (4n * (m + n - 2) + (m - 3) * (2m - 1)) / D(m, n)
        end
    end
end

@inline function B(m::Int, n::Int, force::Bool = false)::Float64
    @assert m > 0
    if m <= QTYPE_PRECOMP && n < QTYPE_PRECOMP && !force
        return @inbounds PRECOMP_B[m, n + 1]
    else
        if n === 1 && m === 1
            return -8.0 / 3.0
        elseif n === 0
            if m === 1
                return -1.0
            else
                return 2.0 * (1.0 - m)
            end
        else
            return -2.0 * (2n - 1) * (m + 2n - 3) * (m + 2n - 2) * (m + 2n - 1) / D(m, n)
        end
    end
end

@inline function C(m::Int, n::Int, force::Bool = false)::Float64
    @assert m > 0
    if n === 0
        return NaN
    elseif m <= QTYPE_PRECOMP && n < QTYPE_PRECOMP && !force
        return @inbounds PRECOMP_C[m, n + 1]
    else
        if m === 1 && n === 1
            return -11.0 / 3.0
        elseif m === 1 && n === 2
            return 0.0
        else
            return n * (2n - 3) * (m + 2n - 1) * (2m + 2n - 3) / D(m, n)
        end
    end
end

@inline function D(m::Int, n::Int)::Float64
    @assert m > 0
    return Float64((4n^2 - 1) * (m + n - 2) * (m + 2n - 3))
end

const PRECOMP_g = Matrix{Float64}(reshape(collect(g(m, n, true) for n in 0:(QTYPE_PRECOMP - 1) for m in 1:QTYPE_PRECOMP), (QTYPE_PRECOMP, QTYPE_PRECOMP)))
const PRECOMP_f = Matrix{Float64}(reshape(collect(f(m, n, true) for n in 0:(QTYPE_PRECOMP - 1) for m in 1:QTYPE_PRECOMP), (QTYPE_PRECOMP, QTYPE_PRECOMP)))
const PRECOMP_A = Matrix{Float64}(reshape(collect(A(m, n, true) for n in 0:(QTYPE_PRECOMP - 1) for m in 1:QTYPE_PRECOMP), (QTYPE_PRECOMP, QTYPE_PRECOMP)))
const PRECOMP_B = Matrix{Float64}(reshape(collect(B(m, n, true) for n in 0:(QTYPE_PRECOMP - 1) for m in 1:QTYPE_PRECOMP), (QTYPE_PRECOMP, QTYPE_PRECOMP)))
const PRECOMP_C = Matrix{Float64}(reshape(collect(C(m, n, true) for n in 0:(QTYPE_PRECOMP - 1) for m in 1:QTYPE_PRECOMP), (QTYPE_PRECOMP, QTYPE_PRECOMP)))

"""
    S(coeffs::SVector{NP1,T}, m::Int x::T) -> T

Evaluates ``\\sum_{n=0}^{N}c_n^mQ_n^m(x)`` where ``c_n^m`` is either an ``\\alpha`` or ``\\beta`` QType coefficient and ``m \\gt 0``.
"""
function S(coeffs::SVector{NP1,T}, m::Int, x::T)::T where {T<:Real,NP1}
    @assert m > 0
    if all(iszero.(coeffs))
        return zero(T)
    end
    N = NP1 - 1 # offset for indexing
    if N === 0
        return coeffs[1]
    end
    αₙ₊₂ = zero(T)
    αₙ₊₁ = zero(T)
    αₙ = zero(T)
    @inbounds for n in N:-1:0
        αₙ = coeffs[n + 1] + (A(m, n) + B(m, n) * x) * αₙ₊₁ - C(m, n + 1) * αₙ₊₂
        if n > 0
            αₙ₊₂ = αₙ₊₁
            αₙ₊₁ = αₙ
        end
    end
    if m === 1 && N > 2
        return αₙ / 2 - 2 / 5 * αₙ₊₂
    else
        return αₙ / 2
    end
end

"""
    dSdx(coeffs::SVector{NP1,T}, x::T) -> T

Evaluates ``\\frac{\\partial}{\\partial x}\\sum_{n=0}^{N}c_n^mQ_n^m(x)`` where ``c_n^m`` is either an ``\\alpha`` or ``\\beta`` QType coefficient and ``m \\gt 0``.
"""
function dSdx(coeffs::SVector{NP1,T}, m::Int, x::T)::T where {T<:Real,NP1}
    @assert m > 0
    if all(iszero.(coeffs))
        return zero(T)
    end
    N = NP1 - 1 # offset for indexing
    if N === 0
        return zero(T)
    end
    αₙ₊₂ = zero(T)
    αₙ₊₁ = zero(T)
    dαₙ₊₂ = zero(T)
    dαₙ₊₁ = zero(T)
    dαₙ = zero(T)
    @inbounds for n in N:-1:0
        # calcualte deriv
        Bmn = B(m, n)
        k = (A(m, n) + Bmn * x)
        Cmnp1 = C(m, n + 1)
        if n < N # derivative calcualted from N-1 to 0, so ignore first iteration
            dαₙ = Bmn * αₙ₊₁ + k * dαₙ₊₁ - Cmnp1 * dαₙ₊₂
            if n > 0
                dαₙ₊₂ = dαₙ₊₁
                dαₙ₊₁ = dαₙ
            end
        end
        # calculate height for use in next iter deriv
        αₙ = coeffs[n + 1] + k * αₙ₊₁ - Cmnp1 * αₙ₊₂
        αₙ₊₂ = αₙ₊₁
        αₙ₊₁ = αₙ
    end
    if m === 1 && N > 2
        return dαₙ / 2 - 2 / 5 * dαₙ₊₂
    else
        return dαₙ / 2
    end
end

# for cases where m == 0 below:

function f0(n::Int, force::Bool = false)::Float64
    if n < QTYPE_PRECOMP && !force
        return @inbounds PRECOMP_f0[n + 1]
    elseif n === 0
        return 2.0
    elseif n === 1
        return sqrt(19) / 2
    else
        return sqrt(n * (n + 1) + 3.0 - g0(n - 1, force)^2 - h0(n - 2, force)^2)
    end
end

function g0(n::Int, force::Bool = false)::Float64
    if n < QTYPE_PRECOMP && !force
        return @inbounds PRECOMP_g0[n + 1]
    elseif n === 0
        return -0.5
    else
        return -(1 + g0(n - 1, force) * h0(n - 1, force)) / f0(n, force)
    end
end

function h0(n::Int, force::Bool = false)::Float64
    if n < QTYPE_PRECOMP && !force
        return @inbounds PRECOMP_h0[n + 1]
    else
        return -(n + 2) * (n + 1) / (2 * f0(n, force))
    end
end

const PRECOMP_g0 = Vector{Float64}(collect(g0(n, true) for n in 0:(QTYPE_PRECOMP - 1)))
const PRECOMP_f0 = Vector{Float64}(collect(f0(n, true) for n in 0:(QTYPE_PRECOMP - 1)))
const PRECOMP_h0 = Vector{Float64}(collect(h0(n, true) for n in 0:(QTYPE_PRECOMP - 1)))

"""
    S0(coeffs::SVector{NP1,T}, x::T) -> T

Evaluates ``\\sum_{n=0}^{N}\\alpha_n^0Q_n^0(x)``.
"""
function S0(coeffs::SVector{NP1,T}, x::T)::T where {T<:Real,NP1}
    if all(iszero.(coeffs))
        return zero(T)
    end
    N = NP1 - 1 # offset for indexing
    if N === 0
        return coeffs[1]
    end
    k = T(2) - 4 * x
    αₙ₊₂ = coeffs[N + 1]
    αₙ₊₁ = coeffs[N] + k * αₙ₊₂
    for n in (N - 2):-1:0
        αₙ = coeffs[n + 1] + k * αₙ₊₁ - αₙ₊₂
        αₙ₊₂ = αₙ₊₁
        αₙ₊₁ = αₙ
    end
    return 2 * (αₙ₊₁ + αₙ₊₂) # a0 and a1
end

"""
    dS0dx(coeffs::SVector{NP1,T}, x::T) -> T

Evaluates ``\\frac{\\partial}{\\partial x}\\sum_{n=0}^{N}\\alpha_n^0Q_n^0(x)``.
"""
function dS0dx(coeffs::SVector{NP1,T}, x::T)::T where {T<:Real,NP1}
    if all(iszero.(coeffs))
        return zero(T)
    end
    N = NP1 - 1 # offset for indexing
    if N == 0
        return zero(T)
    end
    k = T(2) - 4 * x
    αₙ₊₂ = coeffs[N + 1]
    αₙ₊₁ = coeffs[N] + k * αₙ₊₂
    dαₙ₊₂ = zero(T)
    dαₙ₊₁ = -4 * αₙ₊₂
    for n in (N - 2):-1:0
        dαₙ = k * dαₙ₊₁ - dαₙ₊₂ - 4 * αₙ₊₁
        dαₙ₊₂ = dαₙ₊₁
        dαₙ₊₁ = dαₙ
        αₙ = coeffs[n + 1] + k * αₙ₊₁ - αₙ₊₂
        αₙ₊₂ = αₙ₊₁
        αₙ₊₁ = αₙ
    end
    return 2 * (dαₙ₊₁ + dαₙ₊₂) # a0 and a1
end

########################################
# for testing, not used in actual calculation because there is a more efficient method to find the sum directly (as in S and S0)

function P(m::Int, n::Int, x::T)::T where {T<:Real}
    if m === 0
        if n === 0
            return T(2)
        elseif n === 1
            return T(6) - 8 * x
        else
            return (T(2) - 4 * x) * P(0, n - 1, x) - P(0, n - 2, x)
        end
    else
        if n === 0
            return one(T) / 2
        elseif n === 1
            if m === 1
                return one(T) - x / 2
            else
                return m - 1 / 2 - (m - 1) * x
            end
        else
            return (A(m, n - 1) + B(m, n - 1) * x) * P(m, n - 1, x) - C(m, n - 1) * P(m, n - 2, x)
        end
    end
end

function Q(m::Int, n::Int, x::T)::T where {T<:Real}
    if m === 0
        if n === 0
            return one(T)
        elseif n === 1
            return (T(13) - 16x) / sqrt(19)
        else
            return (P(0, n, x) - g0(n - 1) * Q(0, n - 1, x) - h0(n - 2) * Q(0, n - 2, x)) / f0(n)
        end
    else
        if n === 0
            return 1 / (2 * f(m, 0))
        else
            return (P(m, n, x) - g(m, n - 1) * Q(m, n - 1, x)) / f(m, n)
        end
    end
end

function plot!(m::Int, n::Int)
    @assert m >= 0 && n >= 0
    us = 0:0.01:1
    vs = []
    for u in us
        if m === 0
            push!(vs, u^2 * (1 - u^2) * Q(m, n, u^2))
        else
            push!(vs, u^m * Q(m, n, u^2))
        end
    end
    Plots.plot!(us, vs, label = m === 0 ? "u^2 * (1-u^2) * Q($m, $n, u^2)" : "u^m * Q($m, $n, u^2)")
end

end # module QType

#########################################################################################################

"""
    QTypeSurface{T,D,M,N} <: ParametricSurface{T,D}

Surface incorporating the QType polynomials - radius and conic are defined relative to absolute semi-diameter, QType terms are normalized according to the `normradius` parameter.
`T` is the datatype, `D` is the dimensionality, `M` and `N` are the maximum QType terms used.

The surface is centered at the origin and treated as being the cap of an infinite cylinder, thus creating a true half-space.
Outside of 0 <= ρ <= 1 the height of the surface is not necessarily well defined, so NaN may be returned.

```julia
QTypeSurface(semidiameter; radius = Inf, conic = 0.0, αcoeffs = nothing, βcoeffs = nothing, normradius = semidiameter)
```

`αcoeffs` and `βcoeffs` should be a vector of tuples of the form `(m, n, v)` where `v` is the value of the coefficient ``α_n^m`` or ``β_n^m`` respectively.

The sag is defined by the equation

```math
\\begin{aligned}
z(r,\\phi) = & \\frac{cr^2}{1 + \\sqrt{1 - (1+k)c^2r^2}} + \\frac{\\sqrt{1 + kc^2r^2}}{\\sqrt{1-(1+k)c^2r^2}} \\cdot \\\\
             & \\left\\{ \\rho^2(1-\\rho^2)\\sum_{n=0}^{N}\\alpha_n^0 Q_n^0 (\\rho^2) + \\sum_{m=1}^{M}\\rho^m\\sum_{n=0}^N \\left[ \\alpha_n^m\\cos{m\\phi} +\\beta_n^m\\sin{m\\phi}\\right]Q_n^m(\\rho^2) \\right\\}
\\end{aligned}
```

where ``\\rho = \\frac{r}{\\texttt{normradius}}``, ``c = \\frac{1}{\\texttt{radius}}``, ``k = \\texttt{conic}`` and ``Q_n^m`` is the QType polynomial index ``m``, ``n``.
"""
struct QTypeSurface{T,D,M,N} <: ParametricSurface{T,D}
    semidiameter::T
    curvature::T
    conic::T
    boundingcylinder::Cylinder{T,D}
    b0coeffs::SVector{N,T}
    dαcoeffs::SVector{M,SVector{N,T}} # m goes from 1:M while n goes from 0:N-1
    dβcoeffs::SVector{M,SVector{N,T}}
    normradius::T
    maxheight::T

    function QTypeSurface(semidiameter::T; radius::T = typemax(T), conic::T = zero(T), αcoeffs::Union{Nothing,Vector{Tuple{Int,Int,T}}} = nothing, βcoeffs::Union{Nothing,Vector{Tuple{Int,Int,T}}} = nothing, normradius::T = semidiameter) where {T<:Real}
        @assert !isnan(semidiameter) && !isnan(radius) && !isnan(conic)
        @assert semidiameter > zero(T)
        @assert one(T) - (1 / radius)^2 * (conic + one(T)) * semidiameter^2 > 0 "Invalid surface (conic/radius combination: $radius, $conic)"
        # work out maximum coefficient value
        N = -1
        M = 0
        if αcoeffs !== nothing
            αcoeffs = αcoeffs::Vector{Tuple{Int,Int,T}}
            for α in αcoeffs
                m, n, v = α
                if n > N
                    N = n
                end
                if m > M
                    M = m
                end
            end
        end
        if βcoeffs !== nothing
            βcoeffs = βcoeffs::Vector{Tuple{Int,Int,T}}
            for β in βcoeffs
                m, n, v = β
                if n > N
                    N = n
                end
                if m > M
                    M = m
                end
            end
        end
        # process the inputs to get parameter matrices
        α0coeffs = zeros(MVector{N + 1,T})
        αcoeffsproc = zeros(MMatrix{M,N + 1,T})
        βcoeffsproc = zeros(MMatrix{M,N + 1,T})
        if αcoeffs !== nothing
            for α in αcoeffs
                m, n, v = α
                @assert m >= 0 && n >= 0
                if m == 0
                    α0coeffs[n + 1] = v
                else
                    αcoeffsproc[m, n + 1] = v
                end
            end
        end
        if βcoeffs !== nothing
            for β in βcoeffs
                m, n, v = β
                @assert m >= 0 && n >= 0
                βcoeffsproc[m, n + 1] = v
            end
        end
        # calculate the α term coefficients for m = 0
        b0coeffs = zeros(MVector{N + 1,T})
        if N >= 0
            bₙp2 = α0coeffs[N + 1] / QType.f0(N)
            b0coeffs[N + 1] = bₙp2
            if N > 0
                bₙp1 = (α0coeffs[N] - QType.g0(N - 1) * bₙp2) / QType.f0(N - 1)
                b0coeffs[N] = bₙp1
                if N > 1
                    for n in (N - 2):-1:0
                        bₙ = (α0coeffs[n + 1] - QType.g0(n) * bₙp1 - QType.h0(n) * bₙp2) / QType.f0(n)
                        b0coeffs[n + 1] = bₙ
                        bₙp2 = bₙp1
                        bₙp1 = bₙ
                    end
                end
            end
        end
        # m==0 max value is <0.4 for any n, otherwise max is 1 for any m and n, for simplicity just sum everything though it may not be the tightest
        maxheight = zero(T)
        maxheight += 0.4 * sum(abs.(α0coeffs))
        maxheight += 0.4 * sum(abs.(b0coeffs))
        @inbounds @simd for m in 1:M
            maxheight += sum(abs.(αcoeffsproc[m, :]))
            maxheight += sum(abs.(βcoeffsproc[m, :]))
        end
        # calculate the α and β term coefficients
        dαcoeffs = zeros(MVector{M,SVector{N + 1,T}})
        dβcoeffs = zeros(MVector{M,SVector{N + 1,T}})
        for m in 1:M
            thisα = zeros(MVector{N + 1,T})
            thisβ = zeros(MVector{N + 1,T})
            if N >= 0
                fmN = QType.f(m, N)
                lastdαₙ = αcoeffsproc[m, N + 1] / fmN
                thisα[N + 1] = lastdαₙ
                lastdβₙ = βcoeffsproc[m, N + 1] / fmN
                thisβ[N + 1] = lastdβₙ
                for n in (N - 1):-1:0
                    gmn = QType.g(m, n)
                    fmn = QType.f(m, n)
                    dαₙ = (αcoeffsproc[m, n + 1] - gmn * lastdαₙ) / fmn
                    dβₙ = (βcoeffsproc[m, n + 1] - gmn * lastdβₙ) / fmn
                    lastdαₙ = dαₙ
                    lastdβₙ = dβₙ
                    thisα[n + 1] = dαₙ
                    thisβ[n + 1] = dβₙ
                end
            end
            dαcoeffs[m] = SVector{N + 1,T}(thisα)
            dβcoeffs[m] = SVector{N + 1,T}(thisβ)
        end
        NP1 = N + 1
        new{T,3,M,NP1}(semidiameter, 1 / radius, conic, Cylinder(semidiameter, interface = opaqueinterface(T)), SVector{NP1,T}(b0coeffs), SVector{M,SVector{NP1,T}}(dαcoeffs), SVector{M,SVector{NP1,T}}(dβcoeffs), normradius, maxheight) # TODO!! incorrect interface on cylinder
    end
end
export QTypeSurface

uvrange(::Type{QTypeSurface{T,D,M,N}}) where {T<:Real,D,M,N} = ((zero(T), one(T)), (-T(π), T(π))) # ρ and θ

semidiameter(z::QTypeSurface{T}) where {T<:Real} = z.semidiameter

function point(z::QTypeSurface{T,3,M,N}, ρ::T, θ::T)::SVector{3,T} where {T<:Real,M,N}
    rad = z.semidiameter
    r = ρ * rad
    # ρ is normalised [0, 1]
    # r is absolute
    r2 = r^2
    ρ2 = ρ^2
    c2 = z.curvature^2
    q = (one(T) + z.conic) * c2 * r2
    if q > one(T)
        return SVector{3,T}(NaN, NaN, NaN)
    end
    sqrtq = sqrt(one(T) - q)
    h = z.curvature * r2 / (one(T) + sqrtq)
    if N > 0
        p = one(T) + z.conic * c2 * r2
        if p < zero(T)
            return SVector{3,T}(NaN, NaN, NaN)
        end
        u = ρ * rad / z.normradius
        u2 = u^2
        tot = u2 * (one(T) - u2) * QType.S0(z.b0coeffs, u2)
        @inbounds @simd for m in 1:M
            tot += u^m * (cos(m * θ) * QType.S(z.dαcoeffs[m], m, u2)::T + sin(m * θ) * QType.S(z.dβcoeffs[m], m, u2)::T)
        end
        h += (sqrt(p) / sqrtq) * tot
    end
    return SVector{3,T}(r * cos(θ), r * sin(θ), h)
end

function partials(z::QTypeSurface{T,3,M,N}, ρ::T, θ::T)::Tuple{SVector{3,T},SVector{3,T}} where {T<:Real,M,N}
    rad = z.semidiameter
    r = ρ * rad
    # ρ is normalised [0, 1]
    # r is absolute
    r2 = r^2
    ρ2 = ρ^2
    k = z.conic
    c = z.curvature
    c2 = c^2
    q = (one(T) + k) * c2 * r2
    if q > one(T)
        return SVector{3,T}(NaN, NaN, NaN), SVector{3,T}(NaN, NaN, NaN)
    end
    sqrtq = sqrt(one(T) - q)
    dhdρ = -((rad * c * r * sqrtq) / (q - one(T)))
    dhdθ = zero(T)
    if N > 0
        p = one(T) + k * c2 * r2
        if p < zero(T)
            return SVector{3,T}(NaN, NaN, NaN), SVector{3,T}(NaN, NaN, NaN)
        end
        sqrtp = sqrt(p)
        a = sqrtp / sqrtq
        dadρ = (c2 * (2k + one(T)) * rad * r) / (sqrtp * sqrtq^3)
        n = rad / z.normradius
        u = ρ * n
        u2 = u^2
        S0v = QType.S0(z.b0coeffs, u2)::T
        b = u2 * (one(T) - u2) * S0v
        dbdu = 2 * u * (one(T) - 2 * u2) * S0v + 2 * u^3 * (one(T) - u2) * QType.dS0dx(z.b0coeffs, u2)::T
        c = zero(T)
        dcdu = zero(T)
        if M > 0
            @inbounds @simd for m in 1:M
                Sα = QType.S(z.dαcoeffs[m], m, u2)::T
                Sβ = QType.S(z.dβcoeffs[m], m, u2)::T
                sinmθ = sin(m * θ)
                cosmθ = cos(m * θ)
                dhdθ += u^m * m * (-sinmθ * Sα + cosmθ * Sβ)
                c += u^m * (cosmθ * Sα + sinmθ * Sβ)
                dcdu += u^(m - 1) * (m * (cosmθ * Sα + sinmθ * Sβ) + 2 * u2 * (cosmθ * QType.dSdx(z.dαcoeffs[m], m, u2)::T + sinmθ * QType.dSdx(z.dβcoeffs[m], m, u2)::T))
            end
        end
        dhdθ *= a
        dhdρ += dadρ * (b + c) + a * n * (dbdu + dcdu) # want the derivative wrt ρ, not u
    end
    cosθ = cos(θ)
    sinθ = sin(θ)
    pu = SVector{3,T}(rad * cosθ, rad * sinθ, dhdρ)
    pv = SVector{3,T}(r * -sinθ, r * cosθ, dhdθ)
    return pu, pv
end

function normal(z::QTypeSurface{T,D,M,N}, ρ::T, θ::T)::SVector{3,T} where {T<:Real,D,M,N}
    du, dv = partials(z, ρ, θ)
    if ρ == zero(T) && norm(dv) == zero(T)
        # in cases where there is no δθ at ρ = 0 (i.e. anything which is rotationally symetric)
        # then we get some big problems, hardcoding this case solves the problems
        return SVector{3,T}(0, 0, 1)
    end
    return normalize(cross(du, dv))
end

function uv(z::QTypeSurface{T,3,M,N}, p::SVector{3,T})::SVector{2,T} where {T<:Real,M,N}
    # avoid divide by zero for ForwardDiff
    ϕ = NaNsafeatan(p[2], p[1])
    if p[1] == zero(T) && p[2] == zero(T)
        ρ = zero(T)
    else
        ρ = sqrt(p[1]^2 + p[2]^2) / semidiameter(z)
    end
    return SVector{2,T}(ρ, ϕ)
end

function onsurface(surf::QTypeSurface{T,3,M,N}, p::SVector{3,T}) where {T<:Real,M,N}
    ρ, θ = uv(surf, p)
    if ρ > one(T)
        return false
    else
        surfpoint = point(surf, ρ, θ)
        return samepoint(p[3], surfpoint[3])
    end
end

function inside(surf::QTypeSurface{T,3,M,N}, p::SVector{3,T}) where {T<:Real,M,N}
    ρ, θ = uv(surf, p)
    if ρ > one(T)
        return false
    else
        surfpoint = point(surf, ρ, θ)
        return p[3] < surfpoint[3]
    end
end

#########################################################################################################

# Assumes the ray has been transformed into the canonical coordinate frame which has the vertical axis passing through (0,0,0) and aligned with the z axis.
function surfaceintersection(surf::AcceleratedParametricSurface{T,3,QTypeSurface{T,3,M,N}}, r::AbstractRay{T,3}) where {T<:Real,M,N}
    cylint = surfaceintersection(surf.surface.boundingcylinder, r)
    if cylint isa EmptyInterval{T}
        return EmptyInterval(T)
    else
        if doesintersect(surf.triangles_bbox, r) || inside(surf.triangles_bbox, origin(r))
            surfint = triangulatedintersection(surf, r)
            if !(surfint isa EmptyInterval{T})
                return intervalintersection(cylint, surfint)
            end
        end
        # hasn't hit the surface
        if lower(cylint) isa RayOrigin{T} && upper(cylint) isa Infinity{T}
            if inside(surf.surface, origin(r))
                return Interval(RayOrigin(T), Infinity(T))
            else
                return EmptyInterval(T)
            end
            # otherwise check that the intersection is underneath the surface
        else
            p = point(closestintersection(cylint, false))
            ρ, ϕ = uv(surf, p)
            surfpoint = point(surf.surface, ρ, ϕ)
            if p[3] < surfpoint[3]
                return cylint # TODO!! UV (and interface) issues?
            else
                return EmptyInterval(T)
            end
        end
    end
end

function AcceleratedParametricSurface(surf::S, numsamples::Int = 17; interface::NullOrFresnel{T} = NullInterface(T)) where {T<:Real,N,S<:QTypeSurface{T,N}}
    # Zernike users ρ, ϕ uv space so need to modify extension of triangulation
    a = AcceleratedParametricSurface(surf, triangulate(surf, numsamples, true, false, true, false), interface = interface)
    emptytrianglepool!(T)
    return a
end

function BoundingBox(surf::QTypeSurface{T,3}) where {T<:Real}
    xmin = -semidiameter(surf)
    xmax = semidiameter(surf)
    ymin = -semidiameter(surf)
    ymax = semidiameter(surf)
    # curvature only goes one way
    q = one(T) - (one(T) + surf.conic) * surf.curvature^2 * surf.semidiameter^2
    if q < zero(T)
        throw(ErrorException("The surface is invalid, no bounding box can be constructed"))
    end
    hmax = surf.curvature * surf.semidiameter^2 / (one(T) + sqrt(q))
    if hmax > zero(T)
        zmax = hmax + surf.maxheight
        zmin = -surf.maxheight
    else
        zmax = surf.maxheight
        zmin = hmax - surf.maxheight
    end
    return BoundingBox(xmin, xmax, ymin, ymax, zmin, zmax)
end
