struct PowerBasisCurve{P,S,N,M} <: Spline{P,S,N,M}
    controlpolygon::Array{S,2} #not exactly the right name, should be coefficients instead, but consistent with the other splines
    #data stored like this:
    # a0x a1x ... aM+1x
    # a0y a1y ... aM+1x
    # and so on for however many spatial dimensions

    function PowerBasisCurve{P,S,N,M}(coefficients::Array{S,2}) where {P,S,N,M}
        @assert (N, M + 1) == size(coefficients)

        return new{P,S,N,M}(copy(coefficients))
    end
end
export PowerBasisCurve

function coefficients(curve::PowerBasisCurve{P,S,N,M}, spatialindex) where {P,S,N,M}
    return curve.controlpolygon[spatialindex, :]
end

function PowerBasisCurve{P,S,N,M}(curve::BezierCurve{P,S,N,M}) where {P,S,N,M}
    controlpolygon = curve.controlpolygon
    coefficients = Array{S,2}(undef, N, M + 1) # array stores coefficients for all N spatial dimensions

    #this is probably not the most efficient way to convert from Bernstein to Power basis. Optimize later if necessary.
    for k in 0:M
        for i in k:M
            @. coefficients[:, i + 1] += (-1)^(i - k) * binomial(M, i) * binomial(i, k) * controlpolygon[k + 1]
        end
    end
    return PowerBasisCurve{P,S,N,M}(coefficients)
end

function beziertopowerbasis(k, n)
    coefficient = 0
    for i in k:n
        @. coefficients[:, i + 1] += (-1)^(i - k) * binomial(n, i) * binomial(i, k)
    end
end

function point(curve::PowerBasisCurve, u::T) where {T<:Number}
    _, n = size(curve.controlpolygon)
    sum = curve.controlpolygon[:, n]

    for i in (n - 1):-1:1
        @. sum = u * sum + curve.controlpolygon[:, i]
    end

    return sum
end
