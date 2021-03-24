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
