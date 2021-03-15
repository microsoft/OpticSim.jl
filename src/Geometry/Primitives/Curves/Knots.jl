"""
    KnotVector{T<:Number}

Vector to define knots used for [`BSplineCurve`](@ref) and [`BSplineSurface`](@ref).
"""
struct KnotVector{T<:Number}
    knots::Array{T,1}

    function KnotVector{T}(a::Array{T,1}) where {T<:Number}
        temp = a[1]
        for i in 2:lastindex(a) # ensure knot vector is non-decreasing
            @assert a[i] >= temp
            temp = a[i]
        end
        return new(copy(a))
    end
end
export KnotVector

Base.copy(knots::KnotVector{S}) where {S} = KnotVector{S}(knots.knots)

urange(knots::KnotVector) = (knots.knots[1], knots.knots[end])
numknots(knots::KnotVector) = size(knots.knots)[1]

struct InvalidParameterValue <: Exception
    invalid_u::Any
    validrange::Any
end

Base.showerror(io::IO, e::InvalidParameterValue) = print(io, "invalid parameter value: $(e.invalid_u) allowable parameter range: $(e.validrange)")

function findspan(knots::KnotVector, curveorder, u)
    pmin = curveorder + 1
    nknots = numknots(knots)
    pmax = nknots - (curveorder + 1)
    if knots.knots[pmax + 1] == u
        return pmax
    end

    for i in pmin:pmax # if u == ui+1 then the span goes from i to i+1.
        if u <= knots.knots[i + 1]
            return i
        end
    end
    # this might be too extreme. round off error could cause parameter value to go a little above max value. Might need to loosen this later.
    throw(InvalidParameterValue(u, (knots.knots[1], knots.knots[end]))) # return the bad u value and the legal range of u values
end

function basisfunctions(knots::KnotVector{T}, u::T, curveorder) where {T<:Number}
    # computes values of p+1 basis functions where p is the curve order. knotsegment is the index of the knot vector segment, u is the value of the curve parameter
    knotsegment = findspan(knots, curveorder, u)
    U = knots.knots
    tempsize = curveorder + 1
    left = zeros(T, tempsize) # later rewrite this code to have a struct that stores these temps so they aren't being allocated on every curve evaluation
    right = zeros(T, tempsize)
    N = zeros(T, tempsize)

    N[1] = one(T)

    for j in 1:curveorder
        left[j + 1] = u - U[knotsegment + 1 - j]
        right[j + 1] = U[knotsegment + j] - u
        saved = zero(T)
        for r in 0:(j - 1)
            leftindex = j - r + 1
            temp = N[r + 1] / (right[r + 2] + left[leftindex])
            N[r + 1] = saved + right[r + 2] * temp
            saved = left[leftindex] * temp
        end
        N[j + 1] = saved
    end
    return N
end
