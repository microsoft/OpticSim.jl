# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.


"""
    plane_from_points(points::SMatrix{3, N, Float64}}) ->  centroid, normal

Estimate the best fitting plane for a set of points in 3D.
A more efficient version of plane_from_points.
"""
function plane_from_points(points::SMatrix{3, N, Float64} where {N}) 
    center = Statistics.mean(points,dims=2)

    u, _, _ = svd(points .- center)
    normal = u[:,3]             # singular vectors in decending order

    # make sure the normal is pointing consistently to positive Z direction 
    if (dot(normal, unitZ3()) < 0.0)
        normal = normal * -1.0
    end

    return SVector(center), SVector(normal)     # convert from SMatrix to SVector
end

"""only returns real roots"""
function quadraticroots(a::T, b::T, c::T) where {T<:Real}
    temp = b^2 - 4 * a * c

    if temp < zero(T)
        return nothing # no real roots so no ray cylinder intersection
    end

    radical = sqrt(temp)
    if b >= zero(T)
        radicalterm = -b - radical
        x1 = radicalterm / (2a)
        x2 = (2c) / radicalterm
    else
        radicalterm = -b + radical
        x1 = (2c) / radicalterm
        x2 = radicalterm / (2a)
    end

    return SVector{2,T}(x1, x2)
end

replprint(a) = show(IOContext(stdout), "text/plain", a)

@inline samepoint(pt1, pt2) = isapprox(pt1, pt2, atol = 1e-10)

password() = String(collect(rand(('1':'9'..., 'A':'Z'..., 'a':'z'...)) for i in 1:14))

# 1 if val>= 0 and -1 otherwise. Unlike Math.Sign which returns 0 if val==0

@inline sign(val::T) where {T<:Real} = val >= zero(T) ? 1 : -1

@inline function NaNsafeatan(x::T, y::T)::T where {T<:Real}
    if y == zero(T)
        if x == zero(T)
            return zero(T)
        elseif x > zero(T)
            return T(π / 2)
        else
            return -T(π / 2)
        end
    else
        return atan(x, y)
    end
end

@inline function NaNsafeasin(x::T)::T where {T<:Real}
    if x == zero(T)
        return zero(T)
    elseif x == one(T)
        return T(π / 2)
    elseif x == -one(T)
        return -T(π / 2)
    else
        return asin(x)
    end
end

@inline function NaNsafeacos(x::T)::T where {T<:Real}
    if x == zero(T)
        return T(π / 2)
    elseif x == one(T)
        return zero(T)
    elseif x == -one(T)
        return T(π)
    else
        return acos(x)
    end
end


# some place holders for package level function names.
# these names need to exist before any internal module can override them.
function origin end    
