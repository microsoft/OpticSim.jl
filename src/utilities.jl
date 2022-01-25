# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.


"""
    plane_from_points(points::SMatrix{D, N, P}}) ->  centroid, normal, local_to_world transform

Points to be fitted are assumed to be stored by column in the `points` matrix.
Estimate the best fitting plane for a set of points in 3D.
`D` is the dimension of the plane.
`N` is the number of points to fit.
`P` is the number type used to represent points.
"""
function plane_from_points(points::SMatrix{D, N, P}) where {D,N,P<:Real}
    center = Statistics.mean(points,dims=2) #compute average of columns of point matrix.

    u, _, _ = svd(points .- center)
    
    #always want a rotation matrix.
    if det(u) < 0
        u = SMatrix{D,D}(u[:,1:(end-1)]...,-u[:,end]...) #change sign of last singular vector to turn into rotation matrix
    end

    normal = u[:,3]             # The two largest singular vectors lie in the plane that is the best fit to the points, i.e., that accounts for the largest fraction of variance in the set of points. The smallest singular vector is perpendicular to this plane.

    # make sure the normal is pointing consistently to positive Z direction of local coordinate frame
    if dot(normal, unitZ3()) < P(0.0)
        normal = normal * P(-1.0)
    end

    return SVector(center), SVector(normal), u     # convert from SMatrix to SVector, return u matrix to use as local coordinate frame for the plane.
end

function plane_from_points(points::AbstractMatrix{P}) where{P<:Real}
    dims = size(points)
    temp = MMatrix{dims[1],dims[2],P}(points)
    return plane_from_points(SMatrix(temp))
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
