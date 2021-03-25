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

# only returns real roots
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
