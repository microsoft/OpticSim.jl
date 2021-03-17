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

import Base: setindex!, getindex

"""
    HierarchicalImage{T<:Number} <: AbstractArray{T,2}

Image type which dynamically allocated memory for pixels when their value is set, the value of unset pixels is assumed to be zero.

This is used for the detector image of [`OpticalSystem`](@ref)s which can typically be very high resolution, but often have a large proportion of the image blank.
"""
struct HierarchicalImage{T<:Number} <: AbstractArray{T,2} #allow for complex values. May want to expand this even further to allow more complicated detector image types
    toplevel::Array{Array{T,2},2}
    secondlevelrows::Int64
    secondlevelcols::Int64
    apparentrows::Int64
    apparentcols::Int64

    function HierarchicalImage{T}(height, width) where {T<:Number}
        heightfirst, heightsecond = nearestsqrts(height)
        widthfirst, widthsecond = nearestsqrts(width)
        return HierarchicalImage{T}(height, width, heightfirst, widthfirst, heightsecond, widthsecond)
    end
    HierarchicalImage{T}(apparentheight, apparentwidth, height, width, secondlevelheight, secondlevelwidth) where {T<:Number} = new{T}(Array{Array{T,2},2}(undef, height, width), secondlevelheight, secondlevelwidth, apparentheight, apparentwidth)
end
export HierarchicalImage

function nearestsqrts(num)
    rt1 = Int64(ceil(sqrt(num)))
    rt2 = Int64(floor(sqrt(num)))
    if rt1 * rt2 >= num
        return (rt1, rt2)
    else
        return (rt1, rt2 + 1)
    end
end

toplevel(a::HierarchicalImage) = a.toplevel
secondlevelrows(a::HierarchicalImage) = a.secondlevelrows
secondlevelcols(a::HierarchicalImage) = a.secondlevelcols
apparentrows(a::HierarchicalImage) = a.apparentrows
apparentcols(a::HierarchicalImage) = a.apparentcols

topindex(a::HierarchicalImage, index1, index2)::Int = linearindex(size(toplevel(a), 1), 1 + (index1 - 1) ÷ secondlevelrows(a), 1 + (index2 - 1) ÷ secondlevelcols(a))
lowerindex(a::HierarchicalImage, index1, index2)::Int = linearindex(secondlevelrows(a), 1 + mod(index1 - 1, secondlevelrows(a)), 1 + mod(index2 - 1, secondlevelcols(a)))

linearindex(rows::Int, index1::Int, index2::Int)::Int = index1 + rows * (index2 - 1)

function outofbounds(a::HierarchicalImage, indices::Tuple{Int,Int})
    if any(indices .> size(a))
        throw(BoundsError(a, indices))
    end
end

Base.eltype(::HierarchicalImage{T}) where {T<:Number} = T
Base.length(a::HierarchicalImage{T}) where {T<:Number} = reduce(*, size(a))
Base.size(a::HierarchicalImage{T}) where {T<:Number} = (apparentrows(a), apparentcols(a))

function Base.getindex(a::HierarchicalImage{T}, indices::Vararg{Int,2}) where {T<:Number}
    outofbounds(a, indices)

    top = topindex(a, indices[1], indices[2])
    bott = lowerindex(a, indices[1], indices[2])

    return getindexdirect(a, top, bott)
end

function getindexdirect(a::HierarchicalImage{T}, topidx::Int, botidx::Int) where {T<:Real}
    if !isassigned(toplevel(a), topidx) #by default any part of the image that was not written to has the value zero
        return zero(T)
    else
        return toplevel(a)[topidx][botidx]
    end
end

function Base.setindex!(a::HierarchicalImage, v::T, indices::Vararg{Int,2}) where {T<:Number} #if use HierarchicalImage{T} get error "setindex! not defined". No idea why this happens
    outofbounds(a, indices)

    top = topindex(a, indices[1], indices[2])
    bott = lowerindex(a, indices[1], indices[2])

    setindexdirect!(a, v, top, bott)
end

function setindexdirect!(a::HierarchicalImage{T}, v::T, topidx::Int, botidx::Int) where {T<:Real}
    if !isassigned(toplevel(a), LinearIndices(toplevel(a))[topidx]) #by default any part of the image that has not been not written to is set to zero
        toplevel(a)[topidx] = fill(zero(T), secondlevelrows(a), secondlevelcols(a))
    end
    toplevel(a)[topidx][botidx] = v
end

"""
    reset!(a::HierarchicalImage{T})

Resets the pixels in the image to `zero(T)`.
Do this rather than `image .= zero(T)` because that will cause every pixel to be accessed, and therefore allocated.
For large images this can cause huge memory traffic.
"""
function reset!(a::HierarchicalImage{T}) where {T}
    for index in 1:length(toplevel(a))
        if isassigned(toplevel(a), index)
            temp = toplevel(a)[index]
            for k in 1:length(temp)
                temp[k] = zero(T)
            end
        end
    end
end

"""
    sum!(a::HierarchicalImage{T}, b::HierarchicalImage{T})

Add the contents of `b` to `a` in an efficient way.
"""
function sum!(a::HierarchicalImage{T}, b::HierarchicalImage{T}) where {T}
    @assert size(a) == size(b)
    for index in 1:length(toplevel(b))
        if isassigned(toplevel(b), index)
            temp = toplevel(b)[index]
            for k in 1:length(temp)
                setindexdirect!(a, getindexdirect(a, index, k) + temp[k], index, k)
            end
        end
    end
end

export sum!, reset!
