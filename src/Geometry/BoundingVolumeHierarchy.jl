# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

#=
Notes for CSG bounding boxes
    The intersection nodes could compute bounding boxes as the CSGTree is being constructed because the bbox can only be equal or smaller than either of the children. Subtrees of union nodes should be processed as a whole so they can be properly subdivided or the bounding boxes are likely to be very inefficient.
    Need to cluster nodes geometrically to subdivide well but the order in which the unions are called by the user's code is unlikely to have the best geometric distribution.
    This makes the recursion a bit more complicated than conventional BVH because of the interleaved union and intersection operations. 
=#

"""if fewer than this many primitives in a bounding box then do no subdivide further"""
const minprims = 4

abstract type BVHNode{T<:Real,S<:Primitive{T}} end

struct SplitBVH{T<:Real,S<:Primitive{T}} <: BVHNode{T,S}
    left::BVHNode{T,S}
    right::BVHNode{T,S}
    axis::Int64
end

struct LeafBVH{T<:Real,S<:Primitive{T}} <: BVHNode{T,S}
    primitives::Vector{S}
end

struct BVHData{T<:Real,S<:Primitive{T}}
    primitive::S
    centroid::SVector{3,T}
    boundingbox::BoundingBox{T}
end

centroid(a::BVHData) = a.centroid
boundingbox(a::BVHData) = a.boundingbox
primitive(a::BVHData) = a.primitive

struct PrimitiveData{T<:Real,S<:Primitive{T}}
    data::Vector{BVHData{T,S}}
    boundingbox::BoundingBox{T}

    function PrimitiveData(primitives::Vector{S}) where {T<:Real,S<:Primitive{T}}
        combined = Vector{BVHData{T,S}}(undef, length(primitives))
        totalbox = nothing

        for i in eachindex(primitives)
            box = BoundingBox(primitives[i])
            combined[i] = BVHData{T,S}(primitives[i], centroid(primitives[i]), box)
            totalbox = union(totalbox, box)
        end

        return new{T,S}(combined, totalbox)
    end
end

function computebvh(a::PrimitiveData)
    subdivide(a.boundingbox, a)

end

function subdivide(parentbounds::BoundingBox{T}, primitives::AbstractVector{S}) where {T<:Real,S<:Primitive{T}}
    if length(primitives <= minprims)
        return LeafBVH(primitives)
    else
        # test subdivision in 3 axes.
        #partition primitives
        # compute costs using parentbounds
        # return SplitBVH(subdivide(reduce(union,leftprims),leftprims),subdivide(reduce(union,rightprims,rightprims))
    end
end

cost(enclosingbox::BoundingBox{T}, childboxes::AbstractVector{BoundingBox{T}}) where {T<:Real} = sum(area.(childboxes)) / area(enclosingbox)

cost(boxes::Vector{BVHData{T,S}}) where {T,S} = sum(area.(boundingbox.(boxes)))

function partition!(a::Vector{S}, valfunc::Function, value::T) where {S,T<:Real}
    lower = 0
    upper = lastindex(a) + 1

    while true
        while upper - 1 >= 1 && valfunc(a[upper - 1]) > value
            upper = upper - 1
        end
        while lower + 1 <= lastindex(a) && valfunc(a[lower + 1]) <= value
            lower = lower + 1
        end

        if upper - lower == 1
            if lower == 0
                return (nothing, a)
            else
                if upper == lastindex(a) + 1
                    return (a, nothing)
                else
                    return (view(a, 1:lower), view(a, upper:lastindex(a)))
                end
            end
        end

        temp = a[upper - 1]
        a[upper - 1] = a[lower + 1]
        a[lower + 1] = temp
        upper = upper - 1
        lower = lower + 1
    end
end

"""compute compacted version of bounding volume hierarchy that can be traversed more quickly than the original tree structure, without recursion"""
function linearizebvh(bvhtree::BVHNode{T,S}) where {T<:Real,S<:Primitive{T}} end

"""traverses BVH and computes surface intersection"""
function surfaceintersection(a::BVHNode{T,S}, r::Ray{T,N}) where {T<:Real,S<:Primitive{T},N} end

# function testpartition()
#     split = .5

#     for i in 1:10000000
#         a = rand(5)
#         lower,upper = partition!(a,(x)->x,split)

#         if lower !== nothing
#             for i in lower
#                 @assert i < split
#             end
#         end

#         if upper !== nothing
#             for i in upper
#                 @assert i > split
#             end
#         end
#     end
# end

export PrimitiveData
