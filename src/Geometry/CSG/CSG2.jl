# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

import Base: ∩, ∪, -
using LinearAlgebra

export CSGTree2

# TODO CSGComplement, CSGDifference
abstract type CSGType end
abstract type CSGLeaf <: CSGType end
abstract type CSGUnion <: CSGType end
abstract type CSGIntersection <: CSGType end

"""
    CSGData{T<:Real}

A flexible struct that can hold the data for any kind of CSG node. Leaf nodes define both a `transform` and `surface`;
non-leaf nodes leave the latter undefined. The inverse transform is calculated and stored in the `inverse` field.
"""
struct CSGData{T<:Real}
    transform::Transform{T}
    inverse::Transform{T}
    surface::ParametricSurface{T,3}

    # non-leaf node data
    CSGData(transform::Transform{T}) where T = new{T}(transform, inv(transform))
    CSGData{T}() where T = CSGData(identitytransform(T))

    # leaf node data
    function CSGData(surface::ParametricSurface{T,3}, transform::Transform{T} = identitytransform(T)) where T
        return new{T}(transform, inv(transform), surface)
    end
end
function Base.show(io::IO, data::CSGData)
    print(io, "transform = $(data.transform ≈ I ? "I" : data.transform)")
    isdefined(data, :surface) && print(io, ", surface = $(data.surface)")
end

"""
    CSGTree{N<:CSGType,T<:Real}

A binary tree representing a CSG structure. `N` defines the root node's type: Leaf, Union, Intersection or Difference.

The leaf nodes of a `CSGTree` contain parametric surfaces. Non-leaf nodes define CSG operations between their children.
At each node, it is possible to store a transform that propagates down the tree at evaluation.
"""
mutable struct CSGTree2{N<:CSGType,T<:Real} <: Primitive{T}
    data::CSGData{T}
    left::CSGTree2{<:CSGType,T}
    right::CSGTree2{<:CSGType,T}
    parent::CSGTree2{<:CSGType,T}

    function CSGTree2(surface::ParametricSurface{T,3}, transform::Transform{T} = identitytransform(T)) where T
        return new{CSGLeaf,T}(CSGData(surface, transform))
    end

    function CSGTree2{N}(
        left::CSGTree2{<:CSGType,T}, right::CSGTree2{<:CSGType,T}, transform::Transform{T} = identitytransform(T)
    ) where {N,T}
        return new{N,T}(CSGData(transform), left, right)
    end
end
include("abstracttrees.jl")

# create and return a parent node whilst simultaneously setting the parent field for the left and right child nodes
csg_op!(::Type{N}, l::CSGTree2, r::CSGTree2) where N = l.parent = r.parent = CSGTree2{N}(l, r)
∪(l::CSGTree2, r::CSGTree2) = csg_op!(CSGUnion, l, r)
∩(l::CSGTree2, r::CSGTree2) = csg_op!(CSGIntersection, l, r)

# find the global transform/inverse of a node by traversing the path up to the root node (recursive)
transform(node::CSGTree2) = isdefined(node, :parent) ? transform(node.parent) * node.data.transform : node.data.transform
inverse(node::CSGTree2) = isdefined(node, :parent) ? node.data.inverse * inverse(node.parent) : node.data.inverse

# apply a local transform value to a node
function transform!(node::CSGTree2, transform::Transform)
    node.data = CSGData(transform * node.data.transform)
    return node
end

BoundingBox(node::CSGTree2{CSGLeaf}) = BoundingBox(node.data.surface, transform(node))
# etc...
