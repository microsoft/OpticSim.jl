# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

# optional add-in for AbstractTrees support
# https://github.com/JuliaCollections/AbstractTrees.jl/blob/master/examples/binarytree_infer.jl
using AbstractTrees

# Implement iteration over the immediate children of a node
function Base.iterate(node::CSGTree2)
    isdefined(node, :left) && return (node.left, false)
    isdefined(node, :right) && return (node.right, true)
    return nothing
end
function Base.iterate(node::CSGTree2, state::Bool)
    state && return nothing
    isdefined(node, :right) && return (node.right, true)
    return nothing
end
Base.IteratorSize(::Type{CSGTree2}) = Base.SizeUnknown()
Base.eltype(::Type{CSGTree2{N,T}}) where {N,T} = CSGTree2{N,T}

## Things we need to define to leverage the native iterator over children for the purposes of AbstractTrees.
# Set the traits of this kind of tree
Base.eltype(::Type{<:TreeIterator{CSGTree2{N,T}}}) where {N,T} = CSGTree2{N,T}
Base.IteratorEltype(::Type{<:TreeIterator{CSGTree2}}) = Base.HasEltype()
AbstractTrees.parentlinks(::Type{CSGTree2}) = AbstractTrees.StoredParents()
AbstractTrees.siblinglinks(::Type{CSGTree2}) = AbstractTrees.StoredSiblings()
# Use the native iteration for the children
AbstractTrees.children(node::CSGTree2) = node

Base.parent(::CSGTree2, node::CSGTree2) = isdefined(node, :parent) ? node.parent : nothing

function AbstractTrees.nextsibling(::CSGTree2, child::CSGTree2)
    isdefined(child, :parent) || return nothing
    p = child.parent
    if isdefined(p, :right)
        child === p.right && return nothing
        return p.right
    end
    return nothing
end

# We also need `pairs` to return something sensible.
# If you don't like integer keys, you could do, e.g.,
#   Base.pairs(node::BinaryNode) = BinaryNodePairs(node)
# and have its iteration return, e.g., `:left=>node.left` and `:right=>node.right` when defined.
# But the following is easy:
Base.pairs(node::CSGTree2) = enumerate(node)

AbstractTrees.printnode(io::IO, node::CSGTree2{N,T}) where {N,T} = print(io, "$N: $(node.data)")
