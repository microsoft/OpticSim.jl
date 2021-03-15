"""Abstract type representing any evaluated CSG structure."""
abstract type CSGTree{T} <: Primitive{T} end
export CSGTree

"""
    ComplementNode{T,C<:CSGTree{T}} <: CSGTree{T}

An evaluated complement node within the CSG tree, must be the second child of a [`IntersectionNode`](@ref) forming a subtraction.
"""
struct ComplementNode{T,C<:CSGTree{T}} <: CSGTree{T}
    child::C

    function ComplementNode(child::C) where {T<:Real,C<:CSGTree{T}}
        return new{T,C}(child)
    end
end
Base.show(io::IO, a::ComplementNode{T}) where {T} = print(io, "Complement($(a.child))")
BoundingBox(a::ComplementNode{T}) where {T<:Real} = BoundingBox(a.child)

"""
    UnionNode{T,L<:CSGTree{T},R<:CSGTree{T}} <: CSGTree{T}

An evaluated union node within the CSG tree.
"""
struct UnionNode{T,L<:CSGTree{T},R<:CSGTree{T}} <: CSGTree{T}
    leftchild::L
    rightchild::R
    bbox::BoundingBox{T}

    function UnionNode(a::L, b::R) where {T<:Real,L<:CSGTree{T},R<:CSGTree{T}}
        # union should never contain a complement so should be fine
        return new{T,L,R}(a, b, union(BoundingBox(a), BoundingBox(b)))
    end
end
Base.show(io::IO, a::UnionNode{T}) where {T} = print(io, "Union($(a.leftchild), $(a.rightchild))")

"""
    IntersectionNode{T,L<:CSGTree{T},R<:CSGTree{T}} <: CSGTree{T}

An evaluated intersection node within the CSG tree.
"""
struct IntersectionNode{T,L<:CSGTree{T},R<:CSGTree{T}} <: CSGTree{T}
    leftchild::L
    rightchild::R
    bbox::BoundingBox{T}

    function IntersectionNode(a::L, b::R) where {T<:Real,L<:CSGTree{T},R<:CSGTree{T}}
        # normal intersection is fine for most nodes
        return new{T,L,R}(a, b, intersection(BoundingBox(a), BoundingBox(b)))
    end

    function IntersectionNode(a::L, b::R) where {T<:Real,L<:CSGTree{T},R<:ComplementNode{T}}
        # this is a subtraction so just take the original unclipped bounding box for simplicity
        return new{T,L,R}(a, b, BoundingBox(a))
    end
end
Base.show(io::IO, a::IntersectionNode{T}) where {T} = print(io, "Intersection($(a.leftchild), $(a.rightchild))")

"""
    LeafNode{T,S<:ParametricSurface{T}} <: CSGTree{T}

An evaluated leaf node in the CSG tree, `geometry` attribute which contains a [`ParametricSurface`](@ref) of type `S`.
The leaf node also has a transform associated which is the composition of all nodes above it in the tree.
As such, transforming points from the geometry using this transform puts them in world space, and transforming rays by the inverse transform puts them in object space.
"""
struct LeafNode{T,S<:ParametricSurface{T,3}} <: CSGTree{T}
    geometry::S
    transform::RigidBodyTransform{T}
    invtransform::RigidBodyTransform{T}
    bbox::BoundingBox{T}

    function LeafNode(a::S, transform::RigidBodyTransform{T}) where {T<:Real,S<:ParametricSurface{T,3}}
        # store the transformed bounding box so nodes higher in the tree have correct global space bounding boxes
        return new{T,S}(a, transform, inv(transform), BoundingBox(a, transform))
    end
end
Base.show(io::IO, a::LeafNode{T}) where {T} = a.transform != identitytransform(T) ? print(io, "Leaf($(a.geometry), $(a.transform))") : print(io, "Leaf($(a.geometry))")

BoundingBox(a::CSGTree{T}) where {T<:Real} = a.bbox

"""
    CSGGenerator{T<:Real}

This is the type you should use when making CSG objects.
This type allows for the construction of [`CSGTree`](@ref) objects with different transforms.
When the generator is evaluated, all transforms are propagated down to the [`LeafNode`](@ref)s and stored there.

# Example
```julia
a = Cylinder(1.0,1.0)
b = Plane([0.0,0.0,1.0], [0.0,0.0,0.0])
generator = csgintersection(a,b)
# now make a csg object that can be ray traced
csgobj = generator(RigidBodyTransform(1.0,1.0,2.0))
```
"""
struct CSGGenerator{T<:Real}
    f::Function
end
export CSGGenerator

function (a::CSGGenerator{T})(transform::RigidBodyTransform{T})::CSGTree{T} where {T<:Real}
    a.f(transform)
end
function (a::CSGGenerator{T})()::CSGTree{T} where {T<:Real}
    a.f(identitytransform(T))
end

"""
    leaf(surf::ParametricSurface{T}, transform::RigidBodyTransform{T} = identitytransform(T)) -> CSGGenerator{T}

Create a leaf node from a parametric surface with a given transform.
"""
function leaf(surf::S, transform::RigidBodyTransform{T} = identitytransform(T))::CSGGenerator{T} where {T<:Real,S<:ParametricSurface{T}}
    return CSGGenerator{T}((parenttransform) -> LeafNode(surf, parenttransform * transform))
end
"""
    leaf(surf::CSGGenerator{T}, transform::RigidBodyTransform{T} = identitytransform(T)) -> CSGGenerator{T}

Create a (pseudo) leaf node from another CSGGenerator, this is useful if you want multiple copies of a premade CSG structure with different transforms, for example in an MLA.
"""
function leaf(n::CSGGenerator{T}, transform::RigidBodyTransform{T} = identitytransform(T))::CSGGenerator{T} where {T<:Real}
    return CSGGenerator{T}((parenttransform) -> n(parenttransform * transform))
end
export leaf

# CSG objects are trees, not graphs, since it makes no sense to reuse intermediate results. Hence no need to keep track of common subexpressions.
# Static arrays are much faster than mutable arrays so want node transforms to be static. Unfortunately the transformations cascade from the root to the leaves but the nodes have to be created from the leaves to the roots.
# Wrap the csg operations in function that delays evaluation until the transform has been computed.

"""
    csgintersection(a::CSGGenerator{T} b::CSGGenerator{T}, transform::RigidBodyTransform{T} = identitytransform(T)) -> CSGGenerator{T}

Create a binary node in the CSG tree representing an intersection between a and b.
A shortcut method for `a` and `b` as [`ParametricSurface`](@ref)s is also available.

![Intersect Image](https://upload.wikimedia.org/wikipedia/commons/0/0b/Boolean_intersect.PNG)
"""
csgintersection(a::CSGGenerator{T}, b::CSGGenerator{T}, transform::RigidBodyTransform{T} = identitytransform(T)) where {T<:Real} = CSGGenerator{T}((parenttransform) -> IntersectionNode(a(parenttransform * transform), b(parenttransform * transform)))
csgintersection(a::ParametricSurface{T}, b::ParametricSurface{T}, transform::RigidBodyTransform{T} = identitytransform(T)) where {T<:Real} = csgintersection(leaf(a), leaf(b), transform)
csgintersection(a::ParametricSurface{T}, b::CSGGenerator{T}, transform::RigidBodyTransform{T} = identitytransform(T)) where {T<:Real} = csgintersection(leaf(a), b, transform)
csgintersection(a::CSGGenerator{T}, b::ParametricSurface{T}, transform::RigidBodyTransform{T} = identitytransform(T)) where {T<:Real} = csgintersection(a, leaf(b), transform)

export csgintersection

"""
    csgunion(a::CSGGenerator{T}, b::CSGGenerator{T}, transform::RigidBodyTransform{T} = identitytransform(T)) -> CSGGenerator{T}

Create a binary node in the CSG tree representing a union between a and b.
A shortcut method for `a` and `b` as [`ParametricSurface`](@ref)s is also available.

![Union Image](https://upload.wikimedia.org/wikipedia/commons/4/4a/Boolean_union.PNG)
"""
csgunion(a::CSGGenerator{T}, b::CSGGenerator{T}, transform::RigidBodyTransform{T} = identitytransform(T)) where {T<:Real} = CSGGenerator{T}((parenttransform) -> UnionNode(a(parenttransform * transform), b(parenttransform * transform)))
csgunion(a::ParametricSurface{T}, b::ParametricSurface{T}, transform::RigidBodyTransform{T} = identitytransform(T)) where {T<:Real} = csgunion(leaf(a), leaf(b), transform)
csgunion(a::ParametricSurface{T}, b::CSGGenerator{T}, transform::RigidBodyTransform{T} = identitytransform(T)) where {T<:Real} = ccsgunion(leaf(a), b, transform)
csgunion(a::CSGGenerator{T}, b::ParametricSurface{T}, transform::RigidBodyTransform{T} = identitytransform(T)) where {T<:Real} = csgunion(a, leaf(b), transform)

export csgunion

"""
    csgdifference(a::CSGGenerator{T}, b::CSGGenerator{T}, transform::RigidBodyTransform{T} = identitytransform(T)) -> CSGGenerator{T}

Create a binary node in the CSG tree representing the difference of a and b, essentially a - b.
A shortcut method for `a` and `b` as [`ParametricSurface`](@ref)s is also available.

![Difference Image](https://upload.wikimedia.org/wikipedia/commons/8/86/Boolean_difference.PNG)
"""
csgdifference(a::CSGGenerator{T}, b::CSGGenerator{T}, transform::RigidBodyTransform{T} = identitytransform(T)) where {T<:Real} = CSGGenerator{T}((parenttransform) -> IntersectionNode(a(parenttransform * transform), ComplementNode(b(parenttransform * transform))))
csgdifference(a::CSGGenerator{T}, b::ParametricSurface{T}, transform::RigidBodyTransform{T} = identitytransform(T)) where {T<:Real} = csgdifference(a, leaf(b), transform)
csgdifference(a::ParametricSurface{T}, b::CSGGenerator{T}, transform::RigidBodyTransform{T} = identitytransform(T)) where {T<:Real} = csgdifference(leaf(a), b, transform)
csgdifference(a::ParametricSurface{T}, b::ParametricSurface{T}, transform::RigidBodyTransform{T} = identitytransform(T)) where {T<:Real} = csgdifference(leaf(a), leaf(b), transform)

export csgdifference

function evalcsg(a::UnionNode{T}, ray::AbstractRay{T,N}, normalreverse::Bool = false)::Union{EmptyInterval{T},DisjointUnion{T},Interval{T}} where {T<:Real,N}
    if !doesintersect(a.bbox, ray)
        return EmptyInterval(T)
    end
    l = evalcsg(a.leftchild, ray, normalreverse)
    if l isa Interval{T}
        if lower(l) isa RayOrigin{T} && upper(l) isa Infinity{T}
            return l
        end
    end
    r = evalcsg(a.rightchild, ray, normalreverse)
    if l isa EmptyInterval{T}
        return r
    end
    if r isa Interval{T}
        if lower(r) isa RayOrigin{T} && upper(r) isa Infinity{T}
            return r
        end
    end
    if r isa EmptyInterval{T}
        return l
    end
    return intervalunion(l, r)
end

function evalcsg(a::IntersectionNode{T}, ray::AbstractRay{T,N}, normalreverse::Bool = false)::Union{EmptyInterval{T},DisjointUnion{T},Interval{T}} where {T<:Real,N}
    if !doesintersect(a.bbox, ray)
        return EmptyInterval(T)
    end
    l = evalcsg(a.leftchild, ray, normalreverse)
    if l isa EmptyInterval{T}
        return EmptyInterval(T)
    else
        r = evalcsg(a.rightchild, ray, normalreverse)
        if r isa EmptyInterval{T}
            return EmptyInterval(T)
        else
            return intervalintersection(l, r)
        end
    end
end

function evalcsg(a::ComplementNode{T}, ray::AbstractRay{T,N}, normalreverse::Bool = false)::Union{EmptyInterval{T},DisjointUnion{T},Interval{T}} where {T<:Real,N}
    return intervalcomplement(evalcsg(a.child, ray, !normalreverse))
end

function evalcsg(a::LeafNode{T}, ray::AbstractRay{T,N}, normalreverse::Bool = false)::Union{EmptyInterval{T},DisjointUnion{T},Interval{T}} where {T<:Real,N}
    # the bounding box is in global coordinates so use the un-transformed ray
    if !doesintersect(a.bbox, ray)
        return EmptyInterval(T)
    end
    intscts = surfaceintersection(a.geometry, a.invtransform * ray)
    if !(intscts isa EmptyInterval{T})
        if normalreverse
            intscts = reversenormal(intscts)
        end
        if intscts isa Interval{T}
            return a.transform * intscts
        elseif intscts isa DisjointUnion{T}
            for i in eachindex(intervals(intscts))
                intervals(intscts)[i] = a.transform * intervals(intscts)[i]
            end
            return intscts
        else
            throw(ErrorException("EvalCSG error - returned type not an interval or disjoint union"))
        end
    else
        return EmptyInterval(T)
    end
end

"""
    surfaceintersection(obj::CSGTree{T}, r::AbstractRay{T,N})

Calculates the intersection of `r` with CSG object, `obj`.

Returns an [`EmptyInterval`](@ref) if there is no intersection, an [`Interval`](@ref) if there is one or two interesections and a [`DisjointUnion`](@ref) if there are more than two intersections.

The ray is intersected with the [`LeafNode`](@ref)s that make up the CSG object and the resulting `Interval`s and `DisjointUnion`s are composed with the same boolean operations to give a final result.
The ray is transformed by the inverse of the transform associated with the leaf node to put it in _object space_ for that node before the intersection is carried out, typically this _object space_ is centered at the origin, but may differ for each primitive.

Some intersections are culled without actually evaluating them by first checking if the ray intersects the [`BoundingBox`](@ref) of each node in the [`CSGTree`](@ref), this can substantially improve performance in some cases.
"""
surfaceintersection(element::CSGTree{T}, r::AbstractRay{T,N}) where {T<:Real,N} = evalcsg(element, r)

point(::EmptyInterval) = nothing
normal(::EmptyInterval) = nothing
point(a::DisjointUnion) = point(closestintersection(a))
normal(a::DisjointUnion) = normal(closestintersection(a))
point(a::Interval) = point(closestintersection(a))
normal(a::Interval) = point(closestintersection(a))

"""
    onsurface(obj::CSGTree{T}, point::SVector{3,T}) -> Bool
    onsurface(obj::CSGTree{T}, x::T, y::T, z::T) -> Bool

Tests whether a 3D point in world space is _on_ the surface (i.e. shell) of `obj`.
"""
onsurface(a::CSGTree{T}, x::T, y::T, z::T) where {T<:Real} = onsurface(a, SVector{3,T}(x, y, z))
onsurface(a::ComplementNode{T}, point::SVector{3,T}) where {T<:Real} = onsurface(a.child, point)
onsurface(a::IntersectionNode{T}, point::SVector{3,T}) where {T<:Real} = (onsurface(a.leftchild, point) && inside(a.rightchild, point)) || (onsurface(a.rightchild, point) && inside(a.leftchild, point))
onsurface(a::UnionNode{T}, point::SVector{3,T}) where {T<:Real} = (onsurface(a.leftchild, point) && !inside(a.rightchild, point)) || (onsurface(a.rightchild, point) && !inside(a.leftchild, point))
onsurface(a::LeafNode{T}, point::SVector{3,T}) where {T<:Real} = onsurface(a.geometry, a.invtransform * point)

"""
    inside(obj::CSGTree{T}, point::SVector{3,T}) -> Bool
    inside(obj::CSGTree{T}, x::T, y::T, z::T) -> Bool

Tests whether a 3D point in world space is _inside_ `obj`.
"""
inside(a::CSGTree{T}, x::T, y::T, z::T) where {T<:Real} = inside(a, SVector{3,T}(x, y, z))
inside(a::ComplementNode{T}, point::SVector{3,T}) where {T<:Real} = !inside(a.child, point)
inside(a::IntersectionNode{T}, point::SVector{3,T}) where {T<:Real} = inside(a.bbox, point) && (inside(a.leftchild, point) && inside(a.rightchild, point))
inside(a::UnionNode{T}, point::SVector{3,T}) where {T<:Real} = inside(a.bbox, point) && (inside(a.leftchild, point) || inside(a.rightchild, point))
inside(a::LeafNode{T}, point::SVector{3,T}) where {T<:Real} = inside(a.bbox, point) && inside(a.geometry, a.invtransform * point)

########################################################################################################################

struct TrianglePool{T<:Real}
    allocated::Vector{Vector{Triangle{T}}}
    unallocated::Vector{Vector{Triangle{T}}}

    TrianglePool{T}() where {T<:Real} = new{T}(Vector{Vector{Triangle{T}}}(), Vector{Vector{Triangle{T}}}())
end

function allocate!(a::TrianglePool{T}) where {T<:Real}
    if length(a.unallocated) === 0
        push!(a.unallocated, Vector{Triangle{T}}()) #this will extend the array with a new Vector of Triangles.
    end

    temp = pop!(a.unallocated)
    Base.empty!(temp) #resets count in array but doesn't reclaim space
    push!(a.allocated, temp)
    return temp
end

function empty!(a::TrianglePool{T}) where {T<:Real}
    for i in 1:length(a.allocated)
        push!(a.unallocated, pop!(a.allocated))
    end
end

const threadedtrianglepool = [Dict{DataType,TrianglePool}([Float64 => TrianglePool{Float64}()]) for _ in 1:Threads.nthreads()]

function newintrianglepool!(::Type{T} = Float64, tid::Int = Threads.threadid())::Vector{Triangle{T}} where {T<:Real}
    if T ∉ keys(threadedtrianglepool[tid])
        # if the type of the interval pool has changed then we need to refill it with the correct type
        threadedtrianglepool[tid][T] = TrianglePool{T}()
    end
    return allocate!(threadedtrianglepool[tid][T])
end

function emptytrianglepool!(::Type{T} = Float64, tid::Int = Threads.threadid()) where {T}
    if T ∈ keys(threadedtrianglepool[tid])
        Optics.empty!(threadedtrianglepool[tid][T])
    end
end

########################################################################################################################

# currrently we assume that no triangles span 2 sides of a shape, i.e. any triangle only crosses one shape boundary and so splitting is trivial
function uniontri!(csg::CSGTree{T}, tri::Triangle{T}, triangles::Vector{Triangle{T}}, thisforcoplanar = false) where {T<:Real}
    v1 = vertex(tri, 1)
    v2 = vertex(tri, 2)
    v3 = vertex(tri, 3)
    in1 = inside(csg, v1)
    in2 = inside(csg, v2)
    in3 = inside(csg, v3)
    # second condition for if the vertex is in the plane of the other object, in which case we want to keep one copy of the two coplanar surfaces
    if !thisforcoplanar
        # if we are NOT using this surface for coplanar faces then we want to count the verts on the surface as
        # being inside too in case we have coplanar faces
        in1 |= onsurface(csg, v1)
        in2 |= onsurface(csg, v2)
        in3 |= onsurface(csg, v3)
    end
    if !(in1 || in2 || in3) || (thisforcoplanar && onsurface(csg, v1) && onsurface(csg, v2) && onsurface(csg, v3))
        # whole triangle is valid
        push!(triangles, tri)
    elseif !in1 && !in2 && in3
        splittri2out!(csg, v1, v2, v3, triangles)
    elseif !in2 && !in3 && in1
        splittri2out!(csg, v2, v3, v1, triangles)
    elseif !in1 && !in3 && in2
        splittri2out!(csg, v3, v1, v2, triangles)
    elseif !in1 && in2 && in3
        splittri1out!(csg, v1, v2, v3, triangles)
    elseif !in2 && in1 && in3
        splittri1out!(csg, v2, v3, v1, triangles)
    elseif !in3 && in1 && in2
        splittri1out!(csg, v3, v1, v2, triangles)
    end
end

function intersecttri!(csg::CSGTree{T}, tri::Triangle{T}, triangles::Vector{Triangle{T}}, thisforcoplanar = false) where {T<:Real}
    v1 = vertex(tri, 1)
    v2 = vertex(tri, 2)
    v3 = vertex(tri, 3)
    in1 = inside(csg, v1)
    in2 = inside(csg, v2)
    in3 = inside(csg, v3)
    if thisforcoplanar
        # if we are allowing dupes then we want to count the verts on the surface as being inside too
        # in case we have coplanar faces
        in1 |= onsurface(csg, v1)
        in2 |= onsurface(csg, v2)
        in3 |= onsurface(csg, v3)
    end
    if (in1 && in2 && in3)
        # whole triangle is valid
        push!(triangles, tri)
    elseif in1 && in2 && !in3
        splittri2out!(csg, v1, v2, v3, triangles)
    elseif in2 && in3 && !in1
        splittri2out!(csg, v2, v3, v1, triangles)
    elseif in1 && in3 && !in2
        splittri2out!(csg, v3, v1, v2, triangles)
    elseif in1 && !in2 && !in3
        splittri1out!(csg, v1, v2, v3, triangles)
    elseif in2 && !in1 && !in3
        splittri1out!(csg, v2, v3, v1, triangles)
    elseif in3 && !in1 && !in2
        splittri1out!(csg, v3, v1, v2, triangles)
    end
end

function splittri1out!(csg::CSGTree{T}, outv::SVector{3,T}, inv1::SVector{3,T}, inv2::SVector{3,T}, triangles::Vector{Triangle{T}}, recurse::Int = 0) where {T<:Real}
    if !validtri(outv, inv1, inv2) || recurse > VIS_RECURSION_LIMIT
        # if this triangle is invalid then just stop
        return
    end
    n1 = norm(outv - inv1)
    n2 = norm(outv - inv2)
    n3 = norm(inv1 - inv2)
    if min(n1, n2, n3) < MIN_VIS_TRI_SIZE
        # triangle is too small to bother splitting so just add it as is
        push!(triangles, Triangle(outv, inv1, inv2))
        return
    end
    n1 *= 1.0 + MESH_PRECISION
    n2 *= 1.0 + MESH_PRECISION
    int1 = closestintersection(evalcsg(csg, Ray(outv, inv1 - outv)), false)
    int2 = closestintersection(evalcsg(csg, Ray(outv, inv2 - outv)), false)
    if !(int1 === nothing || int2 === nothing) && validtri(outv, point(int1), point(int2)) && α(int1) <= n1 && α(int2) <= n2
        mid = (inv1 + inv2) / 2
        test = closestintersection(evalcsg(csg, Ray(outv, mid - outv)), false)
        if test !== nothing
            if !isapprox(point(test), mid, atol = MESH_PRECISION)
                splittri1out!(csg, outv, point(int1), point(test), triangles, recurse + 1)
                splittri1out!(csg, outv, point(test), point(int2), triangles, recurse + 1)
            else
                push!(triangles, Triangle(outv, point(int1), point(int2)))
            end
        else
            push!(triangles, Triangle(outv, point(int1), point(int2)))
        end
    end
end

function splittri2out!(csg::CSGTree{T}, outv1::SVector{3,T}, outv2::SVector{3,T}, inv::SVector{3,T}, triangles::Vector{Triangle{T}}, flip::Bool = false, recurse::Int = 0) where {T<:Real}
    if !validtri(outv1, outv2, inv) || recurse > VIS_RECURSION_LIMIT
        # if this triangle is invalid then just stop
        return
    end
    n1 = norm(outv1 - inv)
    n2 = norm(outv2 - inv)
    n3 = norm(outv1 - outv2)
    if min(n1, n2, n3) < MIN_VIS_TRI_SIZE
        # triangle is too small to bother splitting so just add it as is
        push!(triangles, Triangle(outv1, outv2, inv))
        return
    end
    n1 *= 1.0 + MESH_PRECISION
    n2 *= 1.0 + MESH_PRECISION
    int1 = closestintersection(evalcsg(csg, Ray(outv1, inv - outv1)), false)
    int2 = closestintersection(evalcsg(csg, Ray(outv2, inv - outv2)), false)
    if int2 !== nothing && α(int2) <= n2
        # int2 is valid
        if int1 !== nothing && α(int1) <= n1
            if isapprox(point(int1), point(int2), atol = MESH_PRECISION)
                # we just need one triangle in this case as inv lies on the boundary of the shape
                if validtri(outv1, outv2, point(int2))
                    if flip
                        push!(triangles, Triangle(outv1, point(int2), outv2))
                    else
                        push!(triangles, Triangle(outv1, outv2, point(int2)))
                    end
                end
            else
                if flip
                    splittri1out!(csg, outv1, point(int1), point(int2), triangles)
                else
                    splittri1out!(csg, outv1, point(int2), point(int1), triangles)
                end
                splittri2out!(csg, outv2, outv1, point(int2), triangles, !flip, recurse + 1)
            end
        else
            if validtri(outv1, outv2, point(int2))
                if flip
                    push!(triangles, Triangle(outv1, point(int2), outv2))
                else
                    push!(triangles, Triangle(outv1, outv2, point(int2)))
                end
            end
        end
    elseif int1 !== nothing && α(int1) <= n1
        if validtri(outv1, outv2, point(int1))
            if flip
                push!(triangles, Triangle(outv1, point(int1), outv2))
            else
                push!(triangles, Triangle(outv1, outv2, point(int1)))
            end
        end
    end
end

# visualization functions from CSG objects, most actualy work happens in the above functions
function makemeshi(a::UnionNode{T}, subdivisions::Int = 30)::Vector{Triangle{T}} where {T<:Real}
    # do the two children in parallel
    th1 = Threads.@spawn makemesh(a.leftchild, subdivisions)
    th2 = Threads.@spawn makemesh(a.rightchild, subdivisions)
    ltris = (fetch(th1)::TriangleMesh{T}).triangles
    rtris = (fetch(th2)::TriangleMesh{T}).triangles
    triangles = newintrianglepool!(T)
    @inbounds for tri in ltris
        uniontri!(a.rightchild, tri, triangles, true)
    end
    @inbounds for tri in rtris
        uniontri!(a.leftchild, tri, triangles)
    end
    return triangles
end

function makemeshi(a::IntersectionNode{T}, subdivisions::Int = 30)::Vector{Triangle{T}} where {T<:Real}
    # do the two children in parallel
    th1 = Threads.@spawn makemesh(a.leftchild, subdivisions)
    th2 = Threads.@spawn makemesh(a.rightchild, subdivisions)
    ltris = (fetch(th1)::TriangleMesh{T}).triangles
    rtris = (fetch(th2)::TriangleMesh{T}).triangles
    triangles = newintrianglepool!(T)
    @inbounds for tri in ltris
        intersecttri!(a.rightchild, tri, triangles, true)
    end
    @inbounds for tri in rtris
        intersecttri!(a.leftchild, tri, triangles)
    end
    return triangles
end

function makemeshi(a::ComplementNode{T}, subdivisions::Int = 30)::Vector{Triangle{T}} where {T<:Real}
    child_triangles = makemeshi(a.child, subdivisions)
    # flip the normals
    triangles = newintrianglepool!(T)
    @inbounds for i in 1:length(child_triangles)
        tri = child_triangles[i]
        push!(triangles, Triangle(vertex(tri, 3), vertex(tri, 2), vertex(tri, 1)))
    end
    return triangles
end

function makemeshi(a::LeafNode{T}, subdivisions::Int = 30)::Vector{Triangle{T}} where {T<:Real}
    child_triangles = triangulate(a.geometry, subdivisions)
    triangles = newintrianglepool!(T)
    @inbounds for i in 1:length(child_triangles)
        push!(triangles, a.transform * child_triangles[i])
    end
    return triangles
end

function makemesh(c::CSGTree{T}, subdivisions::Int = 30)::TriangleMesh{T} where {T<:Real}
    m = TriangleMesh(copy(makemeshi(c, subdivisions))) # need to copy with multithreading
    emptytrianglepool!(T)
    emptyintervalpool!(T)
    return m
end
