# CSG

## CSG Operations

There are three binary csg operations which can construct extremely complex objects from very simple primitives: union (``\cup``), intersection (``\cap``) and subtraction (i.e. difference).

This diagram shows the basic idea:
![CSG Tree visualization](https://upload.wikimedia.org/wikipedia/commons/8/8b/Csg_tree.png)

The code for this in our system would look this this:

```@example
using OpticSim # hide
cyl = Cylinder(0.7)
cyl_cross = cyl ∪ leaf(cyl, Geometry.rotationd(90, 0, 0)) ∪ leaf(cyl, Geometry.rotationd(0, 90, 0))

cube = Cuboid(1.0, 1.0, 1.0)
sph = Sphere(1.3)
rounded_cube = cube ∩ sph

result = rounded_cube - cyl_cross
Vis.draw(result, numdivisions=100)

Vis.save("assets/csg_ex.png"); nothing # hide
```

![CSG code example image](assets/csg_ex.png)

```@docs
leaf
∪
∩
-
```

## Pre-made CSG Shapes

There are also some shortcut methods available to create common CSG objects more easily:

```@docs
BoundedCylinder
Cuboid
HexagonalPrism
RectangularPrism
TriangularPrism
Spider
```

## CSG Types

These are the types of the primary CSG elements, i.e. the nodes in the CSG tree.

```@docs
OpticSim.CSGTree
OpticSim.CSGGenerator
OpticSim.ComplementNode
OpticSim.UnionNode
OpticSim.IntersectionNode
OpticSim.LeafNode
```

## Additional Functions and Types

These are the internal types and functions used for geometric/CSG operations.

### Functions

```@docs
surfaceintersection(::CSGTree{T}, ::AbstractRay{T,N}) where {T<:Real,N}
inside(a::CSGTree{T}, x::T, y::T, z::T) where {T<:Real}
onsurface(a::CSGTree{T}, x::T, y::T, z::T) where {T<:Real}
```

### Intervals

```@docs
Interval
EmptyInterval
DisjointUnion
OpticSim.isemptyinterval
OpticSim.ispositivehalfspace
OpticSim.israyorigininterval
OpticSim.halfspaceintersection
OpticSim.closestintersection
OpticSim.IntervalPool
```

### Intersections

```@docs
OpticSim.IntervalPoint
RayOrigin
Infinity
Intersection
OpticSim.isinfinity
OpticSim.israyorigin
```
