# Primitives

All geometry is built up from a small(ish) number of primitives and a number of constructive solid geometry (CSG) operations (see [CSG](@ref)).
Primitives are split into two types, `Surface`s and `ParametricSurface`s, the latter being a subset of the former.
`Surface`s are standalone surfaces which cannot be used in CSG operations, e.g. an aperture or rectangle.
`ParametricSurface`s are valid csg objects and can be composed into very complex structures.

## Surfaces

A surface can be any surface in 3D space, it can be bounded and not create a half-space (i.e. not partition space into _inside_ and _outside_).

```@docs
Surface
```

### Basic Shapes

These are the simple shapes with are provided already, they act only as standalone objects and cannot be used in CSG objects.
Adding a new `Surface` is easy, the new structure must simply follow the interface defined above.

```@docs
Ellipse
Circle
Rectangle
Hexagon
Triangle
TriangleMesh
```

### Stops

A number of simple occlusive apertures are provided as constructing such objects using CSG can be inefficient and error-prone.

```@docs
InfiniteStop
FiniteStop
RectangularAperture
CircularAperture
Annulus
```

## Parametric Surfaces

A parametric surface must partition space into two valid half-spaces, i.e. _inside_ and _outside_.
The surface must also be parameterized by two variables, nominally `u` and `v`.
Typically these surfaces cannot be intersected with a ray analytically and so must be triangulated and an iterative solution found.

```@docs
ParametricSurface
AcceleratedParametricSurface
```

### Parametric Surface Types

These are the available types of parametric surfaces which are already implemented, all of which can be used in the creation of CSG objects.
New `ParametricSurface`s can be added with relative ease providing they follow the interface defined above.

```@docs
Cylinder
Plane
Sphere
SphericalCap
ZernikeSurface
BezierSurface
BSplineSurface
QTypeSurface
ChebyshevSurface
GridSagSurface
```

## Functions

These are some useful functions related to `Surface` objects.

```@docs
point(::ParametricSurface{T}, ::T, ::T) where {T<:Real}
normal
partials(::ParametricSurface{T}, ::T, ::T) where {T<:Real}
uvrange
uv
Opticks.uvtopix
inside(s::ParametricSurface{T,3}, x::T, y::T, z::T) where {T<:Real}
onsurface(s::ParametricSurface{T,3}, x::T, y::T, z::T) where {T<:Real}
interface
surfaceintersection(surf::AcceleratedParametricSurface{S,N}, r::AbstractRay{S,N}) where {S,N}
samplesurface
triangulate
makemesh
```

## Bounding Boxes

Bounding boxes are mostly used internally for efficiency, but are also exposed to the user for visualization (and any other) purposes.
All bounding boxes are axis aligned.

```@docs
BoundingBox
doesintersect
surfaceintersection(::BoundingBox{T}, ::AbstractRay{T,3}) where {T<:Real}
```
