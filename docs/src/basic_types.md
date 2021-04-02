# Basic Types

The following are basic geometric types and utilities, such as vectors and transforms, that are used all acoross the package. Using these types requires, in addition to the `using OpticSim` statement to also use the OpticSim.Geoemtry module like:

```@example
using OpticSim, OpticSim.Geometry
```

## Vec3

Representing a 3D vector.

```@docs
Geometry.Vec3
Geometry.unitX3
Geometry.unitY3
Geometry.unitZ3
Geometry.zero3
Geometry.one3
```

## Vec4

Representing a 4D vector 

```@docs
Geometry.Vec4
Geometry.unitX4
Geometry.unitY4
Geometry.unitZ4
Geometry.unitW4
Geometry.zero4
Geometry.one4
```

## Transform

Representing a general 3D transform (4x4 matrix). Currently only used as a rigid-body transform.
Transforms are used to position surfaces within the CSG tree, position emitters in 3D, etc. 

```@docs
Geometry.Transform
Geometry.identitytransform
Geometry.rotationX
Geometry.rotationY
Geometry.rotationZ
Geometry.rotation
Geometry.rotationd
Geometry.rotate
Geometry.translation
Geometry.rotmat
Geometry.rotmatd
Geometry.rotmatbetween
```

