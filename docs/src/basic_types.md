# Basic Types

The following are basic geometric types and utilities, such as vectors and transforms, that are used all acoross the package. Using these types requires, in addition to the `using OpticSim` statement to also use the OpticSim.Geometry module like:

```@example
using OpticSim, OpticSim.Geometry
```

## Vec3

Representing a 3D vector.

```@docs
OpticSim.Geometry.Vec3
OpticSim.Geometry.unitX3
OpticSim.Geometry.unitY3
OpticSim.Geometry.unitZ3
```

## Vec4

Representing a 4D vector 

```@docs
OpticSim.Geometry.Vec4
OpticSim.Geometry.unitX4
OpticSim.Geometry.unitY4
OpticSim.Geometry.unitZ4
OpticSim.Geometry.unitW4
```

## Transform

Representing a general 3D transform (4x4 matrix). Currently only used as a rigid-body transform.
Transforms are used to position surfaces within the CSG tree, position emitters in 3D, etc. 

```@docs
OpticSim.Geometry.Transform
OpticSim.Geometry.identitytransform
OpticSim.Geometry.rotationX
OpticSim.Geometry.rotationY
OpticSim.Geometry.rotationZ
OpticSim.Geometry.rotation
OpticSim.Geometry.rotationd
OpticSim.Geometry.rotate
OpticSim.Geometry.translation
OpticSim.Geometry.rotmat
OpticSim.Geometry.rotmatd
OpticSim.Geometry.rotmatbetween
```
