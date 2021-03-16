# Emitters

Emitters create rays in a certain pattern, usually controlled by some parameters.
The are constructed in a modular way, e.g.

```@example
using OpticSim, StaticArrays # hide
Vis.draw(UniformOpticalSource(CollimatedSource(HexapolarOriginPoints(4, 1.0, 1.0, position = SVector(0.0, 0.0, 10.0), direction = SVector(0.0, sind(30), -cosd(30)))), 0.55))
Vis.save("assets/source1.png") # hide
nothing #hide
```

![Emitter example 1 image](assets/source1.png)

```@example
using OpticSim, StaticArrays # hide
Vis.draw(UniformOpticalSource(GridSource(OriginPoint{Float64}(1, position = SVector(0.0, 0.0, 10.0), direction = SVector(0.0, 0.0, -1.0)), 5, 5, π / 4, π / 4), 0.65))
Vis.save("assets/source2.png") # hide
nothing #hide
```

![Emitter example 2 image](assets/source2.png)

```@example
using OpticSim, StaticArrays # hide
Vis.draw(CosineOpticalSource(RandomSource(OriginPoint{Float64}(200, direction = SVector(0.0, 0.0, 1.0))), 1.0, 0.45))
Vis.save("assets/source3.png") # hide
nothing #hide
```

![Emitter example 3 image](assets/source3.png)

## Rays

```@docs
Ray
OpticalRay
OpticSim.generateray
```

## Origin Points

```@docs
OpticSim.RayOriginGenerator
RandomRectOriginPoints
GridRectOriginPoints
RandomEllipseOriginPoints
HexapolarOriginPoints
OriginPoint
OpticSim.genorigin
```

## Directions

```@docs
GeometricRayGenerator
CollimatedSource
GridSource
RandomSource
OpticSim.gendirection
```

## Power

```@docs
OpticalRayGenerator
UniformOpticalSource
CosineOpticalSource
GaussianOpticalSource
```

## Compound

```@docs
PixelSource
OpticalSourceArray
BasicDisplayPanel
OpticalSourceGroup
RayListSource
```

## Fields

```@docs
HexapolarField
GridField
```
