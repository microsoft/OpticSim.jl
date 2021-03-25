# Visualization

There are a number of powerful visualization tools available, we primarily rely on 3D visualization of systems using [Makie](http://makie.juliaplots.org/stable/).

There are a number of helper methods, as well as the ability to draw objects, surfaces, points, rays and more individually. For example, looking at rays passing through a system in 3D and 2D:

```julia
Vis.drawtracerays(Examples.cooketriplet(), trackallrays=true, test=true, numdivisions=100)
```

```@setup base
using OpticSim
```

```@example base
Vis.drawtracerays(Examples.cooketriplet(), trackallrays=true, test=true, numdivisions=100)
Vis.save("assets/vis_ex_3d.png") # hide
Vis.drawtracerays(Examples.cooketriplet(), trackallrays=true, test=true, numdivisions=100, drawsys=true, resolution = (1000, 700))
Vis.make2dy()
Vis.save("assets/vis_ex_2d.png"); nothing #hide
```

![3D visualization example](assets/vis_ex_3d.png)
![2D visualization example](assets/vis_ex_2d.png)

And the image on the detector for a trace of a system:

```julia
Vis.drawtraceimage(Examples.cooketriplet(), test=true)
```

```@example base
using Images # hide
im = Vis.drawtraceimage(Examples.cooketriplet(Float64, 400), test=true)
save("assets/vis_ex_im.png", colorview(Gray, real.(im ./ maximum(im)))); nothing # hide
```

![detector image example](assets/vis_ex_im.png)

## Basic Drawing

These methods are all you need to build up a visualization piece by piece.
For example:

```@example base
obj = csgintersection(Sphere(0.5), Plane(0.0, 1.0, 0.0, 0.0, 0.1, 0.0))()
ray1 = Ray([0.0, -0.1, 1.0], [0.0, 0.0, -1.0])
ray2 = Ray([0.8, 0.0, 0.0], [-1.0, 0.0, 0.0])
Vis.draw(obj)
Vis.draw!(ray1, rayscale=0.2)
Vis.draw!(ray2, rayscale=0.2, color=:blue)
Vis.draw!(surfaceintersection(obj, ray1), color=:red)
Vis.draw!(surfaceintersection(obj, ray2), color=:green)
Vis.save("assets/vis_ex_3d_parts.png"); nothing # hide
```

![basic drawing example](assets/vis_ex_3d_parts.png)

```@docs
OpticSim.Vis.scene
OpticSim.Vis.draw
OpticSim.Vis.draw!(::Any; kwargs...)
OpticSim.Vis.save
```

## Helper Methods

These are the helper methods to provide common visualizations more easily, as used above. Another example:

```julia
Vis.surfacesag(AcceleratedParametricSurface(TestData.zernikesurface2()), (256, 256), (1.55, 1.55))
```

```@example base
using Plots; include("../../test/TestData/TestData.jl") # hide
p = Vis.surfacesag(AcceleratedParametricSurface(TestData.zernikesurface2()), (256, 256), (1.55, 1.55))
Plots.savefig(p, "assets/surface_sag.svg"); nothing # hide
```

![surface sag example](assets/surface_sag.svg)

```@docs
OpticSim.Vis.drawtracerays
OpticSim.Vis.drawtraceimage
OpticSim.Vis.spotdiag
OpticSim.Vis.surfacesag
OpticSim.Vis.eyebox_eval_eye
OpticSim.Vis.eyebox_eval_planar
```

## Complete Drawing Functions

As mentioned above, `Vis.draw!` can be used to draw a large variety of objects, each with their own additional arguments.
Here is a full list of the available drawing function and their associated options.

```@docs
OpticSim.Vis.draw!
```

## Known Issues

If the Makie plot is printing to the console rather than showing properly in a separate window then call `Vis.AbstractPlotting.__init__()`. This will only occur when using a system image including the OpticSim package.