# Emitters

!!! warning
    **The old emitter implementation is deprecated as of v0.5!** See below for the new API.

Emitters create rays in a certain pattern, usually controlled by some parameters. Emitters are defined by Pixels and Spatial Layouts, and have a spectrum and an optical power distribution over the hemisphere. These are intrinsic physical properties of the emitter.

The **basic emitter** ([`Source`](@ref sources)) is constructed as a combination of 4 basic elements and a 3D [`Transform`](@ref). The basic elements include:
- [`Spectrum`](@ref spectrum)
- [`Angular Power Distribution`](@ref angular_power_distribution)
- [`Rays Origins Distribution`](@ref rays_origins_distribution)
- [`Rays Directions Distribution`](@ref rays_directions_distribution)

The [`OpticSim`](index.html) package comes with various implementations of each of these basic elements:
- Spectrum - the **generate** interface returns a tuple (power, wavelength)
  * [`Emitters.Spectrum.Uniform`](@ref) - A flat spectrum bounded (default: from 450nm to 680nm). the range is sampled uniformly.
  * [`Emitters.Spectrum.DeltaFunction`](@ref) - Constant wave length.
  * [`Emitters.Spectrum.Measured`](@ref) - measured spectrum to compute emitter power and wavelength (created by reading CSV files – more details will follow).
- Angular Power Distribution - the interface **apply** returns an OpticalRay with modified power
  * [`Emitters.AngularPower.Lambertian`](@ref)
  * [`Emitters.AngularPower.Cosine`](@ref)
  * [`Emitters.AngularPower.Gaussian`](@ref)
- Rays Origins Distribution - the interface **length** returns the number of samples, and **generate** returns the n'th sample.
  * [`Emitters.Origins.Point`](@ref) - a single point
  * [`Emitters.Origins.RectUniform`](@ref) - a uniformly sampled rectangle with user defined number of samples
  * [`Emitters.Origins.RectGrid`](@ref) - a rectangle sampled in a grid fashion
  * [`Emitters.Origins.Hexapolar`](@ref) - a circle (or an ellipse) sampled in an hexapolar fashion (rings)
- Rays Directions Distribution - the interface **length** returns the number of samples, and **generate** returns the n'th sample.
  * [`Emitters.Directions.Constant`](@ref)
  * [`Emitters.Directions.RectGrid`](@ref)
  * [`Emitters.Directions.UniformCone`](@ref)
  * [`Emitters.Directions.HexapolarCone`](@ref)

## [Examples of Basic Emitters](@id basic_emitters)

**Note**: All of the examples on this page assume that the following statement was executed:

```@example 
dummy = "switching to WGLMakie in order to produce interactive figures with Makie"  # hide
using OpticSim, OpticSim.Geometry, OpticSim.Emitters                                # hide
OpticSim.NotebooksUtils.SetDocsBackend("Web")                                       # hide
```

```@example
using OpticSim, OpticSim.Geometry, OpticSim.Emitters
```

```@example
OpticSim.Examples
draw_gestel()
```

### Simple functions for creating commonly used emitters
Many optical systems by convention have their optical axis parallel to the z axis. These utility functions provide a simple interface to the Emitters package to create emitters that direct their rays in the negative z direction, toward the entrance of the optical system.

```@docs
Emitters.pointemitter
```

```@example
using OpticSim, OpticSim.Geometry, OpticSim.Emitters # hide
pt = Emitters.pointemitter([0.0,0.0,.5],.3)
Vis.draw(pt, debug=true, resolution = (800, 600))
```

```@docs
Emitters.collimatedemitter
```

```@example
using OpticSim, OpticSim.Geometry, OpticSim.Emitters # hide
pt = Emitters.collimatedemitter([0.0,0.0,1.0],.5)
Vis.draw(pt, debug=true, resolution = (800, 600))
```

### Point origin with various Direction distributions
```@example
using OpticSim, OpticSim.Geometry, OpticSim.Emitters # hide
src = Sources.Source(origins=Origins.Point(), directions=Directions.RectGrid(π/4, π/4, 15, 15))
Vis.draw(src, debug=true, resolution = (800, 600))
```

```@example
using OpticSim, OpticSim.Geometry, OpticSim.Emitters # hide
src = Sources.Source(origins=Origins.Point(), directions=Directions.UniformCone(π/6, 1000))
Vis.draw(src, debug=true, resolution = (800, 600))
```


```@example
using OpticSim, OpticSim.Geometry, OpticSim.Emitters # hide
src = Sources.Source(origins=Origins.Point(), directions=Directions.HexapolarCone(π/6, 10))
Vis.draw(src, debug=true, resolution = (800, 600))
```


### Various origins distributions

```@example
using OpticSim, OpticSim.Geometry, OpticSim.Emitters # hide
src = Sources.Source(origins=Origins.RectGrid(1.0, 1.0, 10, 10), directions=Directions.Constant())
Vis.draw(src, debug=true, resolution = (800, 600))
```


```@example
using OpticSim, OpticSim.Geometry, OpticSim.Emitters # hide
src = Sources.Source(origins=Origins.Hexapolar(5, 1.0, 2.0), directions=Directions.Constant())
Vis.draw(src, debug=true, resolution = (800, 600))
```


### Examples of Angular Power Distribution

In these example, the arrow width is proportional to the ray power.

```@example
using OpticSim, OpticSim.Geometry, OpticSim.Emitters # hide
src = Sources.Source(
    origins=Origins.Hexapolar(1, 8.0, 8.0),             # Hexapolar Origins
	directions=Directions.RectGrid(π/6, π/6, 15, 15),   # RectGrid Directions
	power=AngularPower.Cosine(10.0)                     # Cosine Angular Power 
)
Vis.draw(src, debug=true, resolution = (800, 600))
```

```@example
using OpticSim, OpticSim.Geometry, OpticSim.Emitters # hide
src = Sources.Source(
	origins=Origins.RectGrid(1.0, 1.0, 3, 3),           # RectGrid Origins  
	directions=Directions.HexapolarCone(π/6, 10),       # HexapolarCone Directions
	power=AngularPower.Gaussian(2.0, 2.0)               # Gaussian Angular Power 
)
Vis.draw(src, debug=true, resolution = (800, 600))
```

### Composite Sources - Display Example

```@example
using OpticSim, OpticSim.Geometry, OpticSim.Emitters # hide

# construct the emitter's basic components
S = Spectrum.Uniform()
P = AngularPower.Lambertian()
O = Origins.RectGrid(1.0, 1.0, 3, 3)
D = Directions.HexapolarCone(deg2rad(5.0), 3)	
	
# construct the source. in this example a "pixel" source will contain only one source as we are simulating a "b/w" display. 
# for RGB displays we can combine 3 sources to simulate "a pixel".
Tr = Transform(Vec3(0.5, 0.5, 0.0))
source1 = Sources.Source(Tr, S, O, D, P)
	
# create a list of pixels - each one is a composite source
pixels = Vector{Sources.CompositeSource{Float64}}(undef, 0)
for y in 1:5 # image_height
    for x in 1:10 # image_width
        # pixel position relative to the display's origin
        local pixel_position = Vec3((x-1) * 1.1, (y-1) * 1.5, 0.0)
        local Tr = Transform(pixel_position)

        # constructing the "pixel"
        pixel = Sources.CompositeSource(Tr, [source1])

        push!(pixels, pixel)
    end
end
	
Tr = Transform(Vec3(0.0, 0.0, 0.0))
my_display = Sources.CompositeSource(Tr, pixels)

Vis.draw(my_display)                                                # render the display - nothing but the origins primitives
rays = AbstractArray{OpticalRay{Float64, 3}}(collect(my_display))   # collect the rays generated by the display
Vis.draw!(rays)                                                     # render the rays 
```


```@example 
dummy = "switching Back to the GLMakie to allow the rest of the pages to work with static figures"  # hide
using OpticSim, OpticSim.Geometry, OpticSim.Emitters                                                # hide
OpticSim.NotebooksUtils.SetDocsBackend("Static")                                                    # hide
```


## [Spectrum](@id spectrum)

```@docs
Emitters.Spectrum.Uniform
Emitters.Spectrum.DeltaFunction
Emitters.Spectrum.Measured
Emitters.Spectrum.spectrumpower
```

## [Angular Power Distribution](@id angular_power_distribution)

```@docs
Emitters.AngularPower.Lambertian
Emitters.AngularPower.Cosine
Emitters.AngularPower.Gaussian
```

## [Rays Origins Distribution](@id rays_origins_distribution)

```@docs
Emitters.Origins.Point
Emitters.Origins.RectUniform
Emitters.Origins.RectGrid
Emitters.Origins.Hexapolar
```

## [Rays Directions Distribution](@id rays_directions_distribution)

```@docs
Emitters.Directions.Constant
Emitters.Directions.RectGrid
Emitters.Directions.UniformCone
Emitters.Directions.HexapolarCone
```

## [Sources](@id sources)

```@docs
Emitters.Sources.Source
Emitters.Sources.CompositeSource
```






