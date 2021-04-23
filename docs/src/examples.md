# Examples

## Pluto Notebooks

The [`OpticSim`](index.html) package comes with several [`Pluto`](https://github.com/fonsp/Pluto.jl) notebooks (code snippets are coming soon) that allow the user to change and run sample code and view the results in real-time. We highly recommend for you to try these out.
The notebooks are located in the **samples** folder, and you can try them by running:

```julia
import OpticSim.NotebooksUtils as NB    # the **as** option was added in Julia v1.6

NB.run_sample("EmittersIntro.jl")
```

The **run_sample** method will **copy** the notebook to your current folder (if it does not exist) and launch Pluto to run the notebook in the browser.

## Cooke Triplet
```@example
using CodeTracking, OpticSim.Examples # hide
print(@code_string cooketripletlensonly()) # hide
```

```@example
using OpticSim.Examples
sys = cooketripletlensonly()

using OpticSim.Geometry, OpticSim.Emitters
origins = Origins.Hexapolar(8, 15.0, 15.0)
directions = Directions.Constant(0.0, 0.0, -1.0)
s1 = Sources.Source(; origins, directions, sourcenum=1)
s2 = Sources.Source(; transform=Transform(rotmatd(10, 0, 0), unitZ3()), origins, directions, sourcenum=2)
raygenerator = Sources.CompositeSource(Transform(), [s1, s2])

using OpticSim.Vis
Vis.drawtracerays(sys; raygenerator, test=true, trackallrays=true, colorbysourcenum=true, resolution=(1000, 700))
Vis.make2dy(); Vis.save("assets/cooke.png") # hide
sys # hide
```

![Cooke triplet visualization](assets/cooke.png)

## Zoom Lens
```@example
using CodeTracking, OpticSim.Examples # hide
print(@code_string zoom_lens(1)) # hide
```

```@example
using OpticSim.Examples
@show zoom_lens(1)

using OpticSim.Geometry, OpticSim.Emitters
transform = translation(0.0, 0.0, 10.0)
origins = Origins.Hexapolar(8, 10.0, 10.0)
directions = Directions.Constant(0.0, 0.0, -1.0)
raygenerator = Sources.Source(; transform, origins, directions)

using OpticSim.Vis
for i = 1:3
    Vis.drawtracerays(zoom_lens(i); raygenerator, test=true, trackallrays=true, numdivisions=50, resolution=(1200, 600))
    Vis.make2dy(); Vis.save("assets/zoom$i.png") # hide
end
nothing # hide
```

![Zoom position 1 visualization](assets/zoom1.png)
![Zoom position 2 visualization](assets/zoom2.png)
![Zoom position 3 visualization](assets/zoom3.png)

## Schmidt Cassegrain Telescope
```@example
using CodeTracking, OpticSim.Examples # hide
print(@code_string SchmidtCassegrainTelescope()) # hide
```

```@example
using OpticSim.Examples
sys = SchmidtCassegrainTelescope()

using OpticSim.Geometry, OpticSim.Emitters
transform = translation(0.0, 0.0, 10.0)
origins = Origins.Hexapolar(8, 20.0, 20.0)
directions = Directions.Constant(0.0, 0.0, -1.0)
raygenerator = Sources.Source(; transform, origins, directions)

using OpticSim.Vis
Vis.drawtracerays(sys; raygenerator, trackallrays = true, colorbynhits = true, test = true, numdivisions = 100)
Vis.save("assets/tele.png") # hide
sys # hide
```

![Schmidt Cassegrain Telescope visualization](assets/tele.png)

## Lens Construction

```@example
using OpticSim, OpticSim.Geometry
using StaticArrays

topsurface = leaf(AcceleratedParametricSurface(QTypeSurface(9.0, radius = -25.0, conic = 0.3, αcoeffs = [(1, 0, 0.3), (1, 1, 1.0)], βcoeffs = [(1, 0, -0.1), (2, 0, 0.4), (3, 0, -0.6)], normradius = 9.5), interface = FresnelInterface{Float64}(OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air)), translation(0.0, 0.0, 5.0))
botsurface = leaf(Plane(0.0, 0.0, -1.0, 0.0, 0.0, -5.0, vishalfsizeu = 9.5, vishalfsizev = 9.5, interface = FresnelInterface{Float64}(OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air)))
barrel = leaf(Cylinder(9.0, 20.0, interface = FresnelInterface{Float64}(OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, reflectance = zero(Float64), transmission = zero(Float64))))
lens = csgintersection(barrel, csgintersection(topsurface, botsurface))(Transform{Float64}(0.0, Float64(π), 0.0, 0.0, 0.0, -5.0))
sys = CSGOpticalSystem(LensAssembly(lens), Rectangle(15.0, 15.0, [0.0, 0.0, 1.0], [0.0, 0.0, -67.8], interface = opaqueinterface()))
Vis.drawtracerays(sys, test = true, trackallrays = true)
Vis.save("assets/qtype.png") # hide
nothing # hide
```

![lens construction example](assets/qtype.png)

## HOEs

### Focusing
```@example
using CodeTracking # hide
using OpticSim.Examples: HOEfocus # hide
print(@code_string HOEfocus()) # hide
```

```@example
using OpticSim.Examples: HOEfocus
HOEfocus()
using OpticSim.Vis; Vis.save("assets/hoe_f.png") # hide
nothing # hide
```

![Focusing HOE example](assets/hoe_f.png)

### Collimating
```@example
using CodeTracking # hide
using OpticSim.Examples: HOEcollimate # hide
print(@code_string HOEcollimate()) # hide
```

```@example
using OpticSim.Examples: HOEcollimate
HOEcollimate()
using OpticSim.Vis; Vis.save("assets/hoe_c.png") # hide
nothing # hide
```

![Collimating HOE example](assets/hoe_c.png)

### Multi
```@example
using CodeTracking # hide
using OpticSim.Examples: multiHOE # hide
print(@code_string multiHOE()) # hide
```

```@example
using OpticSim.Examples: multiHOE
multiHOE()
using OpticSim.Vis; Vis.save("assets/hoe_m.png") # hide
nothing # hide
```

![Multi-HOE example](assets/hoe_m.png)

## Deterministic Raytracing

```@example
using OpticSim
using OpticSim.GlassCat
using OpticSim.Geometry
using StaticArrays

function stacked_beamsplitters(interfacemode)
    bs_1 = leaf(
             leaf(
               Cuboid(10.0, 20.0, 2.0; 
                      interface = FresnelInterface{Float64}(SCHOTT.N_BK7, Air; 
                                                   reflectance=0.5, transmission=0.5, 
                                                   interfacemode=interfacemode)),
               rotationX(pi/4)),
             translation(0.0, 0.0, -30.0-2*sqrt(2))) 
    l1 = leaf(SphericalLens(SCHOTT.N_BK7, -70.0, 30.0, Inf, 5.0, 10.0), translation(0.0, -1.34, 0.0))
    bs_2 = leaf(
             leaf(
               Cuboid(10.0, 20.0, 2.0; 
                      interface = FresnelInterface{Float64}(SCHOTT.N_BK7, Air; 
                                                   reflectance=0.5, transmission=0.5, 
                                                   interfacemode=interfacemode)),
               rotationX(pi/4)),
              translation(0.0, 40.0, -30.0+2*sqrt(2)))
    l2 = leaf(SphericalLens(SCHOTT.N_BK7, -70.0, 30.0, Inf, 5.0, 10.0), translation(0.0, 40.0, 0.0))
    la = LensAssembly(bs_1(), l1(), bs_2(), l2())
    detector = Rectangle(20.0, 40.0, SVector(0.0, 0.0, 1.0), SVector(0.0, 20.0, -130.0); interface = opaqueinterface())
    CSGOpticalSystem(la, detector)
end

# nondeterministic
Vis.drawtracerays(stacked_beamsplitters(ReflectOrTransmit); trackallrays=true, rayfilter=nothing, colorbynhits=true)
Vis.save("assets/deterministic_trace_1.png") # hide
# deterministic, all beamsplitters transmissive
Vis.drawtracerays(stacked_beamsplitters(Transmit); trackallrays=true, rayfilter=nothing, colorbynhits=true)
Vis.save("assets/deterministic_trace_2.png") # hide
# deterministic, all beamsplitters reflective
Vis.drawtracerays(stacked_beamsplitters(Reflect); trackallrays=true, rayfilter=nothing, colorbynhits=true)
Vis.save("assets/deterministic_trace_3.png") # hide
nothing # hide
```

![Nondeterministic Raytrace](assets/deterministic_trace_1.png)
![Transmission only](assets/deterministic_trace_2.png)
![Reflection only](assets/deterministic_trace_3.png)