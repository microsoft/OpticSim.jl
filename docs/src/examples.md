# Examples

## Cooke Triplet

```@example
using OpticSim

using DataFrames
sys = AxisymmetricOpticalSystem(
    DataFrame(Surface = [:Object, 1, 2, 3, :Stop, 5, 6, :Image],
              Radius = [Inf, 26.777, 66.604, -35.571, 35.571, 35.571, -26.777, Inf],
              Thickness = [Inf, 4.0, 2.0, 4.0, 2.0, 4.0, 44.748, missing],
              Material = [OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_SK16, OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_SF2, OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_SK16, OpticSim.GlassCat.Air, missing],
              SemiDiameter = [Inf, 8.580, 7.513, 7.054, 6.033, 7.003, 7.506, 15.0]))
@show sys
f1 = HexapolarField(sys, collimated = true, samples = 4, sourcenum = 1)
f2 = HexapolarField(sys, collimated = true, samples = 4, sourceangle = -10 / 180 * π, sourcenum = 2)
Vis.drawtracerays(sys, raygenerator = f1, test = true, trackallrays = true, colorbysourcenum = true, resolution = (1000, 700))
Vis.drawtracerays!(sys, raygenerator = f2, test = true, trackallrays = true, colorbysourcenum = true)
Vis.make2dy() # hide
Vis.save("assets/cooke.png") # hide
nothing # hide
```

![Cooke triplet visualization](assets/cooke.png)

## Schmidt Cassegrain Telescope

```@example
using OpticSim, OpticSim.Geometry
using StaticArrays


# glass entrance lens on telescope
topsurf = Plane(SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, 0.0), interface = FresnelInterface{Float64}(OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air), vishalfsizeu = 12.00075, vishalfsizev = 12.00075)
botsurf = AcceleratedParametricSurface(ZernikeSurface(12.00075, radius = -1.14659768e+4, aspherics = [(4, 3.68090959e-7), (6, 2.73643352e-11), (8, 3.20036892e-14)]), 17, interface = FresnelInterface{Float64}(OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air))
coverlens = csgintersection(leaf(Cylinder(12.00075, 1.4)), csgintersection(leaf(topsurf), leaf(botsurf, Transform(OpticSim.rotmatd(0, 180, 0), Vec3(0.0, 0.0, -0.65)))))
# big mirror with a hole in it
bigmirror = ConicLens(OpticSim.GlassCat.SCHOTT.N_BK7, -72.65, -95.2773500000134, 0.077235, Inf, 0.0, 0.2, 12.18263, frontsurfacereflectance = 1.0)
bigmirror = csgdifference(bigmirror, leaf(Cylinder(4.0, 0.3, interface = opaqueinterface()), translation(0.0, 0.0, -72.75)))
# small mirror supported on a spider
smallmirror = SphericalLens(OpticSim.GlassCat.SCHOTT.N_BK7, -40.65, Inf, -49.6845, 1.13365, 4.3223859, backsurfacereflectance = 1.0)
obscuration1 = Circle(4.5, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, -40.649), interface = opaqueinterface())
obscurations2 = Spider(3, 0.5, 12.0, SVector(0.0, 0.0, -40.65))
# put it together with the detector
la = LensAssembly(coverlens(), bigmirror(), smallmirror(), obscuration1, obscurations2...)
det = Circle(3.0, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, -92.4542988), interface = opaqueinterface())
tele = CSGOpticalSystem(la, det)

Vis.drawtracerays(tele, raygenerator = UniformOpticalSource(CollimatedSource(GridRectOriginPoints(5, 5, 10.0, 10.0, position = SVector(0.0, 0.0, 20.0))), 0.55), trackallrays = true, colorbynhits = true, test = true)
Vis.save("assets/tele.png") # hide
nothing # hide
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
using OpticSim
using StaticArrays

rect = Rectangle(5.0, 5.0, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, 0.0))
int = HologramInterface(SVector(0.0, -3.0, -20.0), ConvergingBeam, SVector(0.0, 0.0, -1.0), CollimatedBeam, 0.55, 9.0, OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, OpticSim.GlassCat.Air, OpticSim.GlassCat.Air, 0.05, false)
obj = HologramSurface(rect, int)
sys = CSGOpticalSystem(LensAssembly(obj), Rectangle(10.0, 10.0, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, -25.0), interface = opaqueinterface()))
Vis.drawtracerays(sys; raygenerator = UniformOpticalSource(CollimatedSource(GridRectOriginPoints(5, 5, 3.0, 3.0, position = SVector(0.0, 0.0, 10.0), direction = SVector(0.0, 0.0, -1.0))), 0.55), trackallrays = true, rayfilter = nothing, test = true)
Vis.save("assets/hoe_f.png") # hide
nothing # hide
```

![Focusing HOE example](assets/hoe_f.png)

### Collimating

```@example
using OpticSim
using StaticArrays

rect = Rectangle(5.0, 5.0, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, 0.0))
int = HologramInterface(SVector(0.1, -0.05, -1.0), CollimatedBeam, SVector(0.0, 0.0, 10), DivergingBeam, 0.55, 9.0, OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, OpticSim.GlassCat.Air, OpticSim.GlassCat.Air, 0.05, false)
obj = HologramSurface(rect, int)
sys = CSGOpticalSystem(LensAssembly(obj), Rectangle(10.0, 10.0, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, -25.0), interface = opaqueinterface()))
Vis.drawtracerays(sys; raygenerator = UniformOpticalSource(GridSource(OriginPoint{Float64}(1, position = SVector(0.0, 0.0, 10.0), direction = SVector(0.0, 0.0, -1.0)), 5, 5, π / 4, π / 4), 0.55), trackallrays = true, rayfilter = nothing, test = true)
Vis.save("assets/hoe_c.png") # hide
nothing # hide
```

![Collimating HOE example](assets/hoe_c.png)

### Multi

```@example
using OpticSim
using StaticArrays

rect = Rectangle(5.0, 5.0, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, 0.0))
int1 = HologramInterface(SVector(-5.0, 0.0, -20.0), ConvergingBeam, SVector(0.0, -1.0, -1.0), CollimatedBeam, 0.55, 100.0, OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, OpticSim.GlassCat.Air, OpticSim.GlassCat.Air, 0.05, false)
int2 = HologramInterface(SVector(5.0, 0.0, -20.0), ConvergingBeam, SVector(0.0, 1.0, -1.0), CollimatedBeam, 0.55, 100.0, OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, OpticSim.GlassCat.Air, OpticSim.GlassCat.Air, 0.05, false)
mint = MultiHologramInterface(int1, int2)
obj = MultiHologramSurface(rect, mint)
sys = CSGOpticalSystem(LensAssembly(obj), Rectangle(10.0, 10.0, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, -20.0), interface = opaqueinterface()))
s1 = UniformOpticalSource(CollimatedSource(RandomRectOriginPoints(500, 3.0, 3.0, position = SVector(0.0, 3.0, 3.0), direction = SVector(0.0, -1.0, -1.0))), 0.55, sourcenum = 1)
s2 = UniformOpticalSource(CollimatedSource(RandomRectOriginPoints(500, 3.0, 3.0, position = SVector(0.0, -3.0, 3.0), direction = SVector(0.0, 1.0, -1.0))), 0.55, sourcenum = 2)
s3 = UniformOpticalSource(CollimatedSource(RandomRectOriginPoints(500, 3.0, 3.0, position = SVector(0.0, 0.0, 3.0), direction = SVector(0.0, 0.0, -1.0))), 0.55, sourcenum = 3)
Vis.drawtracerays(sys; raygenerator = s1, trackallrays = true, colorbysourcenum = true, rayfilter = nothing)
Vis.drawtracerays!(sys; raygenerator = s2, trackallrays = true, colorbysourcenum = true, rayfilter = nothing, drawgen = true)
Vis.drawtracerays!(sys; raygenerator = s3, trackallrays = true, colorbysourcenum = true, rayfilter = nothing, drawgen = true)
Vis.save("assets/hoe_m.png") # hide
nothing # hide
```

![Multi-HOE example](assets/hoe_m.png)
