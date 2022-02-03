# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

const Examples_N_BK7 = GlassCat.Glass(GlassCat.GlassID(GlassCat.AGF, 885), 2, 1.03961212, 0.00600069867, 0.231792344, 0.0200179144, 1.01046945, 103.560653, 0.0, 0.0, NaN, NaN, 0.3, 2.5, 1.86e-6, 1.31e-8, -1.37e-11, 4.34e-7, 6.27e-10, 0.17, 20.0, -0.0009, 2.3, 1.0, 7.1, 1.0, 1, 1.0, [(0.3, 0.05, 25.0), (0.31, 0.25, 25.0), (0.32, 0.52, 25.0), (0.334, 0.78, 25.0), (0.35, 0.92, 25.0), (0.365, 0.971, 25.0), (0.37, 0.977, 25.0), (0.38, 0.983, 25.0), (0.39, 0.989, 25.0), (0.4, 0.992, 25.0), (0.405, 0.993, 25.0), (0.42, 0.993, 25.0), (0.436, 0.992, 25.0), (0.46, 0.993, 25.0), (0.5, 0.994, 25.0), (0.546, 0.996, 25.0), (0.58, 0.995, 25.0), (0.62, 0.994, 25.0), (0.66, 0.994, 25.0), (0.7, 0.996, 25.0), (1.06, 0.997, 25.0), (1.53, 0.98, 25.0), (1.97, 0.84, 25.0), (2.325, 0.56, 25.0), (2.5, 0.36, 25.0)], 1.5168, 2.3, 0.0, 0, 64.17, 0, 2.51, 0.0)


intersectionat(α::Float64) = Intersection(α, [0.7071067811865475, 0.7071067811865475, 0.0] * α, [0.0, 0.0, 0.0], 0.0, 0.0, NullInterface())

const greenwavelength = Unitful.u"550nm"
const redwavelength = Unitful.u"650nm"
const bluewavelength = Unitful.u"450nm"
const roomtemperature = 20 * Unitful.u"°C"

transmissivethingrating(period, orders) = ThinGratingInterface(SVector(0.0, 1.0, 0.0), period, OpticSim.GlassCat.Air, OpticSim.GlassCat.Air, minorder = -orders, maxorder = orders, reflectance = [0.0 for _ in (-orders):orders], transmission = [0.2 for _ in (-orders):orders])
reflectivethingrating(period, orders) = ThinGratingInterface(SVector(0.0, 1.0, 0.0), period, OpticSim.GlassCat.Air, OpticSim.GlassCat.Air, minorder = -orders, maxorder = orders, transmission = [0.0 for _ in (-orders):orders], reflectance = [0.2 for _ in (-orders):orders])

function multiHOE()
    rect = Rectangle(5.0, 5.0, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, 0.0))
    inbeam = SVector(0.0, 0.0, -1.0), CollimatedBeam
    int1 = HologramInterface(SVector(-5.0, 0.0, -20.0), ConvergingBeam, SVector(0.0, -1.0, -1.0), CollimatedBeam, 0.55, 100.0, OpticSim.GlassCat.Air, Examples_N_BK7, OpticSim.GlassCat.Air, OpticSim.GlassCat.Air, OpticSim.GlassCat.Air, 0.05, false)
    int2 = HologramInterface(SVector(5.0, 0.0, -20.0), ConvergingBeam, SVector(0.0, 1.0, -1.0), CollimatedBeam, 0.55, 100.0, OpticSim.GlassCat.Air, Examples_N_BK7, OpticSim.GlassCat.Air, OpticSim.GlassCat.Air, OpticSim.GlassCat.Air, 0.05, false)
    mint = MultiHologramInterface(int1, int2)
    obj = MultiHologramSurface(rect, mint)
    sys = CSGOpticalSystem(LensAssembly(obj), Rectangle(50.0, 50.0, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, -20.0), interface = opaqueinterface()))
    return sys
end
