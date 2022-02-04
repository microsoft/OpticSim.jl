# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

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
