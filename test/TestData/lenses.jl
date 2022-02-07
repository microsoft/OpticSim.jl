# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

function zernikelens()
    topsurface = leaf(AcceleratedParametricSurface(ZernikeSurface(9.5, zcoeff = [(1, 0.0), (4, -1.0), (7, 0.5)]), interface = FresnelInterface{Float64}(OpticSim.Examples.Examples_N_BK7, OpticSim.GlassCat.Air)), translation(0.0, 0.0, 5.0))
    botsurface = leaf(Plane(0.0, 0.0, -1.0, 0.0, 0.0, -5.0, vishalfsizeu = 9.5, vishalfsizev = 9.5, interface = FresnelInterface{Float64}(OpticSim.Examples.Examples_N_BK7, OpticSim.GlassCat.Air)))
    barrel = leaf(Cylinder(9.0, 20.0, interface = FresnelInterface{Float64}(OpticSim.Examples.Examples_N_BK7, OpticSim.GlassCat.Air, reflectance = zero(Float64), transmission = zero(Float64))))
    return (barrel ∩ topsurface ∩ botsurface)(Transform(0.0, Float64(π), 0.0, 0.0, 0.0, -5.0))
end

function evenasphericlens()
    topsurface = leaf(AcceleratedParametricSurface(AsphericSurface(9.5, radius = -50.0,aspherics = [(2, 0.0), (4, -0.001)]), interface = FresnelInterface{Float64}(OpticSim.Examples.Examples_N_BK7, OpticSim.GlassCat.Air)), Geometry.translation(0.0, 0.0, 5.0))
    botsurface = leaf(Plane(0.0, 0.0, -1.0, 0.0, 0.0, -5.0, vishalfsizeu = 9.5, vishalfsizev = 9.5, interface = FresnelInterface{Float64}(OpticSim.Examples.Examples_N_BK7, OpticSim.GlassCat.Air)))
    barrel = leaf(Cylinder(9.0, 20.0, interface = FresnelInterface{Float64}(OpticSim.Examples.Examples_N_BK7, OpticSim.GlassCat.Air, reflectance = zero(Float64), transmission = zero(Float64))))
    return (barrel ∩ topsurface ∩ botsurface)(Transform(0.0, Float64(π), 0.0, 0.0, 0.0, -5.0))
end

function oddasphericlens()
    topsurface = leaf(AcceleratedParametricSurface(AsphericSurface(9.5, radius = -50.0,aspherics = [(3, 0.001), (5, -0.0001)]), interface = FresnelInterface{Float64}(OpticSim.Examples.Examples_N_BK7, OpticSim.GlassCat.Air)), Geometry.translation(0.0, 0.0, 5.0))
    botsurface = leaf(Plane(0.0, 0.0, -1.0, 0.0, 0.0, -5.0, vishalfsizeu = 9.5, vishalfsizev = 9.5, interface = FresnelInterface{Float64}(OpticSim.Examples.Examples_N_BK7, OpticSim.GlassCat.Air)))
    barrel = leaf(Cylinder(9.0, 20.0, interface = FresnelInterface{Float64}(OpticSim.Examples.Examples_N_BK7, OpticSim.GlassCat.Air, reflectance = zero(Float64), transmission = zero(Float64))))
    return (barrel ∩ topsurface ∩ botsurface)(Transform(0.0, Float64(π), 0.0, 0.0, 0.0, -5.0))
end

function asphericlens()
    topsurface = leaf(AcceleratedParametricSurface(AsphericSurface(9.5, aspherics = [(3, 0.0), (4, -0.001)]), interface = FresnelInterface{Float64}(OpticSim.Examples.Examples_N_BK7, OpticSim.GlassCat.Air)), translation(0.0, 0.0, 5.0))
    botsurface = leaf(Plane(0.0, 0.0, -1.0, 0.0, 0.0, -5.0, vishalfsizeu = 9.5, vishalfsizev = 9.5, interface = FresnelInterface{Float64}(OpticSim.Examples.Examples_N_BK7, OpticSim.GlassCat.Air)))
    barrel = leaf(Cylinder(9.0, 20.0, interface = FresnelInterface{Float64}(OpticSim.Examples.Examples_N_BK7, OpticSim.GlassCat.Air, reflectance = zero(Float64), transmission = zero(Float64))))
    return (barrel ∩ topsurface ∩ botsurface)(Transform(0.0, Float64(π), 0.0, 0.0, 0.0, -5.0))
end

function coniclensA()
    topsurface = leaf(AcceleratedParametricSurface(AsphericSurface(9.5, radius = -15.0, conic = -5.0), interface = FresnelInterface{Float64}(OpticSim.Examples.Examples_N_BK7, OpticSim.GlassCat.Air)), translation(0.0, 0.0, 5.0))
    botsurface = leaf(Plane(0.0, 0.0, -1.0, 0.0, 0.0, -5.0, vishalfsizeu = 9.5, vishalfsizev = 9.5, interface = FresnelInterface{Float64}(OpticSim.Examples.Examples_N_BK7, OpticSim.GlassCat.Air)))
    barrel = leaf(Cylinder(9.0, 20.0, interface = FresnelInterface{Float64}(OpticSim.Examples.Examples_N_BK7, OpticSim.GlassCat.Air, reflectance = zero(Float64), transmission = zero(Float64))))
    return (barrel ∩ topsurface ∩ botsurface)(Transform(0.0, Float64(π), 0.0, 0.0, 0.0, -5.0))
end

function coniclensQ()
    topsurface = leaf(AcceleratedParametricSurface(QTypeSurface(9.5, radius = -15.0, conic = -5.0), interface = FresnelInterface{Float64}(OpticSim.Examples.Examples_N_BK7, OpticSim.GlassCat.Air)), translation(0.0, 0.0, 5.0))
    botsurface = leaf(Plane(0.0, 0.0, -1.0, 0.0, 0.0, -5.0, vishalfsizeu = 9.5, vishalfsizev = 9.5, interface = FresnelInterface{Float64}(OpticSim.Examples.Examples_N_BK7, OpticSim.GlassCat.Air)))
    barrel = leaf(Cylinder(9.0, 20.0, interface = FresnelInterface{Float64}(OpticSim.Examples.Examples_N_BK7, OpticSim.GlassCat.Air, reflectance = zero(Float64), transmission = zero(Float64))))
    return (barrel ∩ topsurface ∩ botsurface)(Transform(0.0, Float64(π), 0.0, 0.0, 0.0, -5.0))
end

function qtypelens()
    topsurface = leaf(AcceleratedParametricSurface(QTypeSurface(9.0, radius = -25.0, conic = 0.3, αcoeffs = [(1, 0, 0.3), (1, 1, 1.0)], βcoeffs = [(1, 0, -0.1), (2, 0, 0.4), (3, 0, -0.6)], normradius = 9.5), interface = FresnelInterface{Float64}(OpticSim.Examples.Examples_N_BK7, OpticSim.GlassCat.Air)), translation(0.0, 0.0, 5.0))
    botsurface = leaf(Plane(0.0, 0.0, -1.0, 0.0, 0.0, -5.0, vishalfsizeu = 9.5, vishalfsizev = 9.5, interface = FresnelInterface{Float64}(OpticSim.Examples.Examples_N_BK7, OpticSim.GlassCat.Air)))
    barrel = leaf(Cylinder(9.0, 20.0, interface = FresnelInterface{Float64}(OpticSim.Examples.Examples_N_BK7, OpticSim.GlassCat.Air, reflectance = zero(Float64), transmission = zero(Float64))))
    return (barrel ∩ topsurface ∩ botsurface)(Transform(0.0, Float64(π), 0.0, 0.0, 0.0, -5.0))
end

function chebyshevlens()
    topsurface = leaf(AcceleratedParametricSurface(ChebyshevSurface(9.0, 9.0, [(1, 0, -0.081), (2, 0, -0.162), (0, 1, 0.243), (1, 1, 0.486), (2, 1, 0.081), (0, 2, -0.324), (1, 2, -0.405), (2, 2, -0.81)], radius = -25.0), interface = FresnelInterface{Float64}(OpticSim.Examples.Examples_N_BK7, OpticSim.GlassCat.Air)), translation(0.0, 0.0, 5.0))
    botsurface = leaf(Plane(0.0, 0.0, -1.0, 0.0, 0.0, -5.0, vishalfsizeu = 9.5, vishalfsizev = 9.5, interface = FresnelInterface{Float64}(OpticSim.Examples.Examples_N_BK7, OpticSim.GlassCat.Air)))
    barrel = leaf(Cylinder(9.0, 20.0, interface = FresnelInterface{Float64}(OpticSim.Examples.Examples_N_BK7, OpticSim.GlassCat.Air, reflectance = zero(Float64), transmission = zero(Float64))))
    return (barrel ∩ topsurface ∩ botsurface)(Transform(0.0, Float64(π), 0.0, 0.0, 0.0, -5.0))
end

function gridsaglens()
    topsurface = leaf(AcceleratedParametricSurface(GridSagSurface{Float64}(joinpath(@__DIR__, "../../test/test.GRD"), radius = -25.0, conic = 0.4, interpolation = GridSagBicubic), interface = FresnelInterface{Float64}(OpticSim.Examples.Examples_N_BK7, OpticSim.GlassCat.Air)), translation(0.0, 0.0, 5.0))
    botsurface = leaf(Plane(0.0, 0.0, -1.0, 0.0, 0.0, -5.0, vishalfsizeu = 9.5, vishalfsizev = 9.5, interface = FresnelInterface{Float64}(OpticSim.Examples.Examples_N_BK7, OpticSim.GlassCat.Air)))
    barrel = leaf(Cylinder(9.0, 20.0, interface = FresnelInterface{Float64}(OpticSim.Examples.Examples_N_BK7, OpticSim.GlassCat.Air, reflectance = zero(Float64), transmission = zero(Float64))))
    return (barrel ∩ topsurface ∩ botsurface)(Transform(0.0, Float64(π), 0.0, 0.0, 0.0, -5.0))
end
