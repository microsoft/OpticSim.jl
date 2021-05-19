# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

using OpticSim
using OpticSim.GlassCat
using OpticSim.Geometry

using StaticArrays
using AbstractTrees

top = CSGTree2(
    AcceleratedParametricSurface(
        QTypeSurface(
            9.0,
            radius = -25.0,
            conic = 0.3,
            αcoeffs = [(1, 0, 0.3), (1, 1, 1.0)],
            βcoeffs = [(1, 0, -0.1), (2, 0, 0.4), (3, 0, -0.6)],
            normradius = 9.5),
        interface = FresnelInterface{Float64}(SCHOTT.N_BK7, Air)
    ),
    translation(0.0, 0.0, 5.0)
)
bot = CSGTree2(
    Plane(
        SVector(0.0, 0.0, -1.0),
        SVector(0.0, 0.0, -5.0),
        vishalfsizeu = 9.5,
        vishalfsizev = 9.5,
        interface = FresnelInterface{Float64}(SCHOTT.N_BK7, Air)
    )
)
barrel = CSGTree2(
    Cylinder(9.0, 20.0, interface = FresnelInterface{Float64}(SCHOTT.N_BK7, Air, reflectance=0.0, transmission=0.0))
)
lens = barrel ∩ top ∩ bot
OpticSim.transform!(lens, Transform{Float64}(0.0, Float64(π), 0.0, 0.0, 0.0, -5.0))
print_tree(lens)
