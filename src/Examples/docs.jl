# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

# Group examples that are used in the docs (examples.md)
export cooketripletlensonly, SchmidtCassegrainTelescope, zoom_lens, HOEfocus, HOEcollimate, multiHOE

cooketripletlensonly(::Type{T} = Float64) where {T<:Real} = AxisymmetricOpticalSystem{T}(DataFrame(
    Surface      = [:Object, 1,             2,      3,            :Stop,  5,             6,       :Image ],
    Radius       = [Inf,     26.777,        66.604, -35.571,      35.571, 35.571,        -26.777, Inf    ],
    Thickness    = [Inf,     4.0,           2.0,    4.0,          2.0,    4.0,           44.748,  missing],
    Material     = [Air,     SCHOTT.N_SK16, Air,    SCHOTT.N_SF2, Air,    SCHOTT.N_SK16, Air,     missing],
    SemiDiameter = [Inf,     8.580,         7.513,  7.054,        6.033,  7.003,         7.506,   15.0   ]))

function SchmidtCassegrainTelescope()
    # glass entrance lens on telescope
    topsurf = Plane(
        SVector(0.0, 0.0, 1.0),
        SVector(0.0, 0.0, 0.0),
        interface = FresnelInterface{Float64}(SCHOTT.N_BK7, Air),
        vishalfsizeu = 12.00075,
        vishalfsizev = 12.00075)
    botsurf = AcceleratedParametricSurface(ZernikeSurface(
        12.00075,
        radius = -1.14659768e+4,
        aspherics = [(4, 3.68090959e-7), (6, 2.73643352e-11), (8, 3.20036892e-14)]),
        17,
        interface = FresnelInterface{Float64}(SCHOTT.N_BK7, Air))
    coverlens = csgintersection(
        leaf(Cylinder(12.00075, 1.4)),
        csgintersection(
            leaf(topsurf), 
            leaf(botsurf, Transform(rotmatd(0, 180, 0), Vec3(0.0, 0.0, -0.65)))))

    # big mirror with a hole in it
    bigmirror = csgdifference(
        ConicLens(SCHOTT.N_BK7, -72.65, -95.2773500000134, 0.077235, Inf, 0.0, 0.2, 12.18263, frontsurfacereflectance = 1.0),
        leaf(Cylinder(4.0, 0.3, interface = opaqueinterface()), translation(0.0, 0.0, -72.75)))

    # small mirror supported on a spider
    smallmirror = SphericalLens(SCHOTT.N_BK7, -40.65, Inf, -49.6845, 1.13365, 4.3223859, backsurfacereflectance = 1.0)

    obscuration1 = Circle(4.5, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, -40.649), interface = opaqueinterface())
    obscurations2 = Spider(3, 0.5, 12.0, SVector(0.0, 0.0, -40.65))

    # put it together with the detector
    la = LensAssembly(coverlens(), bigmirror(), smallmirror(), obscuration1, obscurations2...)
    det = Circle(3.0, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, -92.4542988), interface = opaqueinterface())

    return CSGOpticalSystem(la, det)
end

function zoom_lens(pos)
    stops = [2.89, 3.99, 4.90]
    zooms = [9.48, 4.48, 2.00]
    dists = [4.46970613, 21.21, 43.81]
    return AxisymmetricOpticalSystem{Float64}(DataFrame(
        Surface = [:Object, :Stop, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, :Image],
        Radius = [Inf64, Inf64, -1.6202203499676E+01, -4.8875855327468E+01, 1.5666614444619E+01, -4.2955326460481E+01, 1.0869565217391E+02, 2.3623907394283E+01, -1.6059097478722E+01, -4.2553191489362E+02, -3.5435861091425E+01, -1.4146272457208E+01, -2.5125628140704E+02, -2.2502250225023E+01, -1.0583130489999E+01, -4.4444444444444E+01, Inf64],
        Aspherics = [missing, missing, missing, missing, missing, [(4, 1.0386E-04), (6, 1.4209E-07), (8, -8.8495E-09), (10, 1.2477E-10), (12, -1.0367E-12), (14, 3.6556E-15)], missing, missing, [(4, 4.2721E-05), (6, 1.2484E-07), (8, 9.7079E-09), (10, -1.8444E-10), (12, 1.8644E-12), (14, -7.7975E-15)], [(4, 1.1339E-04), (6, 4.8165E-07), (8, 1.8778E-08), (10, -5.7571E-10), (12, 8.9994E-12), (14, -4.6768E-14)], missing, missing, missing, missing, missing, missing, missing],
        Thickness = [Inf64, 0.0, 5.18, 0.10, 4.40, 0.16, 1.0, 4.96, zooms[pos], 4.04, 1.35, 1.0, 2.80, 3.0, 1.22, dists[pos], missing],
        Material = [Air, Air, OHARA.S_LAH66, Air, NIKON.LLF6, Air, OHARA.S_TIH6, OHARA.S_FSL5, Air, OHARA.S_FSL5, Air, OHARA.S_LAL8, OHARA.S_FSL5, Air, OHARA.S_LAH66, Air, missing], # S_FL4 -> S_FSL5>
        SemiDiameter = [Inf64, stops[pos], 3.85433218451, 3.85433218451, 4.36304692871, 4.36304692871, 4.72505505439, 4.72505505439, 4.72505505439, 4.45240784026, 4.45240784026, 4.50974054117, 4.50974054117, 4.50974054117, 4.76271114409, 4.76271114409, 15.0]))
end

function HOEfocus(; kwargs...)
    rect = Rectangle(5.0, 5.0, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, 0.0))
    int = HologramInterface(
        SVector(0.0, -3.0, -20.0), ConvergingBeam,
        SVector(0.0, 0.0, -1.0), CollimatedBeam,
        0.55, 9.0, Air, SCHOTT.N_BK7, Air, Air, Air, 0.05, false)
    obj = HologramSurface(rect, int)
    sys = CSGOpticalSystem(
        LensAssembly(obj),
        Rectangle(10.0, 10.0, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, -25.0),
        interface = opaqueinterface()))
    raygenerator = Source(;
        transform = translation(0.0, 0.0, 10.0),
        spectrum = DeltaFunction(0.55),
        origins = Origins.RectGrid(3.0, 3.0, 5, 5),
        directions = Constant(0.0, 0.0, -1.0))
    Vis.drawtracerays(sys; raygenerator, trackallrays = true, rayfilter = nothing, test = true, kwargs...)
end

function HOEcollimate(; kwargs...)
    rect = Rectangle(5.0, 5.0, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, 0.0))
    int = HologramInterface(
        SVector(0.1, -0.05, -1.0), CollimatedBeam,
        SVector(0.0, 0.0, 10), DivergingBeam,
        0.55, 9.0, Air, SCHOTT.N_BK7, Air, Air, Air, 0.05, false)
    obj = HologramSurface(rect, int)
    sys = CSGOpticalSystem(
        LensAssembly(obj),
        Rectangle(10.0, 10.0, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, -25.0),
        interface = opaqueinterface()))
    raygenerator = Source(
        transform = Transform(rotmatd(180, 0, 0), Vec3(0.0, 0.0, 10.0)),
        spectrum = DeltaFunction(0.55),
        origins = Origins.Point(),
        directions = Directions.RectGrid(π/4, π/4, 8, 8))
    Vis.drawtracerays(sys; raygenerator, trackallrays = true, rayfilter = nothing, test = true, kwargs...)
end

function multiHOE(; kwargs...)
    rect = Rectangle(5.0, 5.0, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, 0.0))
    int1 = HologramInterface(
        SVector(-5.0, 0.0, -20.0), ConvergingBeam,
        SVector(0.0, -1.0, -1.0), CollimatedBeam,
        0.55, 100.0, Air, SCHOTT.N_BK7, Air, Air, Air, 0.05, false)
    int2 = HologramInterface(
        SVector(5.0, 0.0, -20.0), ConvergingBeam,
        SVector(0.0, 1.0, -1.0), CollimatedBeam,
        0.55, 100.0, Air, SCHOTT.N_BK7, Air, Air, Air, 0.05, false)
    mint = MultiHologramInterface(int1, int2)
    obj = MultiHologramSurface(rect, mint)
    sys = CSGOpticalSystem(
        LensAssembly(obj),
        Rectangle(10.0, 10.0, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, -20.0), interface = opaqueinterface()))
    spectrum = DeltaFunction(0.55)
    origins = Origins.RectUniform(3.0, 3.0, 500)
    directions = Constant(0.0, 0.0, -1.0)
    s1 = Source(; spectrum, origins, directions, sourcenum = 1,
        transform = Transform(rotmatd(-45, 0, 0), Vec3(0.0, 3.0, 3.0)))
    s2 = Source(; spectrum, origins, directions, sourcenum = 2,
        transform = Transform(rotmatd(45, 0, 0), Vec3(0.0, -3.0, 3.0)))
    s3 = Source(; spectrum, origins, directions, sourcenum = 3,
        transform = translation(0.0, 0.0, 3.0))
    Vis.drawtracerays(sys; raygenerator = s1, trackallrays = true, colorbysourcenum = true, rayfilter = nothing, kwargs...)
    Vis.drawtracerays!(sys; raygenerator = s2, trackallrays = true, colorbysourcenum = true, rayfilter = nothing, drawgen = true, kwargs...)
    Vis.drawtracerays!(sys; raygenerator = s3, trackallrays = true, colorbysourcenum = true, rayfilter = nothing, drawgen = true, kwargs...)
end
