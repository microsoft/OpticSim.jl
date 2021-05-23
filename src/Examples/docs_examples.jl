# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

# Group examples that are used in the docs (examples.md)
export draw_cooketriplet, draw_schmidtcassegraintelescope, draw_lensconstruction, draw_zoomlenses, draw_HOEfocus, draw_HOEcollimate, draw_multiHOE, draw_stackedbeamsplitters

function draw_cooketriplet(filename::Union{Nothing,AbstractString} = nothing)
    g1, g2 = SCHOTT.N_SK16, SCHOTT.N_SF2
    sys = AxisymmetricOpticalSystem{Float64}(DataFrame(
        SurfaceType  = ["Object", "Standard", "Standard", "Standard", "Stop", "Standard", "Standard", "Image"],
        Radius       = [Inf,      26.777,     66.604,     -35.571,    35.571, 35.571,     -26.777,    Inf    ],
        Thickness    = [Inf,      4.0,        2.0,        4.0,        2.0,    4.0,        44.748,     missing],
        Material     = [Air,      g1,         Air,        g2,         Air,    g1,         Air,        missing],
        SemiDiameter = [Inf,      8.580,      7.513,      7.054,      6.033,  7.003,      7.506,      15.0   ],
    ))

    origins = Origins.Hexapolar(8, 15.0, 15.0)
    directions = Directions.Constant(0.0, 0.0, -1.0)
    s1 = Sources.Source(; origins, directions, sourcenum=1)

    transform = Transform(rotmatd(10, 0, 0), unitZ3())
    s2 = Sources.Source(; transform, origins, directions, sourcenum=2)

    raygenerator = Sources.CompositeSource(Transform(), [s1, s2])

    trackallrays = test = colorbysourcenum = true; resolution = (1000, 700)
    Vis.drawtracerays(sys; raygenerator, trackallrays, test, colorbysourcenum, resolution)
    Vis.make2dy(); Vis.save(filename)
    return sys
end

function draw_zoomlenses(filenames::Vector{<:Union{Nothing,AbstractString}} = repeat([nothing], 3))
    stops = [2.89, 3.99, 4.90]
    zooms = [9.48, 4.48, 2.00]
    dists = [4.46970613, 21.21, 43.81]

    transform = translation(0.0, 0.0, 10.0)
    origins = Origins.Hexapolar(8, 10.0, 10.0)
    directions = Directions.Constant(0.0, 0.0, -1.0)
    raygenerator = Sources.Source(; transform, origins, directions)

    aspherics = [
        ["4" => 1.0386E-04, "6" => 1.4209E-07, "8" => -8.8495E-09, "10" => 1.2477E-10, "12" => -1.0367E-12, "14" => 3.6556E-15],
        ["4" => 4.2721E-05, "6" => 1.2484E-07, "8" => 9.7079E-09, "10" => -1.8444E-10, "12" => 1.8644E-12, "14" => -7.7975E-15],
        ["4" => 1.1339E-04, "6" => 4.8165E-07, "8" => 1.8778E-08, "10" => -5.7571E-10, "12" => 8.9994E-12, "14" => -4.6768E-14],
    ]
    syss = [
        AxisymmetricOpticalSystem{Float64}(DataFrame(
            SurfaceType = ["Object", "Stop", "Standard", "Standard", "Standard", "Aspheric", "Standard", "Standard", "Aspheric", "Aspheric", "Standard", "Standard", "Standard", "Standard", "Standard", "Standard", "Image"],
            Radius = [Inf64, Inf64, -1.6202203499676E+01, -4.8875855327468E+01, 1.5666614444619E+01, -4.2955326460481E+01, 1.0869565217391E+02, 2.3623907394283E+01, -1.6059097478722E+01, -4.2553191489362E+02, -3.5435861091425E+01, -1.4146272457208E+01, -2.5125628140704E+02, -2.2502250225023E+01, -1.0583130489999E+01, -4.4444444444444E+01, Inf64],
            Parameters = [missing, missing, missing, missing, missing, aspherics[1], missing, missing, aspherics[2], aspherics[3], missing, missing, missing, missing, missing, missing, missing],
            Thickness = [Inf64, 0.0, 5.18, 0.10, 4.40, 0.16, 1.0, 4.96, zoom, 4.04, 1.35, 1.0, 2.80, 3.0, 1.22, dist, missing],
            Material = [Air, Air, OHARA.S_LAH66, Air, NIKON.LLF6, Air, OHARA.S_TIH6, OHARA.S_FSL5, Air, OHARA.S_FSL5, Air, OHARA.S_LAL8, OHARA.S_FSL5, Air, OHARA.S_LAH66, Air, missing],
            SemiDiameter = [Inf64, stop, 3.85433218451, 3.85433218451, 4.36304692871, 4.36304692871, 4.72505505439, 4.72505505439, 4.72505505439, 4.45240784026, 4.45240784026, 4.50974054117, 4.50974054117, 4.50974054117, 4.76271114409, 4.76271114409, 15.0]))
        for (stop, zoom, dist) in zip(stops, zooms, dists)]

    for (sys, filename) in zip(syss, filenames)
        Vis.drawtracerays(sys; raygenerator, trackallrays=true, test=true, numdivisions=50, resolution=(1200, 600))
        Vis.make2dy(); Vis.save(filename)
    end
    return syss
end

function draw_schmidtcassegraintelescope(filename::Union{Nothing,AbstractString} = nothing)
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
    coverlens = Cylinder(12.00075, 1.4) ∩ topsurf ∩ leaf(botsurf, Transform(rotmatd(0, 180, 0), Vec3(0.0, 0.0, -0.65)))

    # big mirror with a hole in it
    frontsurfacereflectance = 1.0
    bigmirror = (
        ConicLens(SCHOTT.N_BK7, -72.65, -95.2773500000134, 0.077235, Inf, 0.0, 0.2, 12.18263; frontsurfacereflectance) -
        leaf(Cylinder(4.0, 0.3, interface = opaqueinterface()), translation(0.0, 0.0, -72.75))
    )

    # small mirror supported on a spider
    backsurfacereflectance = 1.0
    smallmirror = SphericalLens(SCHOTT.N_BK7, -40.65, Inf, -49.6845, 1.13365, 4.3223859; backsurfacereflectance)

    obscuration1 = Circle(4.5, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, -40.649), interface = opaqueinterface())
    obscurations2 = Spider(3, 0.5, 12.0, SVector(0.0, 0.0, -40.65))

    # put it together with the detector
    la = LensAssembly(coverlens(), bigmirror(), smallmirror(), obscuration1, obscurations2...)
    det = Circle(3.0, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, -92.4542988), interface = opaqueinterface())
    sys = CSGOpticalSystem(la, det)

    # define ray generator
    transform = translation(0.0, 0.0, 10.0)
    origins = Origins.Hexapolar(8, 20.0, 20.0)
    directions = Directions.Constant(0.0, 0.0, -1.0)
    raygenerator = Sources.Source(; transform, origins, directions)

    # draw and output
    Vis.drawtracerays(sys; raygenerator, trackallrays = true, colorbynhits = true, test = true, numdivisions = 100, drawgen = false)
    Vis.save(filename)
    return nothing
end

function draw_lensconstruction(filename::Union{Nothing,AbstractString} = nothing)
    topsurface = leaf(
        AcceleratedParametricSurface(
            QTypeSurface(
                9.0,
                radius = -25.0,
                conic = 0.3,
                αcoeffs = [(1, 0, 0.3), (1, 1, 1.0)],
                βcoeffs = [(1, 0, -0.1), (2, 0, 0.4), (3, 0, -0.6)],
                normradius = 9.5),
            interface = FresnelInterface{Float64}(SCHOTT.N_BK7, Air)),
        translation(0.0, 0.0, 5.0))
    botsurface = Plane(
        SVector(0.0, 0.0, -1.0),
        SVector(0.0, 0.0, -5.0),
        vishalfsizeu = 9.5,
        vishalfsizev = 9.5,
        interface = FresnelInterface{Float64}(SCHOTT.N_BK7, Air))
    barrel = Cylinder(
        9.0, 20.0, interface = FresnelInterface{Float64}(SCHOTT.N_BK7, Air, reflectance=0.0, transmission=0.0)
    )
    lens = (barrel ∩ topsurface ∩ botsurface)(Transform{Float64}(0.0, Float64(π), 0.0, 0.0, 0.0, -5.0))
    detector = Rectangle(15.0, 15.0, [0.0, 0.0, 1.0], [0.0, 0.0, -67.8], interface = opaqueinterface())
    sys = CSGOpticalSystem(LensAssembly(lens), detector)

    Vis.drawtracerays(sys, test = true, trackallrays = true, colorbynhits = true)
    Vis.save(filename)
    return nothing
end

function draw_HOEfocus(filename::Union{Nothing,AbstractString} = nothing)
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

    raygenerator = Sources.Source(;
        transform = translation(0.0, 0.0, 10.0),
        spectrum = Spectrum.DeltaFunction(0.55),
        origins = Origins.RectGrid(3.0, 3.0, 5, 5),
        directions = Directions.Constant(0.0, 0.0, -1.0))

    Vis.drawtracerays(sys; raygenerator, trackallrays = true, rayfilter = nothing, test = true)
    Vis.save(filename)
    return nothing
end

function draw_HOEcollimate(filename::Union{Nothing,AbstractString} = nothing)
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

    raygenerator = Sources.Source(
        transform = Transform(rotmatd(180, 0, 0), Vec3(0.0, 0.0, 10.0)),
        spectrum = Spectrum.DeltaFunction(0.55),
        origins = Origins.Point(),
        directions = Directions.RectGrid(π/4, π/4, 8, 8))

    Vis.drawtracerays(sys; raygenerator, trackallrays = true, rayfilter = nothing, test = true)
    Vis.save(filename)
    return nothing
end

function draw_multiHOE(filename::Union{Nothing,AbstractString} = nothing)
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

    spectrum = Spectrum.DeltaFunction(0.55)
    origins = Origins.RectUniform(3.0, 3.0, 500)
    directions = Directions.Constant(0.0, 0.0, -1.0)
    s1 = Sources.Source(; spectrum, origins, directions, sourcenum = 1,
        transform = Transform(rotmatd(-45, 0, 0), Vec3(0.0, 3.0, 3.0)))
    s2 = Sources.Source(; spectrum, origins, directions, sourcenum = 2,
        transform = Transform(rotmatd(45, 0, 0), Vec3(0.0, -3.0, 3.0)))
    s3 = Sources.Source(; spectrum, origins, directions, sourcenum = 3,
        transform = translation(0.0, 0.0, 3.0))
    raygenerator = Sources.CompositeSource(Transform(), [s1, s2, s3])

    Vis.drawtracerays(sys; raygenerator, trackallrays = true, colorbysourcenum = true, rayfilter = nothing, drawgen = true)
    Vis.save(filename)
    return nothing
end

function draw_stackedbeamsplitters(filenames::Vector{<:Union{Nothing,AbstractString}} = repeat([nothing], 3))
    # ReflectOrTransmit: nondeterministic
    # Transmit: deterministic, all beamsplitters transmissive
    # Reflect: deterministic, all beamsplitters reflective
    interfacemodes = [ReflectOrTransmit, Transmit, Reflect]

    for (interfacemode, filename) in zip(interfacemodes, filenames)
        interface = FresnelInterface{Float64}(SCHOTT.N_BK7, Air; reflectance=0.5, transmission=0.5, interfacemode)
        bs_1 = leaf(
            leaf(
                Cuboid(10.0, 20.0, 2.0; interface),
                rotationX(π/4)),
            translation(0.0, 0.0, -30.0-2*sqrt(2)))
        l1 = leaf(
            SphericalLens(SCHOTT.N_BK7, -70.0, 30.0, Inf, 5.0, 10.0),
            translation(0.0, -1.34, 0.0))
        bs_2 = leaf(
            leaf(
                Cuboid(10.0, 20.0, 2.0; interface),
                rotationX(π/4)),
            translation(0.0, 40.0, -30.0+2*sqrt(2)))
        l2 = leaf(
            SphericalLens(SCHOTT.N_BK7, -70.0, 30.0, Inf, 5.0, 10.0),
            translation(0.0, 40.0, 0.0))
        la = LensAssembly(bs_1(), l1(), bs_2(), l2())
        detector = Rectangle(20.0, 40.0, SVector(0.0, 0.0, 1.0), SVector(0.0, 20.0, -130.0); interface = opaqueinterface())
        sys = CSGOpticalSystem(la, detector)

        Vis.drawtracerays(sys; trackallrays=true, rayfilter=nothing, colorbynhits=true)
        Vis.save(filename)
    end
    return nothing
end
