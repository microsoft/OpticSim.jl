# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

export cooketriplet, doubleconvexlensonly

const Examples_N_BK7 = GlassCat.Glass(GlassCat.GlassID(GlassCat.AGF, 885), 2, 1.03961212, 0.00600069867, 0.231792344, 0.0200179144, 1.01046945, 103.560653, 0.0, 0.0, NaN, NaN, 0.3, 2.5, 1.86e-6, 1.31e-8, -1.37e-11, 4.34e-7, 6.27e-10, 0.17, 20.0, -0.0009, 2.3, 1.0, 7.1, 1.0, 1, 1.0, [(0.3, 0.05, 25.0), (0.31, 0.25, 25.0), (0.32, 0.52, 25.0), (0.334, 0.78, 25.0), (0.35, 0.92, 25.0), (0.365, 0.971, 25.0), (0.37, 0.977, 25.0), (0.38, 0.983, 25.0), (0.39, 0.989, 25.0), (0.4, 0.992, 25.0), (0.405, 0.993, 25.0), (0.42, 0.993, 25.0), (0.436, 0.992, 25.0), (0.46, 0.993, 25.0), (0.5, 0.994, 25.0), (0.546, 0.996, 25.0), (0.58, 0.995, 25.0), (0.62, 0.994, 25.0), (0.66, 0.994, 25.0), (0.7, 0.996, 25.0), (1.06, 0.997, 25.0), (1.53, 0.98, 25.0), (1.97, 0.84, 25.0), (2.325, 0.56, 25.0), (2.5, 0.36, 25.0)], 1.5168, 2.3, 0.0, 0, 64.17, 0, 2.51, 0.0)

const Examples_N_SK16 =GlassCat.Glass(GlassCat.GlassID(GlassCat.AGF, 828), 2, 1.34317774, 0.00704687339, 0.241144399, 0.0229005, 0.994317969, 92.7508526, 0.0, 0.0, NaN, NaN, 0.31, 2.5, -2.37e-8, 1.32e-8, -1.29e-11, 4.09e-7, 5.17e-10, 0.17, 20.0, -0.0011, 3.2, 1.4, 6.3, 4.0, 1, 53.3, [(0.31, 0.02, 25.0), (0.32, 0.11, 25.0), (0.334, 0.4, 25.0), (0.35, 0.7, 25.0), (0.365, 0.86, 25.0), (0.37, 0.89, 25.0), (0.38, 0.93, 25.0), (0.39, 0.956, 25.0), (0.4, 0.97, 25.0), (0.405, 0.974, 25.0), (0.42, 0.979, 25.0), (0.436, 0.981, 25.0), (0.46, 0.984, 25.0), (0.5, 0.991, 25.0), (0.546, 0.994, 25.0), (0.58, 0.994, 25.0), (0.62, 0.993, 25.0), (0.66, 0.994, 25.0), (0.7, 0.996, 25.0), (1.06, 0.995, 25.0), (1.53, 0.973, 25.0), (1.97, 0.88, 25.0), (2.325, 0.54, 25.0), (2.5, 0.26, 25.0)], 1.62041, 3.3, 4.0, 0, 60.32, 0, 3.58, 0.0)


const Examples_N_SF2 = GlassCat.Glass(GlassCat.GlassID(GlassCat.AGF, 815), 2, 1.47343127, 0.0109019098, 0.163681849, 0.0585683687, 1.36920899, 127.404933, 0.0, 0.0, NaN, NaN, 0.365, 2.5, 3.1e-6, 1.75e-8, 6.62e-11, 7.51e-7, 8.99e-10, 0.277, 20.0, 0.0081, 1.0, 1.4, 6.68, 1.0, 1, 1.0, [(0.365, 0.007, 25.0), (0.37, 0.06, 25.0), (0.38, 0.4, 25.0), (0.39, 0.68, 25.0), (0.4, 0.83, 25.0), (0.405, 0.865, 25.0), (0.42, 0.926, 25.0), (0.436, 0.949, 25.0), (0.46, 0.961, 25.0), (0.5, 0.975, 25.0), (0.546, 0.986, 25.0), (0.58, 0.987, 25.0), (0.62, 0.984, 25.0), (0.66, 0.984, 25.0), (0.7, 0.987, 25.0), (1.06, 0.997, 25.0), (1.53, 0.984, 25.0), (1.97, 0.93, 25.0), (2.325, 0.76, 25.0), (2.5, 0.67, 25.0)], 1.64769, 1.2, 0.0, 0, 33.82, 0, 2.718, 0.0)

const Examples_N_SF14 = GlassCat.Glass(GlassCat.GlassID(GlassCat.AGF, 896), 2, 1.69022361, 0.0130512113, 0.288870052, 0.061369188, 1.7045187, 149.517689, 0.0, 0.0, NaN, NaN, 0.365, 2.5, -5.56e-6, 7.09e-9, -1.09e-11, 9.85e-7, 1.39e-9, 0.287, 20.0, 0.013, 1.0, 2.2, 9.41, 1.0, 1, 1.0, [(0.365, 0.004, 25.0), (0.37, 0.04, 25.0), (0.38, 0.33, 25.0), (0.39, 0.61, 25.0), (0.4, 0.75, 25.0), (0.405, 0.79, 25.0), (0.42, 0.87, 25.0), (0.436, 0.91, 25.0), (0.46, 0.938, 25.0), (0.5, 0.964, 25.0), (0.546, 0.983, 25.0), (0.58, 0.987, 25.0), (0.62, 0.987, 25.0), (0.66, 0.987, 25.0), (0.7, 0.985, 25.0), (1.06, 0.998, 25.0), (1.53, 0.98, 25.0), (1.97, 0.88, 25.0), (2.325, 0.64, 25.0), (2.5, 0.57, 25.0)], 1.76182, 1.0, 0.0, 0, 26.53, 0, 3.118, 0.0)
"""
    hemisphere()

Create a geometric hemisphere
"""
function hemisphere()::CSGTree
    sph = Sphere(10.0)
    pln = Plane(0.0, 0.0, -1.0, 0.0, 0.0, 0.0)
    # CSV operations create a csggenerator which instantiates the csg tree after applying a rigid body transformation.
    # This allows you to make as many instances of the object as you want with different transformations. We just want
    # the CSGTree object rather than a generator.
    return (sph ∩ pln)()
end

"""
    opticalhemisphere()

Create an optical hemisphere that has optical material properties so it will reflect and refract light. In the previous
example the hemisphere object had optical properties of Air, which is the default optical interface, so it won't refract
or reflect light.
"""
function opticalhemisphere()::CSGOpticalSystem
    sph = Sphere(10.0, interface = FresnelInterface{Float64}(Examples_N_BK7, Air))
    pln = Plane(0.0, 0.0, -1.0, 0.0, 0.0, 0.0, interface = FresnelInterface{Float64}(Examples_N_BK7, Air))
    assy = LensAssembly{Float64}((sph ∩ pln)())
    return CSGOpticalSystem(assy, Rectangle(1.0, 1.0, SVector{3,Float64}(0.0, 0.0, 1.0), SVector{3,Float64}(0.0, 0.0, -11.0)))
end

function cooketriplet(::Type{T} = Float64, detpix::Int = 1000) where {T<:Real}
    AxisymmetricOpticalSystem{T}(
        DataFrame(
            SurfaceType = ["Object", "Standard", "Standard", "Standard", "Stop", "Standard", "Standard", "Image"],
            Radius = [Inf, 26.777, 66.604, -35.571, 35.571, 35.571, -26.777, Inf],
            # OptimizeRadius = [false, true, true, true, true, true, true, false],
            Thickness = [Inf, 4.0, 2.0, 4.0, 2.0, 4.0, 44.748, missing],
            # OptimizeThickness = [false, true, true, true, true, true, true, false],
            Material = [Air, Examples_N_SK16, Air, Examples_N_SF2, Air, Examples_N_SK16, Air, missing],
            SemiDiameter = [Inf, 8.580, 7.513, 7.054, 6.033, 7.003, 7.506, 15.0]
        ),
        detpix,
        detpix
    )
end

function cooketripletfirstelement(::Type{T} = Float64) where {T<:Real}
    AxisymmetricOpticalSystem(
        DataFrame(
            SurfaceType = ["Object", "Standard", "Standard", "Image"],
            Radius = [Inf, -35.571, 35.571, Inf],
            Thickness = [Inf, 4.0, 44.748, missing],
            Material = [Air, Examples_N_SK16, Air, missing],
            SemiDiameter = [Inf, 7.054, 6.033, 15.0]
        )
    )
end

function convexplano(::Type{T} = Float64) where {T<:Real}
    AxisymmetricOpticalSystem{T}(
        DataFrame(
            SurfaceType = ["Object", "Standard", "Standard", "Image"],
            Radius = [Inf, 60.0, Inf, Inf],
            Thickness = [Inf, 10.0, 57.8, missing],
            Material = [Air, Examples_N_BK7, Air, missing],
            SemiDiameter = [Inf, 9.0, 9.0, 15.0]
        )
    )
end



function doubleconvex(frontradius::T, rearradius::T) where {T<:Real}
    AxisymmetricOpticalSystem{T}(
        DataFrame(
            SurfaceType = ["Object", "Standard", "Standard", "Image"],
            Radius = [convert(T, Inf64), frontradius, rearradius, convert(T, Inf64)],
            # OptimizeRadius = [false, true, true, false],
            Thickness = [convert(T, Inf64), convert(T, 10.0), convert(T, 57.8), missing],
            # OptimizeThickness = [false, false, false, false],
            Material = [Air,Examples_N_BK7, Air, missing],
            SemiDiameter = [convert(T, Inf64), convert(T, 9.0), convert(T, 9.0), convert(T, 15.0)]
        )
    )
end

function doubleconvexconic(::Type{T} = Float64) where {T<:Real}
    AxisymmetricOpticalSystem{T}(
        DataFrame(
            SurfaceType = ["Object", "Standard", "Standard", "Image"],
            Radius = [Inf64, 60, -60, Inf64],
            # OptimizeRadius = [false, true, true, false],
            Thickness = [Inf64, 10.0, 57.8, missing],
            # OptimizeThickness = [false, false, false, false],
            Conic = [missing, 0.01, 0.01, missing],
            # OptimizeConic = [false, true, true, false],
            Material = [Air, Examples_N_BK7, Air, missing],
            SemiDiameter = [Inf64, 9.0, 9.0, 15.0]
        )
    )
end

function doubleconvexlensonly(frontradius::T,rearradius::T) where{T<:Real}
    AxisymmetricLens{T}(
        DataFrame(
            SurfaceType = ["Object", "Standard", "Standard", "Image"],
            Radius = [convert(T, Inf64), frontradius, rearradius, convert(T, Inf64)],
            # OptimizeRadius = [false, true, true, false],
            Thickness = [convert(T, Inf64), convert(T, 10.0), convert(T, 57.8), missing],
            # OptimizeThickness = [false, false, false, false],
            Material = [Air, Examples_N_BK7, Air, missing],
            SemiDiameter = [convert(T, Inf64), convert(T, 9.0), convert(T, 9.0), convert(T, 15.0)]
        )
    )
end

function doubleconvex(
    ::Type{T} = Float64;
    temperature::Unitful.Temperature = GlassCat.TEMP_REF_UNITFUL,
    pressure::T = convert(T, PRESSURE_REF)
) where {T<:Real}
    AxisymmetricOpticalSystem{T}(
        DataFrame(
            SurfaceType = ["Object", "Standard", "Standard", "Image"],
            Radius = [Inf64, 60, -60, Inf64],
            # OptimizeRadius = [false, true, true, false],
            Thickness = [Inf64, 10.0, 57.8, missing],
            # OptimizeThickness = [false, true, true, false],
            Material = [Air, Examples_N_BK7, Air, missing],
            SemiDiameter = [Inf64, 9.0, 9.0, 15.0]
        );
        temperature,
        pressure
    )
end

function doubleconcave(::Type{T} = Float64) where {T<:Real}
    AxisymmetricOpticalSystem{T}(
        DataFrame(
            SurfaceType = ["Object", "Standard", "Standard", "Image"],
            Radius = [Inf64, -41.0, 41.0, Inf64],
            Thickness = [Inf64, 10.0, 57.8, missing],
            Material = [Air, Examples_N_BK7, Air, missing],
            SemiDiameter = [Inf64, 9.0, 9.0, 15.0]
        )
    )
end

function planoconcaverefl(::Type{T} = Float64) where {T<:Real}
    AxisymmetricOpticalSystem{T}(
        DataFrame(
            SurfaceType = ["Object", "Standard", "Standard", "Image"],
            Radius = [Inf64, Inf64, -41.0, Inf64],
            Thickness = [Inf64, 10.0, -57.8, missing],
            Material = [Air, Examples_N_BK7, Air, missing],
            SemiDiameter = [Inf64, 9.0, 9.0, 25.0],
            Reflectance = [missing, missing, 1.0, missing]
        )
    )
end

function concaveplano(::Type{T} = Float64) where {T<:Real}
    AxisymmetricOpticalSystem{T}(
        DataFrame(
            SurfaceType = ["Object", "Standard", "Standard", "Image"],
            Radius = [Inf64, -41.0, Inf64, Inf64],
            Thickness = [Inf64, 10.0, 57.8, missing],
            Material = [Air, Examples_N_BK7, Air, missing],
            SemiDiameter = [Inf64, 9.0, 9.0, 15.0]
        )
    )
end

function planoplano(::Type{T} = Float64) where {T<:Real}
    AxisymmetricOpticalSystem{T}(
        DataFrame(
            SurfaceType = ["Object", 1, 2, "Image"],
            Radius = [Inf64, Inf64, Inf64, Inf64],
            Thickness = [Inf64, 10.0, 57.8, missing],
            Material = [Air, Examples_N_BK7, Air, missing],
            SemiDiameter = [Inf64, 9.0, 9.0, 15.0]
        )
    )
end

function autodrawrays(lens::AxisymmetricOpticalSystem = cooketriplet(), angle = 10; kwargs...)
    f1 = HexapolarField(lens, collimated = true, wavelength = 0.45, sourcenum = 1)
    Vis.drawtracerays(lens, raygenerator = f1, test = true, trackallrays = true, colorbysourcenum = true; kwargs...)
    f2 = HexapolarField(lens, collimated = true, wavelength = 0.45, sourceangle = angle / 180 * π, sourcenum = 2)
    Vis.drawtracerays!(lens, raygenerator = f2, test = true, trackallrays = true, colorbysourcenum = true; kwargs...)
end

function autospotdiag(lens::AxisymmetricOpticalSystem = cooketriplet(); kwargs...)
    f1 = HexapolarField(lens, collimated = true, wavelength = 0.45, sourcenum = 1)
    f2 = HexapolarField(lens, collimated = true, wavelength = 0.45, sourceangle = 5 / 180 * π, sourcenum = 2)
    f3 = HexapolarField(lens, collimated = true, wavelength = 0.45, sourceangle = 10 / 180 * π, sourcenum = 3)
    Vis.spotdiaggrid(lens, [f1, f2, f3]; kwargs...)
end

# Display the spot diagram of a simple cooketriplet lens
function hexapolarspotdiagramexample(lens = cooketriplet(), numrings::Int = 5, angle = 0.0)
    Vis.spotdiag(lens, samples = numrings, sourceangle = angle)
end

function cartesiangridspotdiagramexample(lens = cooketriplet(), numsamples::Int = 5, angle = 0.0)
    Vis.spotdiag(lens, hexapolar = false, samples = numsamples, sourceangle = angle)
end

function prism_refraction()
    # build the triangular prism
    int = FresnelInterface{Float64}(Examples_N_SF14, Air)
    s = 2.0
    prism = (
        Plane(
            SVector(0.0, -1.0, 0.0),
            SVector(0.0, -s, 0.0),
            interface = int,
            vishalfsizeu = 2 * s,
            vishalfsizev = 2 * s
        ) ∩
        Plane(
            SVector(0.0, sind(30), cosd(30)),
            SVector(0.0, s * sind(30), s * cosd(30)),
            interface = int,
            vishalfsizeu = 2 * s,
            vishalfsizev = 2 * s
        ) ∩
        Plane(
            SVector(0.0, sind(30), -cosd(30)),
            SVector(0.0, s * sind(30), -s * cosd(30)),
            interface = int,
            vishalfsizeu = 2 * s,
            vishalfsizev = 2 * s
        )
    )
    sys = CSGOpticalSystem(LensAssembly(prism()), Rectangle(15.0, 15.0, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, -20.0), interface = opaqueinterface()))
    # create some 'white' light
    rays = Vector{OpticalRay{Float64,3}}(undef, 0)
    for i in 0:7
        λ = ((i / 7) * 200 + 450) / 1000
        r = OpticalRay(SVector(0.0, -3.0, 10.0), SVector(0.0, 0.5, -1.0), 1.0, λ)
        push!(rays, r)
    end
    raygen = RayListSource(rays)
    # draw the result
    Vis.drawtracerays(sys, raygenerator = raygen, test = true, trackallrays = true)
end

function fresnel(convex = true; kwargs...)
    lens = FresnelLens(Examples_N_BK7, 0.0, convex ? 15.0 : -15.0, 1.0, 8.0, 0.8, conic = 0.1)
    sys = CSGOpticalSystem(LensAssembly(lens()), Rectangle(15.0, 15.0, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, -25.0), interface = opaqueinterface()))
    Vis.drawtracerays(sys; test = true, trackallrays = true, numdivisions = 30, kwargs...)
end


function eyetrackHOE(nrays = 5000, det = false, showhead = true, zeroorder = false; kwargs...)
    # TODO update for new specs from Chris
    hoehalfwidth = 50.0 #25.0
    hoehalfheight = 50.0 #22.5
    hoecenter = SVector(-8.0 - 25.0, 0.0, -10.0 - 25.0)
    rect = Rectangle(hoehalfheight, hoehalfwidth, SVector(0.0, 1.0, 0.0), hoecenter)
    er = 15.0
    cornea_rad = 7.85
    corneavertex = SVector(0.0, er, 0.0)
    sourceloc = SVector(-33.0, er, 0.0)
    camloc = SVector(20.0, 3.0, -11.0)
    camdir = corneavertex - camloc
    camdir_norm = normalize(camdir)

    interfaces = []

    # offset = SVector(-5.0, 10.0, -10.0)
    # for θ in 0:(π / 6):(2π)
    #     ledloc = SVector(20 * cos(θ) + offset[1], 0 + offset[2], 15 * sin(θ) + offset[3])
    #     int = HologramInterface(ledloc, ConvergingBeam, sourceloc, DivergingBeam, 0.78, 100.0, Air, Examples_N_BK7, Air, Air, Air, 0.05, zeroorder)
    #     push!(interfaces, int)
    # end

    dirs = [SVector(0.7713, 0.6350, -0.0437), SVector(0.5667, 0.8111, -0.1445), SVector(0.3400, 0.9349, -0.1017), SVector(0.1492, 0.9878, 0.0445), SVector(0.0249, 0.9686, 0.2474), SVector(-0.0184, 0.8855, 0.4643), SVector(0.0254, 0.7537, 0.6567), SVector(0.1548, 0.5964, 0.7876), SVector(0.3570, 0.4462, 0.8207), SVector(0.5959, 0.3470, 0.7242), SVector(0.7976, 0.3449, 0.4948), SVector(0.8680, 0.4555, 0.1978)]

    for d in dirs
        int = HologramInterface(normalize(d), CollimatedBeam, sourceloc, DivergingBeam, 0.78, 100.0, Air, Examples_N_BK7, Air, Air, Air, 0.05, zeroorder)
        # int = HologramInterface(corneavertex - 10 * d, ConvergingBeam, sourceloc, DivergingBeam, 0.78, 100.0, Air, Examples_N_BK7, Air, Air, Air, 0.05, zeroorder)
        push!(interfaces, int)
    end

    mint = MultiHologramInterface(interfaces...)
    obj = MultiHologramSurface(rect, mint)
    cornea = leaf(Sphere(cornea_rad, interface = FresnelInterface{Float64}(EYE.CORNEA, Air, reflectance = 1.0, transmission = 0.0)), translation(0.0, er + cornea_rad, 0.0))()

    # cam settings
    fnum = 2.0
    fov = 80
    sensorrad = 1.0
    barrellength = sensorrad / tand(fov / 2)
    aprad = barrellength / fnum / 2
    camrad = max(sensorrad, aprad)

    camap = Annulus(aprad, camrad, camdir_norm, camloc)
    distfromcamtoeye = norm(camdir)
    focallength = 1 / (1 / distfromcamtoeye + 1 / barrellength)
    camlens = ParaxialLensEllipse(focallength, aprad, aprad, -camdir_norm, camloc)
    barrelloc = camloc - barrellength / 2 * camdir_norm
    barreltop = Plane(camdir_norm, camloc)
    barrelbot = Plane(-camdir_norm, camloc - 3 * barrellength * camdir_norm)
    barrelrot = OpticSim.rotmatbetween(SVector(0.0, 0.0, 1.0), camdir_norm)
    cambarrel = (
        barrelbot ∩
        barreltop ∩
        leaf(Cylinder(camrad, barrellength, interface = opaqueinterface(Float64)), Transform(barrelrot, barrelloc))
    )()
    camdet = Circle(sensorrad, camdir_norm, camloc - barrellength * camdir_norm, interface = opaqueinterface(Float64))

    # sourceleft = hoecenter[1] + hoehalfwidth - sourceloc[1]
    # sourceright = hoecenter[1] - hoehalfwidth - sourceloc[1]
    # sourceleftθ = atan(sourceleft, sourceloc[2])
    # sourcerightθ = atan(sourceright, sourceloc[2])
    # midθ = (sourceleftθ + sourcerightθ) / 2
    # sourcedir = normalize(SVector(er * tan(midθ), -er, 0.0))
    # sourceextentθ = abs(midθ - sourcerightθ)
    # source = CosineOpticalSource(RandomSource(OriginPoint{Float64}(1, position = sourceloc, direction = sourcedir), nrays, sourceextentθ), 1.0, 0.78)

    rays = Vector{OpticalRay{Float64,3}}(undef, nrays)
    @simd for i in 1:nrays
        p = point(rect, rand() * 2 - 1, rand() * 2 - 1)
        rays[i] = OpticalRay(sourceloc, p - sourceloc, 1.0, 0.78)
    end
    source = RayListSource(rays)

    sys = CSGOpticalSystem(LensAssembly(obj, cornea, camlens, cambarrel, camap), camdet, 800, 800)
    if det
        Vis.show(OpticSim.traceMT(sys, source))
    else
        Vis.drawtracerays(sys; raygenerator = source, trackallrays = true, kwargs...)
        # for θ in 0:(π / 6):(2π)
        #     ledloc = SVector(20 * cos(θ) + offset[1], 0 + offset[2], 15 * sin(θ) + offset[3])
        #     Vis.draw!(leaf(Sphere(1.0), translation(ledloc...)), color = :red)
        # end
        for d in dirs
            # Vis.draw!(leaf(Sphere(1.0), translation((corneavertex - 10 * d)...)), color = :red)
            Vis.draw!((corneavertex - 50 * d, corneavertex), color = :red)
        end
        if showhead
            Vis.draw!(joinpath(@__DIR__, "../../OBJ/glasses.obj"), scale = 100.0, transform = Transform(OpticSim.rotmatd(90, 0, 0), [27.0, 45.0, -8.0]), color = :black)
            Vis.draw!(joinpath(@__DIR__, "../../OBJ/femalehead.obj"), scale = 13.0, transform = Transform(OpticSim.rotmatd(0, 0, 180), [27.0, 105.0, -148.0]), color = :white)
        end
        Vis.display()
    end
end
