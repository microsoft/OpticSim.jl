# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

export cooketriplet, doubleconvexlensonly




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

            Material = [Air, Examples_N_BK7, Air, missing],
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
            SurfaceType = ["Object", "Standard", "Standard", "Image"],
            Radius = [Inf64, Inf64, Inf64, Inf64],
            Thickness = [Inf64, 10.0, 57.8, missing],
            Material = [Air, Examples_N_BK7, Air, missing],
            SemiDiameter = [Inf64, 9.0, 9.0, 15.0]
        )
    )
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
    raygen = Emitters.Sources.RayListSource(rays)
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
    barrelrot = OpticSim.Geometry.rotmatbetween(SVector(0.0, 0.0, 1.0), camdir_norm)
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
    source = Emitters.Sources.RayListSource(rays)

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
