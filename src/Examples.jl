# MIT License

# Copyright (c) Microsoft Corporation.

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE

"""Contains example usage of the features in the OpticSim.jl package."""
module Examples
using ..OpticSim
using ..OpticSim.Vis
# using ..OpticSim.GlassCat use this if you want to type SCHOTT.N_BK7 rather than OpticSim.GlassCat.SCHOTT.N_BK7
using StaticArrays
using DataFrames
using Images
using Unitful
using Plots
using LinearAlgebra

# Create a geometric hemisphere
function hemisphere()::CSGTree
    sph = Sphere(10.0)
    pln = Plane(0.0, 0.0, -1.0, 0.0, 0.0, 0.0)
    csgintersection(sph, pln)() #csg operations create a csggenerator which instantiates the csg tree after applying a rigid body transformation. This allows you to make as many instances of the object as you want with different transformations. We just want the CSGTree object rather than a generator.
end

# Create an optical hemisphere that has optical material properties so it will reflect and refract light. In the previous example the hemisphere object had optical properties of OpticSim.GlassCat.Air, which is the default optical interface, so it won't refract or reflect light.
function opticalhemisphere()::CSGOpticalSystem
    sph = Sphere(10.0, interface = FresnelInterface{Float64}(OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air))
    pln = Plane(0.0, 0.0, -1.0, 0.0, 0.0, 0.0, interface = FresnelInterface{Float64}(OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air))
    assy = LensAssembly{Float64}(csgintersection(sph, pln)())
    return CSGOpticalSystem(assy, Rectangle(1.0, 1.0, SVector{3,Float64}(0.0, 0.0, 1.0), SVector{3,Float64}(0.0, 0.0, -11.0)))
end

#! format: off

cooketriplet(::Type{T} = Float64, detpix::Int = 1000) where {T<:Real} = AxisymmetricOpticalSystem{T}(
    DataFrame(Surface = [:Object, 1, 2, 3, :Stop, 5, 6, :Image],
              Radius = [Inf, 26.777, 66.604, -35.571, 35.571, 35.571, -26.777, Inf],
              OptimizeRadius = [false,true,true,true,true,true,true,false],
              Thickness = [Inf, 4.0, 2.0, 4.0, 2.0, 4.0, 44.748, missing],
              OptimizeThickness = [false,true,true,true,true,true,true,false],
              Material = [OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_SK16, OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_SF2, OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_SK16, OpticSim.GlassCat.Air, missing],
              SemiDiameter = [Inf, 8.580, 7.513, 7.054, 6.033, 7.003, 7.506, 15.0]), detpix, detpix)
export cooketriplet

#no longer works
cooketripletlensonly(::Type{T} = Float64) where {T<:Real} = AxisymmetricLens{T}(
    DataFrame(Surface = [:Object, 1, 2, 3, :Stop, 5, 6, :Image],
              Radius = [Inf, 26.777, 66.604, -35.571, 35.571, 35.571, -26.777, Inf],
              OptimizeRadius = [false,true,true,true,true,true,true,false],
              Thickness = [Inf, 4.0, 2.0, 4.0, 2.0, 4.0, 44.748, missing],
              OptimizeThickness = [false,true,true,true,true,true,true,false],
              Material = [OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_SK16, OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_SF2, OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_SK16, OpticSim.GlassCat.Air, missing],
              SemiDiameter = [Inf, 8.580, 7.513, 7.054, 6.033, 7.003, 7.506, 15.0]))
export cooketripletlensonly

cooketripletfirstelement(::Type{T} = Float64) where {T<:Real} = AxisymmetricOpticalSystem(
    DataFrame(Surface = [:Object, 1, 2, :Image],
              Radius = [Inf, -35.571, 35.571, Inf],
              Thickness = [Inf, 4.0, 44.748, missing],
              Material = [OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_SK16, OpticSim.GlassCat.Air, missing],
              SemiDiameter = [Inf, 7.054, 6.033, 15.0]))

convexplano(::Type{T} = Float64) where {T<:Real} = AxisymmetricOpticalSystem{T}(
    DataFrame(Surface = [:Object, 1, 2, :Image],
              Radius = [Inf, 60.0, Inf, Inf],
              Thickness = [Inf, 10.0, 57.8, missing],
              Material = [OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, missing],
              SemiDiameter = [Inf, 9.0, 9.0, 15.0]))

doubleconvex(frontradius::T,rearradius::T) where{T<:Real} =
AxisymmetricOpticalSystem{T}(
    DataFrame(Surface = [:Object, 1, 2, :Image],
              Radius = [T(Inf64), frontradius, rearradius, T(Inf64)],
              OptimizeRadius = [false,true,true,false],
              Thickness = [T(Inf64), T(10.0), T(57.8), missing],
              OptimizeThickness = [false,false,false,false],
              Material = [OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, missing],
              SemiDiameter = [T(Inf64), T(9.0), T(9.0), T(15.0)]))

doubleconvexconic(::Type{T} = Float64) where {T<:Real} =
              AxisymmetricOpticalSystem{T}(
                  DataFrame(Surface = [:Object, 1, 2, :Image],
                            Radius = [Inf64, 60, -60, Inf64],
                            OptimizeRadius = [false,true,true,false],
                            Thickness = [Inf64, 10.0, 57.8, missing],
                            OptimizeThickness = [false,false,false,false],
                            Conic = [missing, 0.01, 0.01, missing],
                            OptimizeConic = [false, true, true, false],
                            Material = [OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, missing],
                            SemiDiameter = [Inf64, 9.0, 9.0, 15.0]))

doubleconvexlensonly(frontradius::T,rearradius::T) where{T<:Real} =
AxisymmetricLens{T}(
    DataFrame(Surface = [:Object, 1, 2, :Image],
              Radius = [T(Inf64), frontradius, rearradius, T(Inf64)],
              OptimizeRadius = [false,true,true,false],
              Thickness = [T(Inf64), T(10.0), T(57.8), missing],
              OptimizeThickness = [false,false,false,false],
              Material = [OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, missing],
              SemiDiameter = [T(Inf64), T(9.0), T(9.0), T(15.0)]))
export doubleconvexlensonly

doubleconvexprescription() = DataFrame(
    Surface = [:Object, 1, 2, :Image],
    Radius = [Inf64, 60, -60, Inf64],
    OptimizeRadius = [false,true,true,false],
    Thickness = [Inf64, 10.0, 57.8, missing],
    OptimizeThickness = [false,true,true,false],
    Material = [OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, missing],
    SemiDiameter = [Inf64, 9.0, 9.0, 15.0])

doubleconvex(::Type{T} = Float64; temperature::Unitful.Temperature = OpticSim.GlassCat.TEMP_REF_UNITFUL, pressure::T = T(OpticSim.GlassCat.PRESSURE_REF)) where {T<:Real} = AxisymmetricOpticalSystem{T}(doubleconvexprescription(),temperature = temperature, pressure = pressure)

doubleconcave(::Type{T} = Float64) where {T<:Real} = AxisymmetricOpticalSystem{T}(
    DataFrame(Surface = [:Object, 1, 2, :Image],
              Radius = [Inf64, -41.0, 41.0, Inf64],
              Thickness = [Inf64, 10.0, 57.8, missing],
              Material = [OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, missing],
              SemiDiameter = [Inf64, 9.0, 9.0, 15.0]))

planoconcaverefl(::Type{T} = Float64) where {T<:Real} = AxisymmetricOpticalSystem{T}(
    DataFrame(Surface = [:Object, 1, 2, :Image],
              Radius = [Inf64, Inf64, -41.0, Inf64],
              Thickness = [Inf64, 10.0, -57.8, missing],
              Material = [OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, missing],
              SemiDiameter = [Inf64, 9.0, 9.0, 25.0],
              Reflectance = [missing, missing, 1.0, missing]))

concaveplano(::Type{T} = Float64) where {T<:Real} = AxisymmetricOpticalSystem{T}(
    DataFrame(Surface = [:Object, 1, 2, :Image],
              Radius = [Inf64, -41.0, Inf64, Inf64],
              Thickness = [Inf64, 10.0, 57.8, missing],
              Material = [OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, missing],
              SemiDiameter = [Inf64, 9.0, 9.0, 15.0]))

planoplano(::Type{T} = Float64) where {T<:Real} = AxisymmetricOpticalSystem{T}(
    DataFrame(Surface = [:Object, 1, 2, :Image],
              Radius = [Inf64, Inf64, Inf64, Inf64],
              Thickness = [Inf64, 10.0, 57.8, missing],
              Material = [OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, missing],
              SemiDiameter = [Inf64, 9.0, 9.0, 15.0]))

#! format: on

end #module SphericalLenses




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

function SchmidtCassegrainTelescope()
    # glass entrance lens on telescope
    topsurf = Plane(SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, 0.0), interface = FresnelInterface{Float64}(OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air), vishalfsizeu = 12.00075, vishalfsizev = 12.00075)
    botsurf = AcceleratedParametricSurface(ZernikeSurface(12.00075, radius = -1.14659768e+4, aspherics = [(4, 3.68090959e-7), (6, 2.73643352e-11), (8, 3.20036892e-14)]), 17, interface = FresnelInterface{Float64}(OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air))
    coverlens = csgintersection(leaf(Cylinder(12.00075, 1.4)), csgintersection(leaf(topsurf), leaf(botsurf, RigidBodyTransform(OpticSim.rotmatd(0, 180, 0), SVector(0.0, 0.0, -0.65)))))
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

    return CSGOpticalSystem(la, det)
end

drawSchmidt(; kwargs...) = Vis.drawtracerays(SchmidtCassegrainTelescope(), raygenerator = UniformOpticalSource(CollimatedSource(GridRectOriginPoints(5, 5, 10.0, 10.0, position = SVector(0.0, 0.0, 20.0))), 0.55), trackallrays = true, colorbynhits = true, test = true, numdivisions = 100; kwargs...)

function prism_refraction()
    # build the triangular prism
    int = FresnelInterface{Float64}(OpticSim.GlassCat.SCHOTT.N_SF14, OpticSim.GlassCat.Air)
    s = 2.0
    prism = csgintersection(leaf(Plane(SVector(0.0, -1.0, 0.0), SVector(0.0, -s, 0.0), interface = int, vishalfsizeu = 2 * s, vishalfsizev = 2 * s)), csgintersection(Plane(SVector(0.0, sind(30), cosd(30)), SVector(0.0, s * sind(30), s * cosd(30)), interface = int, vishalfsizeu = 2 * s, vishalfsizev = 2 * s), Plane(SVector(0.0, sind(30), -cosd(30)), SVector(0.0, s * sind(30), -s * cosd(30)), interface = int, vishalfsizeu = 2 * s, vishalfsizev = 2 * s)))
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

function zoom_lens(pos = 1)
    if pos == 0
        stop = 2.89
        zoom = 9.48
        dist = 4.46970613
    elseif pos == 1
        stop = 3.99
        zoom = 4.48
        dist = 21.21
    else
        stop = 4.90
        zoom = 2.00
        dist = 43.81
    end
    #! format: off
    return AxisymmetricOpticalSystem{Float64}(
        DataFrame(Surface = [:Object, :Stop, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, :Image],
                  Radius = [Inf64, Inf64, -1.6202203499676E+01, -4.8875855327468E+01, 1.5666614444619E+01, -4.2955326460481E+01, 1.0869565217391E+02, 2.3623907394283E+01, -1.6059097478722E+01, -4.2553191489362E+02, -3.5435861091425E+01, -1.4146272457208E+01, -2.5125628140704E+02, -2.2502250225023E+01, -1.0583130489999E+01, -4.4444444444444E+01, Inf64],
                  Aspherics = [missing, missing, missing, missing, missing, [(4, 1.03860000000E-04), (6, 1.42090000000E-07), (8, -8.84950000000E-09), (10, 1.24770000000E-10), (12, -1.03670000000E-12), (14, 3.65560000000E-15)], missing, missing, [(4, 4.27210000000E-05), (6, 1.24840000000E-07), (8, 9.70790000000E-09), (10, -1.84440000000E-10), (12, 1.86440000000E-12), (14, -7.79750000000E-15)], [(4, 1.13390000000E-04), (6, 4.81650000000E-07), (8, 1.87780000000E-08), (10, -5.75710000000E-10), (12, 8.99940000000E-12), (14, -4.67680000000E-14)], missing, missing, missing, missing, missing, missing, missing],
                  Thickness = [Inf64, 0.0, 5.18, 0.10, 4.40, 0.16, 1.0, 4.96, zoom, 4.04, 1.35, 1.0, 2.80, 3.0, 1.22, dist, missing],
                  Material = [OpticSim.GlassCat.Air, OpticSim.GlassCat.Air, OpticSim.GlassCat.OHARA.S_LAH66, OpticSim.GlassCat.Air, OpticSim.GlassCat.NIKON.LLF6, OpticSim.GlassCat.Air, OpticSim.GlassCat.OHARA.S_TIH6, OpticSim.GlassCat.OHARA.S_FSL5, OpticSim.GlassCat.Air, OpticSim.GlassCat.OHARA.S_FSL5, OpticSim.GlassCat.Air, OpticSim.GlassCat.OHARA.S_LAL8, OpticSim.GlassCat.SCHOTT.S_FL4, OpticSim.GlassCat.Air, OpticSim.GlassCat.OHARA.S_LAH66, OpticSim.GlassCat.Air, missing],
                  SemiDiameter = [Inf64, stop, 3.85433218451, 3.85433218451, 4.36304692871, 4.36304692871, 4.72505505439, 4.72505505439, 4.72505505439, 4.45240784026, 4.45240784026, 4.50974054117, 4.50974054117, 4.50974054117, 4.76271114409, 4.76271114409, 15.0]))
    #! format: on

end

function fresnel(convex = true; kwargs...)
    lens = FresnelLens(OpticSim.GlassCat.SCHOTT.N_BK7, 0.0, convex ? 15.0 : -15.0, 1.0, 8.0, 0.8, conic = 0.1)
    sys = CSGOpticalSystem(LensAssembly(lens()), Rectangle(15.0, 15.0, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, -25.0), interface = opaqueinterface()))
    Vis.drawtracerays(sys; test = true, trackallrays = true, numdivisions = 30, kwargs...)
end

function grating(; period = 1.0, θ = 0.0, λ = 0.55, kwargs...)
    int = ThinGratingInterface(SVector(0.0, 1.0, 0.0), period, OpticSim.GlassCat.Air, OpticSim.GlassCat.Air, minorder = -2, maxorder = 2, reflectance = [0.0, 0.0, 0.1, 0.0, 0.0], transmission = [0.05, 0.1, 0.4, 0.1, 0.05])
    grating = ThinGratingSurface(Rectangle(5.0, 5.0, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, 0.0)), int)
    back = Rectangle(30.0, 30.0, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, 25.0))
    sys = CSGOpticalSystem(LensAssembly(grating, back), Rectangle(30.0, 30.0, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, -25.0), interface = opaqueinterface()))
    Vis.drawtracerays(sys; raygenerator = UniformOpticalSource(CollimatedSource(OriginPoint{Float64}(100, position = SVector(0.0, 0.0, 10.0), direction = SVector(0.0, sind(θ), -cosd(θ)))), λ), trackallrays = true, rayfilter = nothing, kwargs...)
end

function reflgrating(; period = 1.0, θ = 0.0, λ = 0.55, kwargs...)
    int = ThinGratingInterface(SVector(0.0, 1.0, 0.0), period, OpticSim.GlassCat.Air, OpticSim.GlassCat.Air, minorder = -2, maxorder = 2, transmission = [0.0, 0.0, 0.1, 0.0, 0.0], reflectance = [0.05, 0.1, 0.4, 0.1, 0.05])
    grating = ThinGratingSurface(Rectangle(5.0, 5.0, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, 0.0)), int)
    back = Rectangle(30.0, 30.0, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, -25.0))
    sys = CSGOpticalSystem(LensAssembly(grating, back), Rectangle(30.0, 30.0, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, 25.0), interface = opaqueinterface()))
    Vis.drawtracerays(sys; raygenerator = UniformOpticalSource(CollimatedSource(OriginPoint{Float64}(100, position = SVector(0.0, 0.0, 10.0), direction = SVector(0.0, sind(θ), -cosd(θ)))), λ), trackallrays = true, rayfilter = nothing, kwargs...)
end

function HOE(refl = false, firstorderonly = false; kwargs...)
    rect = Rectangle(5.0, 5.0, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, 0.0))
    if refl
        int = HologramInterface(SVector(0.0, -10.0, 20.0), ConvergingBeam, SVector(0.0, 0.0, -200), ConvergingBeam, 0.55, 9.0, OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, OpticSim.GlassCat.Air, OpticSim.GlassCat.Air, 0.05, !firstorderonly)
    else
        int = HologramInterface(SVector(0.0, -10.0, -20.0), ConvergingBeam, SVector(0.0, 0.0, -200), ConvergingBeam, 0.55, 5.0, OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, OpticSim.GlassCat.Air, OpticSim.GlassCat.Air, 0.05, !firstorderonly)
    end
    obj = HologramSurface(rect, int)
    back = Rectangle(50.0, 50.0, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, 25.0))
    sys = CSGOpticalSystem(LensAssembly(obj, back), Rectangle(50.0, 50.0, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, -25.0), interface = opaqueinterface()))
    Vis.drawtracerays(sys; raygenerator = UniformOpticalSource(GridSource(OriginPoint{Float64}(10, position = SVector(0.0, 0.0, 10.0), direction = SVector(0.0, 0.0, -1.0)), 1, 15, 0.0, π / 6), 0.55), trackallrays = true, rayfilter = nothing, kwargs...)
end

function HOEfocus(; kwargs...)
    rect = Rectangle(5.0, 5.0, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, 0.0))
    int = HologramInterface(SVector(0.0, -3.0, -20.0), ConvergingBeam, SVector(0.0, 0.0, -1.0), CollimatedBeam, 0.55, 9.0, OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, OpticSim.GlassCat.Air, OpticSim.GlassCat.Air, 0.05, false)
    obj = HologramSurface(rect, int)
    sys = CSGOpticalSystem(LensAssembly(obj), Rectangle(10.0, 10.0, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, -25.0), interface = opaqueinterface()))
    Vis.drawtracerays(sys; raygenerator = UniformOpticalSource(CollimatedSource(GridRectOriginPoints(5, 5, 3.0, 3.0, position = SVector(0.0, 0.0, 10.0), direction = SVector(0.0, 0.0, -1.0))), 0.55), trackallrays = true, rayfilter = nothing, test = true, kwargs...)
end

function HOEcollimate(; kwargs...)
    rect = Rectangle(5.0, 5.0, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, 0.0))
    int = HologramInterface(SVector(0.1, -0.05, -1.0), CollimatedBeam, SVector(0.0, 0.0, 10), DivergingBeam, 0.55, 9.0, OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, OpticSim.GlassCat.Air, OpticSim.GlassCat.Air, 0.05, false)
    obj = HologramSurface(rect, int)
    sys = CSGOpticalSystem(LensAssembly(obj), Rectangle(10.0, 10.0, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, -25.0), interface = opaqueinterface()))
    Vis.drawtracerays(sys; raygenerator = UniformOpticalSource(GridSource(OriginPoint{Float64}(1, position = SVector(0.0, 0.0, 10.0), direction = SVector(0.0, 0.0, -1.0)), 5, 5, π / 4, π / 4), 0.55), trackallrays = true, rayfilter = nothing, test = true, kwargs...)
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
    #     int = HologramInterface(ledloc, ConvergingBeam, sourceloc, DivergingBeam, 0.78, 100.0, OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, OpticSim.GlassCat.Air, OpticSim.GlassCat.Air, 0.05, zeroorder)
    #     push!(interfaces, int)
    # end

    dirs = [SVector(0.7713, 0.6350, -0.0437), SVector(0.5667, 0.8111, -0.1445), SVector(0.3400, 0.9349, -0.1017), SVector(0.1492, 0.9878, 0.0445), SVector(0.0249, 0.9686, 0.2474), SVector(-0.0184, 0.8855, 0.4643), SVector(0.0254, 0.7537, 0.6567), SVector(0.1548, 0.5964, 0.7876), SVector(0.3570, 0.4462, 0.8207), SVector(0.5959, 0.3470, 0.7242), SVector(0.7976, 0.3449, 0.4948), SVector(0.8680, 0.4555, 0.1978)]

    for d in dirs
        int = HologramInterface(normalize(d), CollimatedBeam, sourceloc, DivergingBeam, 0.78, 100.0, OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, OpticSim.GlassCat.Air, OpticSim.GlassCat.Air, 0.05, zeroorder)
        # int = HologramInterface(corneavertex - 10 * d, ConvergingBeam, sourceloc, DivergingBeam, 0.78, 100.0, OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, OpticSim.GlassCat.Air, OpticSim.GlassCat.Air, 0.05, zeroorder)
        push!(interfaces, int)
    end

    mint = MultiHologramInterface(interfaces...)
    obj = MultiHologramSurface(rect, mint)
    cornea = leaf(Sphere(cornea_rad, interface = FresnelInterface{Float64}(OpticSim.GlassCat.EYE.CORNEA, OpticSim.GlassCat.Air, reflectance = 1.0, transmission = 0.0)), translation(0.0, er + cornea_rad, 0.0))()

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
    cambarrel = csgintersection(barrelbot, csgintersection(barreltop, leaf(Cylinder(camrad, barrellength, interface = opaqueinterface(Float64)), RigidBodyTransform(barrelrot, barrelloc))))()
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
            Vis.draw!(joinpath(@__DIR__, "../../OBJ/glasses.obj"), scale = 100.0, transform = RigidBodyTransform(OpticSim.rotmatd(90, 0, 0), [27.0, 45.0, -8.0]), color = :black)
            Vis.draw!(joinpath(@__DIR__, "../../OBJ/femalehead.obj"), scale = 13.0, transform = RigidBodyTransform(OpticSim.rotmatd(0, 0, 180), [27.0, 105.0, -148.0]), color = :white)
        end
        Vis.display()
    end
end

function multiHOE(; kwargs...)
    rect = Rectangle(5.0, 5.0, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, 0.0))
    int1 = HologramInterface(SVector(-5.0, 0.0, -20.0), ConvergingBeam, SVector(0.0, -1.0, -1.0), CollimatedBeam, 0.55, 100.0, OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, OpticSim.GlassCat.Air, OpticSim.GlassCat.Air, 0.05, false)
    int2 = HologramInterface(SVector(5.0, 0.0, -20.0), ConvergingBeam, SVector(0.0, 1.0, -1.0), CollimatedBeam, 0.55, 100.0, OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, OpticSim.GlassCat.Air, OpticSim.GlassCat.Air, 0.05, false)
    mint = MultiHologramInterface(int1, int2)
    obj = MultiHologramSurface(rect, mint)
    sys = CSGOpticalSystem(LensAssembly(obj), Rectangle(10.0, 10.0, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, -20.0), interface = opaqueinterface()))
    s1 = UniformOpticalSource(CollimatedSource(RandomRectOriginPoints(500, 3.0, 3.0, position = SVector(0.0, 3.0, 3.0), direction = SVector(0.0, -1.0, -1.0))), 0.55, sourcenum = 1)
    s2 = UniformOpticalSource(CollimatedSource(RandomRectOriginPoints(500, 3.0, 3.0, position = SVector(0.0, -3.0, 3.0), direction = SVector(0.0, 1.0, -1.0))), 0.55, sourcenum = 2)
    s3 = UniformOpticalSource(CollimatedSource(RandomRectOriginPoints(500, 3.0, 3.0, position = SVector(0.0, 0.0, 3.0), direction = SVector(0.0, 0.0, -1.0))), 0.55, sourcenum = 3)
    Vis.drawtracerays(sys; raygenerator = s1, trackallrays = true, colorbysourcenum = true, rayfilter = nothing, kwargs...)
    Vis.drawtracerays!(sys; raygenerator = s2, trackallrays = true, colorbysourcenum = true, rayfilter = nothing, drawgen = true, kwargs...)
    Vis.drawtracerays!(sys; raygenerator = s3, trackallrays = true, colorbysourcenum = true, rayfilter = nothing, drawgen = true, kwargs...)
end
export Examples
