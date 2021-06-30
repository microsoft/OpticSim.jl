# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

"""
Contains all complex data used for testing and benchmarking.
"""
module TestData

using OpticSim
using OpticSim: tobeziersegments # try to use only exported functions so this list should stay short
using OpticSim.Geometry
using OpticSim.Emitters
using StaticArrays
using DataFrames
using Unitful
using Images
using Base: @.
using LinearAlgebra


### CURVES

function bsplinecurve()
    knots = KnotVector{Float64}([0.0, 0.0, 0.0, 0.0, 0.3, 0.6, 1.0, 1.0, 1.0, 1.0])
    points = Array{MVector{2,Float64},1}([(0.05, 0.05), (0.05, 0.8), (0.6, 0.0), (0.7, 0.5), (0.5, 0.95), (0.95, 0.5)])
    curve = BSplineCurve{OpticSim.Euclidean,Float64,2,3}(knots, points)
    return curve
end

function homogeneousbsplinecurve()
    knots = KnotVector{Float64}([0.0, 0.0, 0.0, 0.0, 0.3, 0.6, 1.0, 1.0, 1.0, 1.0])
    points = Array{MVector{3,Float64},1}([(0.05, 0.0, 1.0), (0.05, 0.8, 1.0), (0.6, 0.0, 1.0), (0.7, 0.5, 1.0), (0.5, 0.95, 1.0), (0.95, 0.5, 1.0)])
    curve = BSplineCurve{OpticSim.Rational,Float64,3,3}(knots, points)
    return curve
end

function beziercurve()
    orig = onespanspline()
    return BezierCurve{OpticSim.Euclidean,Float64,2,3}(tobeziersegments(orig)[1])
end

function onespanspline()
    knots = KnotVector{Float64}([0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0])
    points = Array{MVector{2,Float64},1}([(0.05, 0.05), (0.05, 0.95), (0.5, 0.95), (0.95, 0.5)])
    curve = BSplineCurve{OpticSim.Euclidean,Float64,2,3}(knots, points)
    return curve
end

### SURFACES

function bsplinesurface()
    vknots = KnotVector{Float64}([0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0])
    uknots = KnotVector{Float64}([0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0])
    points = map(
        x -> collect(x),
        [
            (0.0, 0.0, 0.0) (0.33, 0.0, 0.0) (0.66, 0.0, 0.0) (1.0, 0.0, 0.0) (1.33, 0.0, 0.0)
            (0.0, 0.33, 0.0) (0.33, 0.33, 1.0) (0.66, 0.33, 1.0) (1.0, 0.33, 0.5) (1.33, 0.33, 0.0)
            (0.0, 0.66, 0.0) (0.33, 0.66, 1.0) (0.66, 0.66, 1.0) (1.0, 0.66, 0.5) (1.33, 0.66, 0.0)
            (0.0, 1.0, 0.0) (0.33, 1.0, 0.0) (0.66, 1.0, 0.0) (1.0, 1.0, 0.0) (1.33, 1.0, 0.0)
        ],
    )
    return BSplineSurface{OpticSim.Euclidean,Float64,3,3}(uknots, vknots, points)
end

function beziersurface()
    # only simple way I could think of to write this array literal as 2D array of 1D arrays.
    # Julia doesn't parse a 2D array of 1D arrays properly. Maybe should use tuples instead though.
    points = map(
        x -> collect(x),
        [
            (0.0, 0.0, 0.0) (0.0, 0.33, 0.0) (0.0, 0.66, 0.0) (0.0, 1.0, 0.0)
            (0.33, 0.0, 0.0) (0.33, 0.33, 1.0) (0.33, 0.66, 1.0) (0.33, 1.0, 0.0)
            (0.66, 0.0, 0.0) (0.66, 0.33, 1.0) (0.66, 0.66, 1.0) (0.66, 1.0, 0.0)
            (1.0, 0.0, 0.0) (1.0, 0.33, 0.0) (1.0, 0.66, 0.0) (1.0, 1.0, 0.0)
        ],
    )
    return BezierSurface{OpticSim.Euclidean,Float64,3,3}(points)
end

function upsidedownbeziersurface()
    points = map(
        x -> collect(x),
        [
            (0.0, 0.0, 0.0) (0.33, 0.0, 0.0) (0.66, 0.0, 0.0) (1.0, 0.0, 0.0)
            (0.0, 0.33, 0.0) (0.33, 0.33, -1.0) (0.66, 0.33, -1.0) (1.0, 0.33, -.5)
            (0.0, 0.66, 0.0) (0.33, 0.66, -1.0) (0.66, 0.66, -1.0) (1.0, 0.66, -.5)
            (0.0, 1.0, 0.0) (0.33, 1.0, 0.0) (0.66, 1.0, 0.0) (1.0, 1.0, 0.0)
        ],
    )
    return BezierSurface{OpticSim.Euclidean,Float64,3,3}(points)
end

function wavybeziersurface()
    points = map(
        x -> collect(x),
        [
            (0.0, 0.0, -0.3) (0.0, 0.33, 1.0) (0.0, 0.66, -1.0) (0.0, 1.0, 0.3)
            (0.33, 0.0, -0.3) (0.33, 0.33, 1.0) (0.33, 0.66, -1.0) (0.33, 1.0, 0.3)
            (0.66, 0.0, -0.3) (0.66, 0.33, 1.0) (0.66, 0.66, -1.0) (0.66, 1.0, 0.3)
            (1.0, 0.0, -0.3) (1.0, 0.33, 1.0) (1.0, 0.66, -1.0) (1.0, 1.0, 0.3)
        ],
    )
    return BezierSurface{OpticSim.Euclidean,Float64,3,3}(points)
end

function verywavybeziersurface()
    points = map(
        x -> collect(x),
        [
            (0.0, 0.0, 0.5) (0.0, 0.2, 1.0) (0.0, 0.4, -3.0) (0.0, 0.6, 3.0) (0.0, 0.8, -1.0) (0.0, 1.0, 0.5)
            (0.2, 0.0, 0.5) (0.2, 0.2, 1.0) (0.2, 0.4, -3.0) (0.2, 0.6, 3.0) (0.2, 0.8, -1.0) (0.2, 1.0, 0.5)
            (0.4, 0.0, 0.0) (0.4, 0.2, 1.0) (0.4, 0.4, -3.0) (0.4, 0.6, 3.0) (0.4, 0.8, -1.0) (0.4, 1.0, 0.5)
            (0.6, 0.0, 0.0) (0.6, 0.2, 1.0) (0.6, 0.4, -3.0) (0.6, 0.6, 3.0) (0.6, 0.8, -1.0) (0.6, 1.0, 0.5)
            (0.8, 0.0, -0.5) (0.8, 0.2, 1.0) (0.8, 0.4, -3.0) (0.8, 0.6, 3.0) (0.8, 0.8, -1.0) (0.8, 1.0, 0.5)
            (1.0, 0.0, -0.5) (1.0, 0.2, 1.0) (1.0, 0.4, -3.0) (1.0, 0.6, 3.0) (1.0, 0.8, -1.0) (1.0, 1.0, 0.5)
        ],
    )
    return BezierSurface{OpticSim.Euclidean,Float64,3,5}(points)
end

qtypesurface1() = QTypeSurface(1.5, radius = 5.0, conic = 1.0, αcoeffs = [(0, 1, -0.3), (5, 4, 0.1)], βcoeffs = [(3, 7, 0.1)], normradius = 1.6)

zernikesurface1() = ZernikeSurface(1.5, radius = 5.0, conic = 1.0, zcoeff = [(1, 0.0), (4, -0.1), (7, 0.05), (33, 0.03)], aspherics = [(2, 0.0), (10, 0.01)])

zernikesurface1a() = ZernikeSurface(1.5, radius = 5.0, conic = 1.0, zcoeff = [(1, 0.0), (4, -0.1), (7, 0.05), (33, 0.03)], aspherics = [(2, 0.0), (10, 0.01)], normradius = 1.6)

zernikesurface2() = ZernikeSurface(1.5, radius = 5.0, conic = 1.0, zcoeff = [(7, 0.2), (9, 0.2)], aspherics = [(10, -0.01)])

zernikesurface3() = ZernikeSurface(9.5, zcoeff = [(1, 0.0), (4, -1.0), (7, 0.5)])

conicsurface() = AcceleratedParametricSurface(ZernikeSurface(9.5, radius = -15.0, conic = -5.0))

function simplesaggrid()
    points = map(
        x -> collect(x),
        [
            (0.1, 0.0, 0.0, 0.0) (0.3, 0.0, 0.0, 0.0) (0.0, 0.0, 0.0, 0.0) (0.3, 0.0, 0.0, 0.0)
            (0.4, 0.0, 0.0, 0.0) (0.1, 0.0, 0.0, 0.0) (0.0, 0.0, 0.0, 0.0) (0.1, 0.0, 0.0, 0.0)
            (-0.3, 0.0, 0.0, 0.0) (0.1, 0.0, 0.0, 0.0) (0.2, 0.0, 0.0, 0.0) (0.3, 0.0, 0.0, 0.0)
            (0.2, 0.0, 0.0, 0.0) (0.3, 0.0, 0.0, 0.0) (-0.1, 0.0, 0.0, 0.0) (-0.1, 0.0, 0.0, 0.0)
            (0.1, 0.0, 0.0, 0.0) (0.3, 0.0, 0.0, 0.0) (0.0, 0.0, 0.0, 0.0) (0.3, 0.0, 0.0, 0.0)
            (0.4, 0.0, 0.0, 0.0) (0.1, 0.0, 0.0, 0.0) (0.0, 0.0, 0.0, 0.0) (0.1, 0.0, 0.0, 0.0)
        ],
    )
    return points
end

gridsagsurfacelinear() = GridSagSurface(ZernikeSurface(1.5, radius = 5.0, conic = 1.0, zcoeff = [(1, 0.0), (4, -0.1), (7, 0.05), (33, 0.03)], aspherics = [(2, 0.0), (10, 0.01)]), simplesaggrid(), interpolation = GridSagLinear)

gridsagsurfacebicubic() = GridSagSurface(ZernikeSurface(1.5, radius = 25.0, conic = 1.0, zcoeff = [(1, 0.0), (4, -0.1), (7, 0.05), (33, 0.03)], aspherics = [(2, 0.0), (10, 0.01)]), simplesaggrid(), interpolation = GridSagBicubic)

gridsagsurfacebicubiccheby() = GridSagSurface(ChebyshevSurface(1.5, 2.0, radius = 25.0, conic = 0.2, [(0, 1, 0.03), (2, 1, -0.01)]), simplesaggrid(), interpolation = GridSagBicubic)

chebyshevsurface1() = ChebyshevSurface(2.0, 2.0, [(1, 1, 0.1), (4, 5, -0.3), (2, 1, 0.0), (22, 23, 0.01)], radius = 25.0, conic = 0.2)

chebyshevsurface2() = ChebyshevSurface(2.0, 2.0, [(1, 1, 0.1), (0, 2, 0.3), (4, 2, -0.1)], radius = 10.0, conic = 0.2)

### LENSES

function zernikelens()
    topsurface = leaf(AcceleratedParametricSurface(ZernikeSurface(9.5, zcoeff = [(1, 0.0), (4, -1.0), (7, 0.5)]), interface = FresnelInterface{Float64}(OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air)), translation(0.0, 0.0, 5.0))
    botsurface = leaf(Plane(0.0, 0.0, -1.0, 0.0, 0.0, -5.0, vishalfsizeu = 9.5, vishalfsizev = 9.5, interface = FresnelInterface{Float64}(OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air)))
    barrel = leaf(Cylinder(9.0, 20.0, interface = FresnelInterface{Float64}(OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, reflectance = zero(Float64), transmission = zero(Float64))))
    return (barrel ∩ topsurface ∩ botsurface)(Transform{Float64}(0.0, Float64(π), 0.0, 0.0, 0.0, -5.0))
end

function asphericlens()
    topsurface = leaf(AcceleratedParametricSurface(ZernikeSurface(9.5, aspherics = [(2, 0.0), (4, -0.001)]), interface = FresnelInterface{Float64}(OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air)), translation(0.0, 0.0, 5.0))
    botsurface = leaf(Plane(0.0, 0.0, -1.0, 0.0, 0.0, -5.0, vishalfsizeu = 9.5, vishalfsizev = 9.5, interface = FresnelInterface{Float64}(OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air)))
    barrel = leaf(Cylinder(9.0, 20.0, interface = FresnelInterface{Float64}(OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, reflectance = zero(Float64), transmission = zero(Float64))))
    return (barrel ∩ topsurface ∩ botsurface)(Transform{Float64}(0.0, Float64(π), 0.0, 0.0, 0.0, -5.0))
end

function coniclensZ()
    topsurface = leaf(AcceleratedParametricSurface(ZernikeSurface(9.5, radius = -15.0, conic = -5.0), interface = FresnelInterface{Float64}(OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air)), translation(0.0, 0.0, 5.0))
    botsurface = leaf(Plane(0.0, 0.0, -1.0, 0.0, 0.0, -5.0, vishalfsizeu = 9.5, vishalfsizev = 9.5, interface = FresnelInterface{Float64}(OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air)))
    barrel = leaf(Cylinder(9.0, 20.0, interface = FresnelInterface{Float64}(OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, reflectance = zero(Float64), transmission = zero(Float64))))
    return (barrel ∩ topsurface ∩ botsurface)(Transform{Float64}(0.0, Float64(π), 0.0, 0.0, 0.0, -5.0))
end

function coniclensQ()
    topsurface = leaf(AcceleratedParametricSurface(QTypeSurface(9.5, radius = -15.0, conic = -5.0), interface = FresnelInterface{Float64}(OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air)), translation(0.0, 0.0, 5.0))
    botsurface = leaf(Plane(0.0, 0.0, -1.0, 0.0, 0.0, -5.0, vishalfsizeu = 9.5, vishalfsizev = 9.5, interface = FresnelInterface{Float64}(OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air)))
    barrel = leaf(Cylinder(9.0, 20.0, interface = FresnelInterface{Float64}(OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, reflectance = zero(Float64), transmission = zero(Float64))))
    return (barrel ∩ topsurface ∩ botsurface)(Transform{Float64}(0.0, Float64(π), 0.0, 0.0, 0.0, -5.0))
end

function qtypelens()
    topsurface = leaf(AcceleratedParametricSurface(QTypeSurface(9.0, radius = -25.0, conic = 0.3, αcoeffs = [(1, 0, 0.3), (1, 1, 1.0)], βcoeffs = [(1, 0, -0.1), (2, 0, 0.4), (3, 0, -0.6)], normradius = 9.5), interface = FresnelInterface{Float64}(OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air)), translation(0.0, 0.0, 5.0))
    botsurface = leaf(Plane(0.0, 0.0, -1.0, 0.0, 0.0, -5.0, vishalfsizeu = 9.5, vishalfsizev = 9.5, interface = FresnelInterface{Float64}(OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air)))
    barrel = leaf(Cylinder(9.0, 20.0, interface = FresnelInterface{Float64}(OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, reflectance = zero(Float64), transmission = zero(Float64))))
    return (barrel ∩ topsurface ∩ botsurface)(Transform{Float64}(0.0, Float64(π), 0.0, 0.0, 0.0, -5.0))
end

function chebyshevlens()
    topsurface = leaf(AcceleratedParametricSurface(ChebyshevSurface(9.0, 9.0, [(1, 0, -0.081), (2, 0, -0.162), (0, 1, 0.243), (1, 1, 0.486), (2, 1, 0.081), (0, 2, -0.324), (1, 2, -0.405), (2, 2, -0.81)], radius = -25.0), interface = FresnelInterface{Float64}(OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air)), translation(0.0, 0.0, 5.0))
    botsurface = leaf(Plane(0.0, 0.0, -1.0, 0.0, 0.0, -5.0, vishalfsizeu = 9.5, vishalfsizev = 9.5, interface = FresnelInterface{Float64}(OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air)))
    barrel = leaf(Cylinder(9.0, 20.0, interface = FresnelInterface{Float64}(OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, reflectance = zero(Float64), transmission = zero(Float64))))
    return (barrel ∩ topsurface ∩ botsurface)(Transform{Float64}(0.0, Float64(π), 0.0, 0.0, 0.0, -5.0))
end

function gridsaglens()
    topsurface = leaf(AcceleratedParametricSurface(GridSagSurface{Float64}(joinpath(@__DIR__, "../../test/test.GRD"), radius = -25.0, conic = 0.4, interpolation = GridSagBicubic), interface = FresnelInterface{Float64}(OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air)), translation(0.0, 0.0, 5.0))
    botsurface = leaf(Plane(0.0, 0.0, -1.0, 0.0, 0.0, -5.0, vishalfsizeu = 9.5, vishalfsizev = 9.5, interface = FresnelInterface{Float64}(OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air)))
    barrel = leaf(Cylinder(9.0, 20.0, interface = FresnelInterface{Float64}(OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, reflectance = zero(Float64), transmission = zero(Float64))))
    return (barrel ∩ topsurface ∩ botsurface)(Transform{Float64}(0.0, Float64(π), 0.0, 0.0, 0.0, -5.0))
end

### SYSTEMS

#! format: off

cooketriplet(::Type{T} = Float64) where {T<:Real} = AxisymmetricOpticalSystem{T}(
    DataFrame(SurfaceType = ["Object", "Standard", "Standard", "Standard", "Stop", "Standard", "Standard", "Image"],
              Radius = [Inf, 26.777, 66.604, -35.571, 35.571, 35.571, -26.777, Inf],
              Thickness = [Inf, 4.0, 2.0, 4.0, 2.0, 4.0, 44.748, missing],
              Material = [OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_SK16, OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_SF2, OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_SK16, OpticSim.GlassCat.Air, missing],
              SemiDiameter = [Inf, 8.580, 7.513, 7.054, 6.033, 7.003, 7.506, 15.0]))

cooketripletfirstelement(::Type{T} = Float64) where {T<:Real} = AxisymmetricOpticalSystem(
    DataFrame(SurfaceType = ["Object", "Standard", "Standard", "Image"],
              Radius = [Inf, -35.571, 35.571, Inf],
              Thickness = [Inf, 4.0, 44.748, missing],
              Material = [OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_SK16, OpticSim.GlassCat.Air, missing],
              SemiDiameter = [Inf, 7.054, 6.033, 15.0]))

convexplano(::Type{T} = Float64) where {T<:Real} = AxisymmetricOpticalSystem{T}(
    DataFrame(SurfaceType = ["Object", "Standard", "Standard", "Image"],
              Radius = [Inf, 60.0, Inf, Inf],
              Thickness = [Inf, 10.0, 57.8, missing],
              Material = [OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, missing],
              SemiDiameter = [Inf, 9.0, 9.0, 15.0]))

doubleconvex(frontradius::T,rearradius::T) where{T<:Real} =
AxisymmetricOpticalSystem{T}(
    DataFrame(SurfaceType = ["Object", "Standard", "Standard", "Image"],
              Radius = [T(Inf64), frontradius, rearradius, T(Inf64)],
              Thickness = [T(Inf64), T(10.0), T(57.8), missing],
              Material = [OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, missing],
              SemiDiameter = [T(Inf64), T(9.0), T(9.0), T(15.0)]),
              temperature = OpticSim.GlassCat.TEMP_REF_UNITFUL, pressure = T(OpticSim.GlassCat.PRESSURE_REF))

doubleconvex(::Type{T} = Float64; temperature::Unitful.Temperature = OpticSim.GlassCat.TEMP_REF_UNITFUL, pressure::Float64 = OpticSim.GlassCat.PRESSURE_REF) where {T<:Real} =
AxisymmetricOpticalSystem{T}(
    DataFrame(SurfaceType = ["Object", "Standard", "Standard", "Image"],
              Radius = [Inf64, 60, -60, Inf64],
              Thickness = [Inf64, 10.0, 57.8, missing],
              Material = [OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, missing],
              SemiDiameter = [Inf64, 9.0, 9.0, 15.0]),
              temperature = temperature, pressure = T(pressure))

doubleconcave(::Type{T} = Float64) where {T<:Real} = AxisymmetricOpticalSystem{T}(
    DataFrame(SurfaceType = ["Object", "Standard", "Standard", "Image"],
              Radius = [Inf64, -41.0, 41.0, Inf64],
              Thickness = [Inf64, 10.0, 57.8, missing],
              Material = [OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, missing],
              SemiDiameter = [Inf64, 9.0, 9.0, 15.0]))

planoconcaverefl(::Type{T} = Float64) where {T<:Real} = AxisymmetricOpticalSystem{T}(
    DataFrame(SurfaceType = ["Object", "Standard", "Standard", "Image"],
              Radius = [Inf64, Inf64, -41.0, Inf64],
              Thickness = [Inf64, 10.0, -57.8, missing],
              Material = [OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, missing],
              SemiDiameter = [Inf64, 9.0, 9.0, 25.0],
              Reflectance = [missing, missing, 1.0, missing]))

concaveplano(::Type{T} = Float64) where {T<:Real} = AxisymmetricOpticalSystem{T}(
    DataFrame(SurfaceType = ["Object", "Standard", "Standard", "Image"],
              Radius = [Inf64, -41.0, Inf64, Inf64],
              Thickness = [Inf64, 10.0, 57.8, missing],
              Material = [OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, missing],
              SemiDiameter = [Inf64, 9.0, 9.0, 15.0]))

planoplano(::Type{T} = Float64) where {T<:Real} = AxisymmetricOpticalSystem{T}(
    DataFrame(SurfaceType = ["Object", "Standard", "Standard", "Image"],
              Radius = [Inf64, Inf64, Inf64, Inf64],
              Thickness = [Inf64, 10.0, 57.8, missing],
              Material = [OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, missing],
              SemiDiameter = [Inf64, 9.0, 9.0, 15.0]))

#! format: on

zernikesystem() = CSGOpticalSystem(LensAssembly(zernikelens()), Rectangle(40.0, 40.0, [0.0, 0.0, 1.0], [0.0, 0.0, -67.8], interface = opaqueinterface()))
asphericsystem() = CSGOpticalSystem(LensAssembly(asphericlens()), Rectangle(25.0, 25.0, [0.0, 0.0, 1.0], [0.0, 0.0, -67.8], interface = opaqueinterface()))
conicsystemZ() = CSGOpticalSystem(LensAssembly(coniclensZ()), Rectangle(25.0, 25.0, [0.0, 0.0, 1.0], [0.0, 0.0, -67.8], interface = opaqueinterface()))
conicsystemQ() = CSGOpticalSystem(LensAssembly(coniclensQ()), Rectangle(25.0, 25.0, [0.0, 0.0, 1.0], [0.0, 0.0, -67.8], interface = opaqueinterface()))
qtypesystem() = CSGOpticalSystem(LensAssembly(qtypelens()), Rectangle(25.0, 25.0, [0.0, 0.0, 1.0], [0.0, 0.0, -67.8], interface = opaqueinterface()))
chebyshevsystem() = CSGOpticalSystem(LensAssembly(chebyshevlens()), Rectangle(40.0, 40.0, [0.0, 0.0, 1.0], [0.0, 0.0, -67.8], interface = opaqueinterface()))
gridsagsystem() = CSGOpticalSystem(LensAssembly(gridsaglens()), Rectangle(40.0, 40.0, [0.0, 0.0, 1.0], [0.0, 0.0, -67.8], interface = opaqueinterface()))

### OTHER

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
    int1 = HologramInterface(SVector(-5.0, 0.0, -20.0), ConvergingBeam, SVector(0.0, -1.0, -1.0), CollimatedBeam, 0.55, 100.0, OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, OpticSim.GlassCat.Air, OpticSim.GlassCat.Air, 0.05, false)
    int2 = HologramInterface(SVector(5.0, 0.0, -20.0), ConvergingBeam, SVector(0.0, 1.0, -1.0), CollimatedBeam, 0.55, 100.0, OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, OpticSim.GlassCat.Air, OpticSim.GlassCat.Air, 0.05, false)
    mint = MultiHologramInterface(int1, int2)
    obj = MultiHologramSurface(rect, mint)
    sys = CSGOpticalSystem(LensAssembly(obj), Rectangle(50.0, 50.0, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, -20.0), interface = opaqueinterface()))
    return sys
end

end # module TestData
export TestData
