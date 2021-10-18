# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

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
