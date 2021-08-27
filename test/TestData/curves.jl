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
