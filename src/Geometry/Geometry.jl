# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

# declaring the Geometry module. in the future we might add or move existing code to it. currently it contains the basic types: Vec3, Vec4 and Transform
module Geometry

include("Transform.jl")

end # module Geometry
export Geometry


using .Geometry

include("Ray.jl")

include("Surface.jl")

# not ideal having this like this, but we need the types for Intersection
include("../Optical/OpticalInterface.jl")

include("CSG/Intersection.jl")
include("CSG/Interval.jl")

include("Primitives/NonCSG/Triangle.jl")

include("BoundingBox.jl")

include("Primitives/Plane.jl")
include("Primitives/Cylinder.jl")
include("Primitives/Sphere.jl")
include("Primitives/SphericalCap.jl")

include("Primitives/NonCSG/Rectangle.jl")
include("Primitives/NonCSG/Hexagon.jl")
include("Primitives/NonCSG/Ellipse.jl")
include("Primitives/NonCSG/Stop.jl")

# include("BoundingVolumeHierarchy.jl")

include("AccelSurface.jl")

include("Primitives/Curves/Knots.jl")
include("Primitives/Curves/Spline.jl")
include("Primitives/Curves/BSpline.jl")
include("Primitives/Curves/Bezier.jl")
include("Primitives/Curves/PowerBasis.jl")

include("CSG/CSG.jl")

include("Primitives/Zernike.jl")
include("Primitives/Qtype.jl")
include("Primitives/Chebyshev.jl")
include("Primitives/GridSag.jl")

include("AnalyticIntersection.jl")

###########################################################################################################

"""
    BoundedCylinder(radius::T, height::T; interface::NullOrFresnel{T} = nullinterface(T)) -> CSGGenerator{T}

Create a cylinder with planar caps on both ends centred at `(0, 0, 0)` with axis `(0, 0, 1)`.
"""
function BoundedCylinder(radius::T, height::T; interface::NullOrFresnel{T} = NullInterface(T)) where {T<:Real}
    barrel = Cylinder(radius, height, interface = interface)
    top = Plane(SVector{3,T}(0.0, 0.0, 1.0), SVector{3,T}(0.0, 0.0, height / 2), vishalfsizeu = radius, vishalfsizev = radius)
    bottom = Plane(SVector{3,T}(0.0, 0.0, -1.0), SVector{3,T}(0.0, 0.0, -height / 2), vishalfsizeu = radius, vishalfsizev = radius)
    return csgintersection(leaf(barrel), csgintersection(top, bottom))
end
export BoundedCylinder

"""
    Cuboid(halfsizex::T, halfsizey::T, halfsizez::T; interface::NullOrFresnel{T} = nullinterface(T)) -> CSGGenerator{T}

Create a cuboid centred at `(0, 0, 0)`.
"""
function Cuboid(halfsizex::T, halfsizey::T, halfsizez::T; interface::NullOrFresnel{T} = NullInterface(T)) where {T<:Real}
    xmin = Plane(SVector{3,T}(-1.0, 0.0, 0.0), SVector{3,T}(-halfsizex, 0.0, 0.0), vishalfsizeu = halfsizez, vishalfsizev = halfsizey, interface = interface)
    xmax = Plane(SVector{3,T}(1.0, 0.0, 0.0), SVector{3,T}(halfsizex, 0.0, 0.0), vishalfsizeu = halfsizez, vishalfsizev = halfsizey, interface = interface)
    ymin = Plane(SVector{3,T}(0.0, -1.0, 0.0), SVector{3,T}(0.0, -halfsizey, 0.0), vishalfsizeu = halfsizez, vishalfsizev = halfsizex, interface = interface)
    ymax = Plane(SVector{3,T}(0.0, 1.0, 0.0), SVector{3,T}(0.0, halfsizey, 0.0), vishalfsizeu = halfsizez, vishalfsizev = halfsizex, interface = interface)
    zmin = Plane(SVector{3,T}(0.0, 0.0, -1.0), SVector{3,T}(0.0, 0.0, -halfsizez), vishalfsizeu = halfsizex, vishalfsizev = halfsizey, interface = interface)
    zmax = Plane(SVector{3,T}(0.0, 0.0, 1.0), SVector{3,T}(0.0, 0.0, halfsizez), vishalfsizeu = halfsizex, vishalfsizev = halfsizey, interface = interface)
    xint = csgintersection(xmin, xmax)
    yint = csgintersection(ymin, ymax)
    zint = csgintersection(zmin, zmax)
    return csgintersection(xint, csgintersection(yint, zint))
end
export Cuboid

"""
    HexagonalPrism(side_length::T, visheight::T = 2.0; interface::NullOrFresnel{T} = nullinterface(T)) -> CSGGenerator{T}

Create an infinitely tall hexagonal prism with axis `(0, 0, 1)`, the longer hexagon diameter is along the x axis.
For visualization `visheight` is used, **note that this does not fully represent the surface**.
"""
function HexagonalPrism(side_length::T, visheight::T = 2.0; interface::NullOrFresnel{T} = NullInterface(T)) where {T<:Real}
    h = side_length * sqrt(3.0) / 2.0
    hcos = h * sqrt(3.0) / 2.0
    hsin = h / 2.0
    side_length *= 1.0001
    s1 = Plane(SVector{3,T}(0.0, 1.0, 0.0), SVector{3,T}(0.0, h, 0.0), vishalfsizeu = visheight / 2.0, vishalfsizev = side_length / 2.0, interface = interface)
    s2 = Plane(SVector{3,T}(0.0, -1.0, 0.0), SVector{3,T}(0.0, -h, 0.0), vishalfsizeu = visheight / 2.0, vishalfsizev = side_length / 2.0, interface = interface)
    s3 = Plane(SVector{3,T}(hcos, hsin, 0.0), SVector{3,T}(hcos, hsin, 0.0), vishalfsizeu = visheight / 2.0, vishalfsizev = side_length / 2.0, interface = interface)
    s4 = Plane(SVector{3,T}(-hcos, -hsin, 0.0), SVector{3,T}(-hcos, -hsin, 0.0), vishalfsizeu = visheight / 2.0, vishalfsizev = side_length / 2.0, interface = interface)
    s5 = Plane(SVector{3,T}(hcos, -hsin, 0.0), SVector{3,T}(hcos, -hsin, 0.0), vishalfsizeu = visheight / 2.0, vishalfsizev = side_length / 2.0, interface = interface)
    s6 = Plane(SVector{3,T}(-hcos, hsin, 0.0), SVector{3,T}(-hcos, hsin, 0.0), vishalfsizeu = visheight / 2.0, vishalfsizev = side_length / 2.0, interface = interface)
    i1 = csgintersection(s1, s2)
    i2 = csgintersection(s3, s4)
    i3 = csgintersection(s5, s6)
    return csgintersection(i1, csgintersection(i2, i3))
end
export HexagonalPrism

"""
    RectangularPrism(halfsizex::T, halfsizey::T, visheight::T=2.0; interface::NullOrFresnel{T} = nullinterface(T)) -> CSGGenerator{T}

Create an infinitely tall rectangular prism with axis `(0, 0, 1)`.
For visualization `visheight` is used, **note that this does not fully represent the surface**.
"""
function RectangularPrism(halfsizex::T, halfsizey::T, visheight::T = 2.0; interface::NullOrFresnel{T} = NullInterface(T)) where {T<:Real}
    s1 = Plane(SVector{3,T}(0.0, 1.0, 0.0), SVector{3,T}(0.0, halfsizey, 0.0), vishalfsizeu = visheight / 2.0, vishalfsizev = halfsizex, interface = interface)
    s2 = Plane(SVector{3,T}(0.0, -1.0, 0.0), SVector{3,T}(0.0, -halfsizey, 0.0), vishalfsizeu = visheight / 2.0, vishalfsizev = halfsizex, interface = interface)
    s3 = Plane(SVector{3,T}(1.0, 0.0, 0.0), SVector{3,T}(halfsizex, 0.0, 0.0), vishalfsizeu = visheight / 2.0, vishalfsizev = halfsizey, interface = interface)
    s4 = Plane(SVector{3,T}(-1.0, 0.0, 0.0), SVector{3,T}(-halfsizex, 0.0, 0.0), vishalfsizeu = visheight / 2.0, vishalfsizev = halfsizey, interface = interface)
    i1 = csgintersection(s1, s2)
    i2 = csgintersection(s3, s4)
    return csgintersection(i1, i2)
end
export RectangularPrism

"""
    TriangularPrism(side_length::T, visheight::T = 2.0; interface::NullOrFresnel{T} = nullinterface(T)) -> CSGGenerator{T}

Create an infinitely tall triangular prism with axis `(0, 0, 1)`.
For visualization `visheight` is used, **note that this does not fully represent the surface**.
"""
function TriangularPrism(side_length::T, visheight::T = 2.0; interface::NullOrFresnel{T} = NullInterface(T)) where {T<:Real}
    side_length = side_length / 2
    p1 = Plane(SVector{3,T}(-1.0, 0.0, 0.0), SVector{3,T}(-side_length / 2, 0.0, 0.0), interface = interface, vishalfsizev = side_length, vishalfsizeu = visheight / 2)
    p2 = Plane(SVector{3,T}(sind(30), cosd(30), 0.0), SVector{3,T}(side_length / 2 * sind(30), side_length / 2 * cosd(30), 0.0), interface = interface, vishalfsizev = side_length, vishalfsizeu = visheight / 2)
    p3 = Plane(SVector{3,T}(sind(30), -cosd(30), 0.0), SVector{3,T}(side_length / 2 * sind(30), -side_length / 2 * cosd(30), 0.0), interface = interface, vishalfsizev = side_length, vishalfsizeu = visheight / 2)
    return csgintersection(leaf(p1), csgintersection(p2, p3))
end
export TriangularPrism

"""
    Spider(narms::Int, armwidth::T, radius::T, origin::SVector{3,T} = SVector{3,T}(0.0, 0.0, 0.0), normal::SVector{3,T} = SVector{3,T}(0.0, 0.0, 1.0)) -> Vector{Rectangle{T}}

Creates a 'spider' obscuration with `narms` rectangular arms evenly spaced around a circle defined by `origin` and `normal`.
Each arm is a rectangle `armwidth`×`radius`.

e.g. for 3 and 4 arms we get:
```
   |         _|_
  / \\         |
```
"""
function Spider(narms::Int, armwidth::T, radius::T, origin::SVector{3,T} = SVector{3,T}(0.0, 0.0, 0.0), normal::SVector{3,T} = SVector{3,T}(0.0, 0.0, 1.0)) where {T<:Real}
    dθ = 2π / narms
    rects = Vector{Rectangle{T}}(undef, 0)
    for i in 0:(narms - 1)
        θ = i * dθ
        r = Rectangle(armwidth / 2, radius / 2, normal, origin + SVector(radius / 2 * cos(θ), radius / 2 * sin(θ), 0.0), rotationvec = SVector(cos(θ), sin(θ), 0.0), interface = opaqueinterface(T))
        push!(rects, r)
    end
    return rects
end
export Spider
