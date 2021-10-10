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

include("Primitives/NonCSG/PlanarShape.jl")
include("Primitives/NonCSG/Rectangle.jl")
include("Primitives/NonCSG/Hexagon.jl")
include("Primitives/NonCSG/ConvexPolygon.jl")
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

include("Primitives/AsphericSurface.jl")
include("Primitives/Zernike.jl")
include("Primitives/Qtype.jl")
include("Primitives/Chebyshev.jl")
include("Primitives/GridSag.jl")

include("AnalyticIntersection.jl")
include("SphericalPolygon.jl")

###########################################################################################################

export BoundedCylinder, Cuboid, HexagonalPrism, RectangularPrism, TriangularPrism, Spider

"""
    BoundedCylinder(radius::T, height::T; interface::NullOrFresnel{T} = nullinterface(T)) -> CSGGenerator{T}

Create a cylinder with planar caps on both ends centred at `(0, 0, 0)` with axis `(0, 0, 1)`.
"""
function BoundedCylinder(radius::T, height::T; interface::NullOrFresnel{T} = NullInterface(T)) where {T<:Real}
    barrel = Cylinder(radius, height, interface = interface)
    top = Plane(SVector{3,T}(0.0, 0.0, 1.0), SVector{3,T}(0.0, 0.0, height / 2), vishalfsizeu = radius, vishalfsizev = radius)
    bottom = Plane(SVector{3,T}(0.0, 0.0, -1.0), SVector{3,T}(0.0, 0.0, -height / 2), vishalfsizeu = radius, vishalfsizev = radius)
    return barrel ∩ top ∩ bottom
end

"""
    Cuboid(halfsizex::T, halfsizey::T, halfsizez::T; interface::NullOrFresnel{T} = nullinterface(T)) -> CSGGenerator{T}

Create a cuboid centred at `(0, 0, 0)`.
"""
function Cuboid(halfsizex::T, halfsizey::T, halfsizez::T; interface::NullOrFresnel{T} = NullInterface(T)) where {T<:Real}
    xmin = Plane(SVector{3,T}(-1.0, 0.0, 0.0), SVector{3,T}(-halfsizex, 0.0, 0.0); vishalfsizeu = halfsizez, vishalfsizev = halfsizey, interface)
    xmax = Plane(SVector{3,T}(1.0, 0.0, 0.0), SVector{3,T}(halfsizex, 0.0, 0.0); vishalfsizeu = halfsizez, vishalfsizev = halfsizey, interface)
    ymin = Plane(SVector{3,T}(0.0, -1.0, 0.0), SVector{3,T}(0.0, -halfsizey, 0.0); vishalfsizeu = halfsizez, vishalfsizev = halfsizex, interface)
    ymax = Plane(SVector{3,T}(0.0, 1.0, 0.0), SVector{3,T}(0.0, halfsizey, 0.0); vishalfsizeu = halfsizez, vishalfsizev = halfsizex, interface)
    zmin = Plane(SVector{3,T}(0.0, 0.0, -1.0), SVector{3,T}(0.0, 0.0, -halfsizez); vishalfsizeu = halfsizex, vishalfsizev = halfsizey, interface)
    zmax = Plane(SVector{3,T}(0.0, 0.0, 1.0), SVector{3,T}(0.0, 0.0, halfsizez); vishalfsizeu = halfsizex, vishalfsizev = halfsizey, interface)
    return xmin ∩ xmax ∩ ymin ∩ ymax ∩ zmin ∩ zmax
end

"""
    HexagonalPrism(side_length::T, visheight::T = 2.0; interface::NullOrFresnel{T} = nullinterface(T)) -> CSGGenerator{T}

Create an infinitely tall hexagonal prism with axis `(0, 0, 1)`, the longer hexagon diameter is along the x axis.
For visualization `visheight` is used, **note that this does not fully represent the surface**.
"""
function HexagonalPrism(side_length::T, visheight::T = 2.0; interface::NullOrFresnel{T} = NullInterface(T)) where {T<:Real}
    h = side_length * sqrt(3.0) / 2.0
    hcos = h * sqrt(3.0) / 2.0
    hsin = h / 2.0
    vishalfsizeu = visheight / 2.0
    vishalfsizev = side_length * 1.0001 / 2.0
    s1 = Plane(SVector{3,T}(0.0, 1.0, 0.0), SVector{3,T}(0.0, h, 0.0); vishalfsizeu, vishalfsizev, interface)
    s2 = Plane(SVector{3,T}(0.0, -1.0, 0.0), SVector{3,T}(0.0, -h, 0.0); vishalfsizeu, vishalfsizev, interface)
    s3 = Plane(SVector{3,T}(hcos, hsin, 0.0), SVector{3,T}(hcos, hsin, 0.0); vishalfsizeu, vishalfsizev, interface)
    s4 = Plane(SVector{3,T}(-hcos, -hsin, 0.0), SVector{3,T}(-hcos, -hsin, 0.0); vishalfsizeu, vishalfsizev, interface)
    s5 = Plane(SVector{3,T}(hcos, -hsin, 0.0), SVector{3,T}(hcos, -hsin, 0.0); vishalfsizeu, vishalfsizev, interface)
    s6 = Plane(SVector{3,T}(-hcos, hsin, 0.0), SVector{3,T}(-hcos, hsin, 0.0); vishalfsizeu, vishalfsizev, interface)
    return s1 ∩ s2 ∩ s3 ∩ s4 ∩ s5 ∩ s6
end

"""
    RectangularPrism(halfsizex::T, halfsizey::T, visheight::T=2.0; interface::NullOrFresnel{T} = nullinterface(T)) -> CSGGenerator{T}

Create an infinitely tall rectangular prism with axis `(0, 0, 1)`.
For visualization `visheight` is used, **note that this does not fully represent the surface**.
"""
function RectangularPrism(halfsizex::T, halfsizey::T, visheight::T = 2.0; interface::NullOrFresnel{T} = NullInterface(T)) where {T<:Real}
    vishalfsizeu = visheight / 2.0
    s1 = Plane(SVector{3,T}(0.0, 1.0, 0.0), SVector{3,T}(0.0, halfsizey, 0.0); vishalfsizeu, vishalfsizev=halfsizex, interface)
    s2 = Plane(SVector{3,T}(0.0, -1.0, 0.0), SVector{3,T}(0.0, -halfsizey, 0.0); vishalfsizeu, vishalfsizev=halfsizex, interface)
    s3 = Plane(SVector{3,T}(1.0, 0.0, 0.0), SVector{3,T}(halfsizex, 0.0, 0.0); vishalfsizeu, vishalfsizev=halfsizey, interface)
    s4 = Plane(SVector{3,T}(-1.0, 0.0, 0.0), SVector{3,T}(-halfsizex, 0.0, 0.0); vishalfsizeu, vishalfsizev=halfsizey, interface)
    return s1 ∩ s2 ∩ s3 ∩ s4
end

"""
    TriangularPrism(side_length::T, visheight::T = 2.0; interface::NullOrFresnel{T} = nullinterface(T)) -> CSGGenerator{T}

Create an infinitely tall triangular prism with axis `(0, 0, 1)`.
For visualization `visheight` is used, **note that this does not fully represent the surface**.
"""
function TriangularPrism(side_length::T, visheight::T = 2.0; interface::NullOrFresnel{T} = NullInterface(T)) where {T<:Real}
    vishalfsizeu = visheight / 2
    vishalfsizev = side_length / 2
    p1 = Plane(
        SVector{3,T}(-1.0, 0.0, 0.0),
        SVector{3,T}(-side_length / 4, 0.0, 0.0);
        interface, vishalfsizeu, vishalfsizev
    )
    p2 = Plane(
        SVector{3,T}(sind(30), cosd(30), 0.0),
        SVector{3,T}(side_length / 4 * sind(30), side_length / 4 * cosd(30), 0.0);
        interface, vishalfsizeu, vishalfsizev
    )
    p3 = Plane(
        SVector{3,T}(sind(30), -cosd(30), 0.0),
        SVector{3,T}(side_length / 4 * sind(30), -side_length / 4 * cosd(30), 0.0);
        interface, vishalfsizeu, vishalfsizev
    )
    return p1 ∩ p2 ∩ p3
end

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
