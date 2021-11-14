# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.


"""The AbstractPlanarShape interface:

`distancefromplane(p::AbstractPlanarShape,point)`  returns distance of the point from the plane the planar shape lies within
`normal(p::AbstractPlanarShape)` returns normal of plane
`interface(p::AbstractPlanarShape)` returns optical interface of plane
`vertices(p::AbstractPlanarShape)` returns vertices of shape. For Ellipse this is an approximation.

There are default functions for plane,normal,interface,vertices which assume each AbstractPlanarShape type has a field of the same name
`  
    plane(a::AbstractPlanarShape) = a.plane
    normal(a::PlanaShape) = a.plane.normal
`
etc.

If your type doesn't have these fields then you should define a more specialized method to handle this.
"""
abstract type AbstractPlanarShape{T} <: Surface{T} end

"""All planar shapes lie on a plane. This function computes the distance from a point to that plane. This is a signed distance. If the point is on the positive side of the plane (the side the normal points toward) the distance will be positive, otherwise negative or 0 if the point lies in the plane."""
distancefromplane(p::AbstractPlanarShape, point::SVector{3}) = distancefromplane(p.plane,point)
normal(p::AbstractPlanarShape) = normal(p.plane)
interface(p::AbstractPlanarShape) = interface(p.plane)
"""The vertices of planar shapes are defined in a plane so they are two dimensional. In the local coordinate frame this is the x,y plane, so the implied z coordinate is 0"""
vertices(p::AbstractPlanarShape) = throw(ErrorException("This function should be defined for any concrete type that is a subtype of AbstractPlanarShape")) #don't think this function can ever be called because you can't instantiate something of PlanarShape type; it's an abstract type.
export vertices

