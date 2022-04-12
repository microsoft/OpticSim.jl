# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.


"""The PlanarShape interface:

`distancefromplane(p::PlanarShape,point)`  returns distance of the point from the plane the planar shape lies within
`normal(p::PlanarShape)` returns normal of plane
`interface(p::PlanarShape)` returns optical interface of plane
`vertices(p::PlanarShape)` returns 2D vertices of shape in the plane of the shape. For Ellipse this is an approximation.

There are default functions for plane,normal,interface,vertices which assume each PlanarShape type has a field of the same name
`  
    plane(a::PlanarShape) = a.plane
    normal(a::PlanaShape) = a.plane.normal
`
etc.

If your type doesn't have these fields then you should define a more specialized method to handle this.
"""
abstract type PlanarShape{T} <: Surface{T} end

"""All planar shapes lie on a plane. This function computes the distance from a point to that plane. This is a signed distance. If the point is on the positive side of the plane (the side the normal points toward) the distance will be positive, otherwise negative or 0 if the point lies in the plane."""
distancefromplane(p::PlanarShape, point::SVector{3}) = distancefromplane(p.plane,point)
normal(p::PlanarShape) = normal(p.plane)
interface(p::PlanarShape) = interface(p.plane)
"""The vertices of planar shapes are defined in a plane so they are two dimensional. In the local coordinate frame this is the x,y plane, so the implied z coordinate is 0"""
vertices(p::PlanarShape) = throw(ErrorException("This function should be defined for any concrete type that is a subtype of PlanarShape")) #don't think this function can ever be called because you can't instantiate something of PlanarShape type; it's an abstract type.
export vertices
plane(p::PlanarShape) = p.plane

