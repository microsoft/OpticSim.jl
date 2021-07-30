
"""The PlanarShapes interface:

`distancefromplane(p::PlanarShapes,point)`  returns distance of the point from the plane the planar shape lies within
`normal(p::PlanarShapes)` returns normal of plane
`interface(p::PlanarShapes)` returns optical interface of plane
`vertices(p::PlanarShapes)` returns vertices of shape. For Ellipse this is an approximation
"""
abstract type PlanarShapes{T} <: Surface{T} end

"""All planar shapes lie on a plane. This function computes the distance from a point to that plane. This is a signed distance. If the point is on the positive side of the plane (the side the normal points toward) the distance will be positive, otherwise negative or 0 if the point lies in the plane."""
distancefromplane(p::PlanarShapes, point::SVector{3}) = distancefromplane(p.plane,point)
normal(p::PlanarShapes) = normal(p.plane)
interface(p::PlanarShapes) = interface(p.plane)
"""The vertices of planar shapes are defined in a plane so they are two dimensional. In the local coordinate frame this is the x,y plane, so the implied z coordinate is 0"""
vertices(p::PlanarShapes) = throw(ErrorException("This function should be defined for any concrete type that is a subtype of PlanarShapes")) #don't think this function can ever be called because you can't instantiate something of PlanarShapes type; it's an abstract type.
export vertices

"""expand 2D vertices defined in the plane to 3D vertices defined with z = 0, i.e., the location of the plane in the local shape coordinate frame"""
vertices3d(a::SMatrix{2,N,T}) where{N,T<:Real} = vcat(a,zeros(SMatrix{1,N,T}))
export vertices3d
