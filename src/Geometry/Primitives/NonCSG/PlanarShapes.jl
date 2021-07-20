
"""The PlanarShapes interface:

`distancefromplane(p::PlanarShapes,point)`  returns distance of the point from the plane the planar shape lies within
"""
abstract type PlanarShapes{T} <: Surface{T} end

"""All planar shapes lie on a plane. This function computes the distance from a point to that plane. This is a signed distance. If the point is on the positive side of the plane (the side the normal points toward) the distance will be positive, otherwise negative or 0 if the point lies in the plane."""
distancefromplane(p::PlanarShapes, point::SVector{3}) = distancefromplane(p.plane,point)