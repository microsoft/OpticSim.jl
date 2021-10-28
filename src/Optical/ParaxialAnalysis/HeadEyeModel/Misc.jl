# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

# based on code from https://www.ilikebigbits.com/2015_03_04_plane_from_points.html
"""
    plane_from_points(points::Vector{SVector{3, Float64}}) ->  centroid, normal

Estimate the best fitting plane for a set of points in 3D.
    
See utilities.jl for a simpler and possibly more numerically stable version.  
"""
function plane_from_points2(points::Vector{SVector{3, Float64}}) 
    if length(points) < 3 
        return nothing      # At least three points required
    end

    sum = Base.sum(points)
    centroid = Vec3(sum * (1.0 / (length(points))))

    # Calc full 3x3 covariance matrix, excluding symmetries:
    xx = 0.0 
    xy = 0.0
    xz = 0.0
    yy = 0.0 
    yz = 0.0
    zz = 0.0

    for p in points 
        r = p - centroid;
        xx += r[1] * r[1]
        xy += r[1] * r[2]
        xz += r[1] * r[3]
        yy += r[2] * r[2]
        yz += r[2] * r[3]
        zz += r[3] * r[3]
    end

    det_x = yy*zz - yz*yz
    det_y = xx*zz - xz*xz
    det_z = xx*yy - xy*xy

    det_max = max(det_x, det_y, det_z)
    if det_max <= 0.0 
        return nothing      # The points don't span a plane
    end

    # Pick path with best conditioning:
    if      (det_max == det_x)
        dir = Vec3(det_x,xz*yz - xy*zz,xy*yz - xz*yy)
    elseif  (det_max == det_y)
        dir = Vec3(xz*yz - xy*zz, det_y, xy*xz - yz*xx)
    else
        dir = Vec3(xy*yz - xz*yy,xy*xz - yz*xx,det_z)
    end
    
    return centroid, normalize(dir) 
end


function csg_sphere(;radius = 10.0)
    sph = Sphere(radius)
    csg = leaf(sph, Transform(Vec3(0.0, 0.0, radius)))
    csg = csg(identitytransform())
    return csg
end

function csg_cylinder(;radius = 10.0, height = 50.0, added_rotation = rotationY(π/2.0))
    cyl = Cylinder(radius, height)
    csg = leaf(cyl, Transform(Vec3(0.0, 0.0, radius)) * added_rotation)
    csg = csg(identitytransform())
    return csg
end

function csg_plane()
    plane = Plane(unitZ3() * -1, zero(Vec3))
    csg = leaf(plane, identitytransform())
    csg = csg(identitytransform())
    return csg
end
