# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.


""" Computes the location of the optical center of the lens that will project the centroid of the display to the centroid of the eyebox. Normally the display centroid will be aligned with the geometric centroid of the lens, rather than the optical center of the lens."""
function compute_optical_center(eyeboxcentroid,displaycenter,lensplane)
    v = displaycenter - eyeboxcentroid
    r = Ray(eyeboxcentroid,v)
    return closestintersection(surfaceintersection(lensplane))
end

function subdivide_fov(eyeboxpolygon)
end

""" eyebox rectangle represented as a 3x4 SMatrix with x,y,z coordinates in rows 1,2,3 respectively."""
function eyeboxpolygon(xsize,ysize)
    xext,yext = (xsize,ysize) ./ 2.0 #divide by 2 to get min and max extensions
    return SMatrix{3,4}([xext;yext;0.0 ;; xext;-yext;0.0 ;; -xext;-yext;0.0 ;; -xext;yext;0.0])
end

subpoints(pts::AbstractVector,subdivisions) = reshape(collect(range(extrema(pts)...,subdivisions)),1,subdivisions)
export subpoints

subdivide(poly::AbstractMatrix,xsubdivisions,ysubdivisions) = subdivide(SMatrix{3,4}(poly),xsubdivisions,ysubdivisions)

"""poly is a two dimensional rectangle represented as a 3x4 SMatrix with x values in row 1 and y values in row 2. Returns a vector of length xsubdivisions*ysubdivisions containing all the subdivided rectangles."""
function subdivide(poly::T, xsubdivisions, ysubdivisions) where{T<:SMatrix{3,4}}
    result = Vector{T}(undef,0) #efficiency not an issue here so don't need to preallocate vector
    
    x = subpoints(poly[1,:],xsubdivisions+1) #range function makes one less subdivision that people will be expecting, e.g., collect(range(1,3,2)) returns [1,3] 
    y = subpoints(poly[2,:],ysubdivisions+1)

    for i in 1:length(x) -1
        for j in 1:length(y) -1
            push!(result,T([x[i];y[j];0.0 ;; x[i+1];y[j];0.0 ;; x[i+1];y[j+1];0.0 ;; x[i];y[j+1];0.0]))
        end
    end
    return result
end
export subdivide
    
function setup_system(subdivisions::NTuple{2,Int})
    eyeballframe = Transform(I(4)) #This establishes the global coordinate frame for the systems. Positive Z axis is assumed to be the forward direction, the default direction of the eye's optical axis when looking directly ahead.
    corneavertex = OpticSim.HumanEye.cornea_to_eyecenter()
    

    eyeboxframe = eyeballframe*OpticSim.translation(0.0,0.0,ustrip(mm,corneavertex))  #unfortunately can't use unitful values in transforms because the rotation and translation components would have different types which is not allowed in a Matrix.

    println(typeof(eyeboxframe))
    eyeboxpoly =  eyeboxpolygon(10mm,9mm) #four corners of the eyebox frame which is assumed centered around the positive Z axis. Transformed to the eyeballframe.
    println(typeof(subdivide(eyeboxpoly,subdivisions...)[1]))
    subdivided_eyeboxpolys = map(x-> eyeboxframe .* x, subdivide(eyeboxpoly,subdivisions...))  
end
export setup_system
