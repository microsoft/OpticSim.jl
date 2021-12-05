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
            push!(result,T([x[i];y[j];0.0mm ;; x[i+1];y[j];0.0mm ;; x[i+1];y[j+1];0.0mm ;; x[i];y[j+1];0.0mm])) #this is klunky but the polygon vertices have units of mm so the 0.0 terms must as well. 
        end
    end
    return result
end
export subdivide
    
function setup_system()
    centroid(verts) = sum(eachcol(verts)) ./ size(verts)[2] 

    eyeballframe = Transform() #This establishes the global coordinate frame for the systems. Positive Z axis is assumed to be the forward direction, the default direction of the eye's optical axis when looking directly ahead.
    corneavertex = OpticSim.HumanEye.cornea_to_eyecenter()

    props = systemproperties(18mm,(10mm,9mm),(55°,45°),4.0mm,.2,11,30)
   
    subdivisions = props[:subdivisions]
    clusterdata = props[:cluster_data]
    cluster = clusterdata[:cluster]

    eyeboxframe = eyeballframe*OpticSim.translation(0.0,0.0,ustrip(mm,corneavertex))  #unfortunately can't use unitful values in transforms because the rotation and translation components would have different types which is not allowed in a Matrix.

    eyeboxpoly =  eyeboxpolygon(10mm,9mm) #four corners of the eyebox frame which is assumed centered around the positive Z axis. Transformed to the eyeballframe.

    polys = subdivide(eyeboxpoly,subdivisions...)
    strippedpolys =  map(x-> ustrip.(mm,x), polys)
    subdivided_eyeboxpolys =[eyeboxframe * x for x in strippedpolys] #these have units of mm which don't interact well with Transform.
    polycentroids = centroid.(subdivided_eyeboxpolys)

    focallength = ustrip(mm,props[:focal_length])
    lenses,coordinates = spherelenslets(Plane(0.0,0.0,1.0,0.0,0.0,12.0),focallength,[0.0,0.0,-1.0],30.0,55°,45°,HexBasis1())
    lensletcolors = [pointcolor(coordinates[:,j],cluster) for j in 1:size(coordinates)[2]]
    
    #compute offset of optical axis for lenslets so they cover the correct part of the eyebox.



    #project eyebox into lenslet display plane and compute bounding box. This is the size of the display for this lenslet

    return (lenses,coordinates,lensletcolors)
end
export setup_system
