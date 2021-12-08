# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

import Statistics

""" Computes the location of the optical center of the lens that will project the centroid of the display to the centroid of the eyebox. Normally the display centroid will be aligned with the geometric centroid of the lens, rather than the optical center of the lens."""
function compute_optical_center(eyeboxcentroid,displaycenter,lensplane)
    v = displaycenter - eyeboxcentroid
    r = Ray(eyeboxcentroid,v)
    return closestintersection(surfaceintersection(lensplane))
end

replace_optical_center(lens,new_optical_center) = ParaxialLensConvexPoly(lens.focallength,lens.convpoly,new_optical_center)

replace_optical_center(eyeboxcentroid,displaycenter,lens) = replace_optical_center(lens,compute_optical_center(eyeboxcentroid,displaycenter,surface(lens)))

function displayplane(lens) 
    center_point = centroid(lens) + normal(lens)
    pln = Plane(-normal(lens), center_point, vishalfsizeu = 5.0, vishalfsizev = 5.0)
    return pln,center_point
export displayplane


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

function eyebox_assignment(tilecoords::NTuple{2,Int64},cluster::R,eyeboxes::T) where{R<:AbstractLatticeCluster,T<:Vector}
    #compute cluster
    println(tilecoords)
    _, tile_index = cluster_coordinates_from_tile_coordinates(cluster,tilecoords)
   return eyeboxes[mod(tile_index-1,length(eyeboxes)) + 1] #use linear index for Matrix.
end




"""computes a vector from the centroid of the eyebox to the center of the display, which is the projection of the geometric center of the lens onto the display plane. This vector is intersected with the lens to compute an optical center which will 
"""
function project_eyebox_to_display(lens,eyebox)
    pln,center = displayplane(lens)
    boxctr = Statistics.mean(eyebox)
    vec = center-boxctr
    r = Ray(center,vec)
    optcenter = closestintersection(surfaceintersection(lens,r))

end

function setup_system()
    centroid(verts) = sum(eachcol(verts)) ./ size(verts)[2] 

    eyeballframe = Transform() #This establishes the global coordinate frame for the systems. Positive Z axis is assumed to be the forward direction, the default direction of the eye's optical axis when looking directly ahead.
    corneavertex = OpticSim.Data.cornea_to_eyecenter()

    props = systemproperties(18mm,(10mm,9mm),(55°,45°),4.0mm,.2,11)
   
    subdivisions = props[:subdivisions]
    clusterdata = props[:cluster_data]
    cluster = clusterdata[:cluster]

    eyeboxframe = eyeballframe*OpticSim.translation(0.0,0.0,ustrip(mm,corneavertex))  #unfortunately can't use unitful values in transforms because the rotation and translation components would have different types which is not allowed in a Matrix.

    eyeboxpoly =  eyeboxpolygon(10mm,9mm) #four corners of the eyebox frame which is assumed centered around the positive Z axis. Transformed to the eyeballframe.

    eyebox_subdivisions = subdivide(eyeboxpoly,subdivisions...)
    strippedpolys =  map(x-> ustrip.(mm,x), eyebox_subdivisions)
    subdivided_eyeboxpolys =[eyeboxframe * x for x in strippedpolys] #these have units of mm which don't interact well with Transform.
    polycentroids = centroid.(subdivided_eyeboxpolys)

    focallength = ustrip(mm,props[:focal_length])
    lenses,coordinates = spherelenslets(Plane(0.0,0.0,1.0,0.0,0.0,12.0),focallength,[0.0,0.0,-1.0],30.0,55°,45°,HexBasis1())
    println(coordinates)

    lensletcolors = pointcolor.(coordinates,Ref(cluster))
    
    lenslet_eyeboxes = eyebox_assignment.(coordinates,Ref(cluster),Ref(subdivided_eyeboxpolys))
    println(lenslet_eyeboxes)

    #project eyebox into lenslet display plane and compute bounding box. This is the size of the display for this lenslet

    return (lenses,coordinates,lensletcolors)
end
export setup_system

function testspherelenslets()
    (;fov, eye_relief,eyebox,display_radius,pupil_diameter,pixel_pitch,minfnumber,mtf,cycles_per_degree,max_display_size) = nominal_system_properties()
    
    computedprops = systemproperties(eye_relief,eyebox,fov,pupil_diameter,mtf,cycles_per_degree,minfnumber = minfnumber,maxdisplaysize = max_display_size,pixelpitch = pixel_pitch)
    focallength = ustrip(mm,computedprops[:focal_length])
    spherelenslets(Plane(0.0,0.0,1.0,0.0,0.0,18.0),focallength,[0.0,0.0,-1.0],ustrip(mm,display_radius),fov[1],fov[2],HexBasis1())
end
export testspherelenslets

