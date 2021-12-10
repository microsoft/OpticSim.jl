# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

import Statistics

""" Computes the location of the optical center of the lens that will project the centroid of the display to the centroid of the eyebox. Normally the display centroid will be aligned with the geometric centroid of the lens, rather than the optical center of the lens."""
function compute_optical_center(eyeboxcentroid,displaycenter,lensplane)
    v = displaycenter - eyeboxcentroid
    r = Ray(eyeboxcentroid,v)
    intsct = point(closestintersection(surfaceintersection(lensplane,r),false)) #this seems complicated when you don't need CSG
    @assert intsct !== nothing
    return intsct
end

replace_optical_center(lens,new_optical_center) = ParaxialLensConvexPoly(focallength(lens),shape(lens),new_optical_center)


function replace_optical_center(eyeboxcentroid,displaycenter,lens)
    transform(lens)
     replace_optical_center(lens,compute_optical_center(eyeboxcentroid,displaycenter,OpticSim.plane(lens)))
end


"""returns display plane represented in world coordinates, and the center point of the display"""
function display_plane(lens) 
    center_point = centroid(lens) + -normal(lens)* OpticSim.focallength(lens)
     pln = Plane(-normal(lens), center_point, vishalfsizeu = 5.0, vishalfsizev = 5.0)
    return pln,center_point
end
export display_plane


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

function eyebox_assignment(tilecoords::NTuple{2,Int64},cluster::R,eyeboxes::AbstractVector) where{R<:AbstractLatticeCluster}
    #compute cluster
    _, tile_index = cluster_coordinates_from_tile_coordinates(cluster,tilecoords)
   return eyeboxes[mod(tile_index-1,length(eyeboxes)) + 1] #use linear index for Matrix.
end

"""given the eye box polygon and how it is to be subdivided computes the subdivided eyeboxes and centroids"""
function compute_lenslet_eyebox_data(eyeboxtransform,eyeboxpoly,subdivisions)
    subdivided_eyeboxpolys = subdivide(eyeboxpoly,subdivisions...)
    strippedpolys =  map(x-> ustrip.(mm,x),  subdivided_eyeboxpolys)
    subdivided_eyeboxpolys =[eyeboxtransform * x for x in strippedpolys] #these have units of mm which don't interact well with Transform.
    polycentroids = Statistics.mean.(subdivided_eyeboxpolys)
    

    return subdivided_eyeboxpolys,polycentroids
end

function setup_system()
    centroid(verts) = sum(eachcol(verts)) ./ size(verts)[2] 

    eyeballframe = Transform() #This establishes the global coordinate frame for the systems. Positive Z axis is assumed to be the forward direction, the default direction of the eye's optical axis when looking directly ahead.
    corneavertex = OpticSim.Data.cornea_to_eyecenter()
    defaulteyebox = (10mm,9mm)

    fov = (90°,60°)
    eyeboxsize = (10mm,8mm)
    eyerelief = 20mm
    pupildiameter = 3.5mm
    displaysphereradius = 125.0mm

    props = systemproperties(eyerelief,eyeboxsize,fov,pupildiameter,.22,11,pixelpitch = .9μm, minfnumber = 2.0)
    subdivisions = props[:subdivisions]
    clusterdata = props[:cluster_data]
    cluster = clusterdata[:cluster]

    eyeboxtransform = eyeballframe*OpticSim.translation(0.0,0.0,ustrip(mm,corneavertex))  #unfortunately can't use unitful values in transforms because the rotation and translation components would have different types which is not allowed in a Matrix.

    eyeboxpoly =  eyeboxpolygon(defaulteyebox...) #four corners of the eyebox frame which is assumed centered around the positive Z axis. Transformed to the eyeballframe.



    focallength = ustrip(mm,props[:focal_length])

    lenses,coordinates = spherelenslets(Plane(0.0,0.0,1.0,0.0,0.0,12.0),eyerelief,focallength,[0.0,0.0,-1.0],displaysphereradius,fov[1],fov[2],HexBasis1())


    temp = display_plane.(lenses)
    displayplanes = [x[1] for x in temp]
    planecenters = [x[2] for x in temp]

    lensletcolors = pointcolor.(coordinates,Ref(cluster))

    subdivided_eyeboxpolys,polycentroids = compute_lenslet_eyebox_data(eyeboxtransform,eyeboxpoly,subdivisions)
    lensleteyeboxes = eyebox_assignment.(coordinates,Ref(cluster),Ref(subdivided_eyeboxpolys))


    lensleteyeboxcenters = [Statistics.mean(eachcol(x)) for x in lensleteyeboxes]

    #make new lenses with optical centers that will cause the centroid of the eyebox assigned to the lens to project to the center of the display plane.
    lenses = replace_optical_center.(lensleteyeboxcenters,planecenters,lenses)

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

function test_eyebox_assignment()
    eyebox_assignment([(0,0),(0,1)],hex19(),)
end
