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

replace_optical_center(lens,new_optical_center) = ParaxialLensConvexPoly(focallength(lens),OpticSim.shape(lens),new_optical_center)

function localframe(a::ParaxialLens)
    if OpticSim.shape(a) isa ConvexPolygon
        return OpticSim.localframe(OpticSim.shape(a))
    else
        return nothing
    end
end

function replace_optical_center(eyeboxcentroid,displaycenter,lens)
    world_optic_center = compute_optical_center(eyeboxcentroid,displaycenter,OpticSim.plane(lens))
    local_optic_center = SVector{2}((inv(localframe(lens)) * world_optic_center)[1:2]) #this is defined in the 2D lens plane so z = 0
    replace_optical_center(lens,local_optic_center)
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

"""System parameters for a typical HMD."""
function systemparameters() 
    return (eye_box = (10.0mm,9.0mm),fov = (90.0°,60.0°),eye_relief = 20.0mm, pupil_diameter = 3.5mm, display_sphere_radius = 125.0mm,min_fnumber = 2.0)
end

"""Coordinate frames for the eye/display system. The origin of this frame is at the geometric center of the eyeball. Positive Z axis is assumed to be the forward direction, the default direction of the eye's optical axis when looking directly ahead."""
function setup_coordinate_frames()
    eyeballframe = Transform()
    corneavertex = OpticSim.Data.cornea_to_eyecenter()
    eyeboxtransform = eyeballframe*OpticSim.translation(0.0,0.0,ustrip(mm,corneavertex))  #unfortunately can't use unitful values in transforms because the rotation and translation components would have different types which is not allowed in a Matrix.

    return (eyeball_frame = eyeballframe,eye_box_frame = eyeboxtransform)
end

"""Create lenslet system that will cover the eyebox and fov"""
function setup_system()
    #All coordinates are ultimately transformed into the eyeball_frame coordinate systems

    (eyeball_frame,eye_box_frame) = setup_coordinate_frames()
    (eye_box,fov,eye_relief,pupil_diameter,display_sphere_radius,min_fnumber) = systemparameters()
    eyeboxz = (eye_box_frame*SVector(0.0,0.0,0.0))[3]
    fov = (90°,60°)

    #get system properties
    props = systemproperties(eye_relief,eye_box,fov,pupil_diameter,.22,11,pixelpitch = .9μm, minfnumber = min_fnumber)
    subdivisions = props[:subdivisions] #tuple representing how the eyebox can be subdivided given the cluster used for the lenslets
    clusterdata = props[:cluster_data] 
    cluster = clusterdata[:cluster] #cluster that is repeated across the display to ensure continuous coverage of the eyebox and fov.
    focallength = ustrip(mm,props[:focal_length]) #strip units off because these don't work well with Transform

    #compute lenslets based on system properties. lattice_coordinates are the (i,j) integer lattice coordinates of the hexagonal lattice making up the display. These coordinates are used to properly assign color and subdivided eyebox to the lenslets.
    lenses,lattice_coordinates = spherelenslets(Plane(0.0,0.0,1.0,0.0,0.0,12.0),eye_relief,focallength,[0.0,0.0,-1.0],display_sphere_radius,fov[1],fov[2],HexBasis1())

    temp = display_plane.(lenses)
    displayplanes = [x[1] for x in temp]
    planecenters = [x[2] for x in temp]

    lensletcolors = pointcolor.(lattice_coordinates,Ref(cluster))

    #compute subdivided eyebox polygons and assign to appropriate lenslets
    eyeboxpoly  =  mm * (eye_box_frame * eyeboxpolygon(ustrip.(mm,eye_box)...)) #four corners of the eyebox frame which is assumed centered around the positive Z axis. Transformed to the eyeballframe. Have to switch back and forth between Unitful and unitless quantities because Transform doesn't work with Unitful values.
    subdivided_eyeboxpolys,polycentroids = compute_lenslet_eyebox_data(eye_box_frame,eyeboxpoly,subdivisions)
    lensleteyeboxes = eyebox_assignment.(lattice_coordinates,Ref(cluster),Ref(subdivided_eyeboxpolys))

    #compute display rectangles that will cover the assigned eyebox polygons
    lensleteyeboxcenters = [Statistics.mean(eachcol(x)) for x in lensleteyeboxes] 
    lenses = replace_optical_center.(lensleteyeboxcenters,planecenters,lenses)  #make new lenses with optical centers that will cause the centroid of the eyebox assigned to the lens to project to the center of the display plane.

    #project eyebox into lenslet display plane and compute bounding box. This is the size of the display for this lenslet
    eyeboxrect = Rectangle(ustrip(mm,eye_box[1]/2),ustrip(mm,eye_box[2]/2),[0.0,0.0,1.0],[0.0,0.0,eyeboxz])
    return (eyeboxrect,lenses,lattice_coordinates,lensletcolors)
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
   eyeboxrect, lenses,coords,colors = setup_system()
    Vis.draw() #clear screen
    Vis.draw!(eyeboxrect)
    for (lens,coord,color) in zip(lenses,coords,colors)
        Vis.draw!(lens,color = color)
    end
end
export test_eyebox_assignment
