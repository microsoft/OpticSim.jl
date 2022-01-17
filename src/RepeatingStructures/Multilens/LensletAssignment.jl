# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

import Statistics

"""
    @unzip xs, ys, ... = us

will expand the assignment into the following code
    xs, ys, ... = map(x -> x[1], us), map(x -> x[2], us), ...
"""
macro unzip(args)
    args.head != :(=) && error("Expression needs to be of form `xs, ys, ... = us`")
    lhs, rhs = args.args
    items = isa(lhs, Symbol) ? [lhs] : lhs.args
    rhs_items = [:(map(x -> x[$i], $rhs)) for i in 1:length(items)]
    rhs_expand = Expr(:tuple, rhs_items...)
    esc(Expr(:(=), lhs, rhs_expand))
end
export @unzip

"""Contains the information defining a lenslet array."""
struct LensletSystem{T<:Real}
    eyebox_rectangle::Rectangle{T}
    subdivisions_of_eyebox::Tuple{Int64,Int64}
    lenses::Vector{ParaxialLens{T}}
    lattice_coordinates::Vector{Tuple{Int64,Int64}}
    lenslet_colors::Vector{String}
    projected_eyeboxes::Vector{SMatrix{3,4,T}}
    displayplanes::Vector{Plane{T,3}}
    lenslet_eye_boxes::Vector{SMatrix{3, 4, T, 12}}
    lenslet_eyebox_numbers::Vector{Int64}
    lenslet_eyebox_centers::Vector{SVector{3,T}}
    system_properties::Dict
    subdivided_eyebox_polys::Vector{SMatrix{3,4,T}}
end


# """ Computes the location of the optical center of the lens that will project the centroid of the display to the centroid of the eyebox. Normally the display centroid will be aligned with the geometric centroid of the lens, rather than the optical center of the lens."""
# function compute_optical_center(eyeboxcentroid,displaycenter,lensplane)
#     v = displaycenter - eyeboxcentroid
#     r = Ray(eyeboxcentroid,v)
#     intsct = point(closestintersection(surfaceintersection(lensplane,r),false)) #this seems complicated when you don't need CSG
#     @assert intsct !== nothing
#     return intsct
# end

""" Computes the location of the optical center of the lens that will project the centroid of the display to the centroid of the eyebox. Normally the display centroid will be aligned with the geometric centroid of the lens, rather than the optical center of the lens."""
function compute_optical_center(eyeboxcentroid,display_center,lens)
    lens_geometric_center = centroid(lens) #geometric and optical center coincide at this point
    v = eyeboxcentroid - lens_geometric_center
    r = Ray(display_center,v)
    #optical center may not lie inside the shape of the lens. surfaceintersection(lens,r) will return nothing in this case which will cause the setup_system code to crash. Instead intersect with the plane of the shape of the lens, which has infinite extent.
    surf = OpticSim.plane(lens)
    intsct = point(closestintersection(surfaceintersection(surf,r),false)) #this seems complicated when you don't need CSG
    @assert intsct !== nothing
    return intsct
end
export compute_optical_center

function localframe(a::ParaxialLens)
    if OpticSim.shape(a) isa ConvexPolygon
        return OpticSim.localframe(OpticSim.shape(a))
    else
        return nothing
    end
end

"""Paraxial lenses are created with a local transform of identity. This is confusing because the lens has a local transform and the shape in the lens has a local transform."""
function replace_optical_center(eyeboxcentroid,displaycenter,lens)
    world_optic_center = compute_optical_center(eyeboxcentroid,displaycenter,lens)
    world_to_local_lens = inv(localframe(lens))
    local_optic_center = world_to_local_lens * world_optic_center
    @assert isapprox(0.0, local_optic_center[3],atol = 1e-12) "z should be zero in the local frame but instead it is $(local_optic_center[3])"
    local_optic_center = SVector{2}((local_optic_center)[1:2]) #this is defined in the 2D lens plane so z = 0
    # need the untransformed convexpoly vertices
    return ParaxialLensConvexPoly(focallength(lens),shape(lens),local_optic_center)
end



"""returns display plane represented in world coordinates, and the center point of the display"""
function display_plane(lens) 
    center_point = centroid(lens) + -OpticSim.normal(lens)* OpticSim.focallength(lens)
    pln = Plane(-OpticSim.normal(lens), center_point, vishalfsizeu = .5, vishalfsizev = .5)
    return pln,center_point
end
export display_plane


""" eyebox rectangle represented as a 3x4 SMatrix with x,y,z coordinates in rows 1,2,3 respectively. By default the eyebox z is 0.0"""
function eyeboxpolygon(xsize,ysize)::SMatrix{3,4}
    xext,yext = (xsize,ysize) ./ 2.0 #divide by 2 to get min and max extensions
    return SMatrix{3,4}([xext;yext;0.0 ;; xext;-yext;0.0 ;; -xext;-yext;0.0 ;; -xext;yext;0.0])
end

subpoints(pts::AbstractVector,subdivisions) = reshape(collect(range(extrema(pts)...,subdivisions)),1,subdivisions)
export subpoints

subdivide(poly::AbstractMatrix,xsubdivisions,ysubdivisions) = subdivide(SMatrix{3,4}(poly),xsubdivisions,ysubdivisions)

"""poly is a two dimensional rectangle represented as a 3x4 SMatrix with x values in row 1 and y values in row 2. Returns a vector of length xsubdivisions*ysubdivisions containing all the subdivided rectangles."""
function subdivide(poly::T, xsubdivisions, ysubdivisions)::Vector{T} where{T<:SMatrix{3,4}}
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

function eyebox_number(tilecoords::NTuple{2,Int64},cluster::R,num_eyeboxes::Integer)::Int64 where {R<:AbstractLatticeCluster}
    _, tile_index = cluster_coordinates_from_tile_coordinates(cluster,tilecoords)
    eyeboxnum = mod(tile_index-1,num_eyeboxes) + 1
    return eyeboxnum
end

function eyebox_assignment(tilecoords::NTuple{2,Int64},cluster::R,eyeboxes) where{R<:AbstractLatticeCluster}
    #compute cluster
    eyeboxnum = eyebox_number(tilecoords,cluster,length(eyeboxes))
    return eyeboxes[eyeboxnum] #use linear index for Matrix.
end

"""given the eye box polygon and how it is to be subdivided computes the subdivided eyeboxes and centroids"""
function compute_lenslet_eyebox_data(eyeboxtransform,eyeboxpoly::SMatrix{3,4,T,12},subdivisions::Tuple{Int64,Int64})::Vector{SMatrix{3,4,Float64}} where{T} #can't figure out how to get the system to accept a type for an SMatrix of a Unitful quantity which is what eyeboxpoly is. Want to extract the number type from the Quantity type but this is not happening.
    subdivided_eyeboxpolys = subdivide(eyeboxpoly,subdivisions...)
    strippedpolys =  map(x-> ustrip.(mm,x),  subdivided_eyeboxpolys)
    subdivided_eyeboxpolys =[eyeboxtransform * x for x in strippedpolys] #these have units of mm which don't interact well with Transform.
    return subdivided_eyeboxpolys
end

function test_compute_lenslet_eyebox_data(eyeboxtransform,eyeboxpoly,subdivisions)
    subpolys = compute_lenslet_eyebox_data(eyeboxtransform,eyeboxpoly,subdivisions)
end

function project_eyebox_to_display_plane(eyeboxpoly::AbstractMatrix{T},lens,displayplane)::SMatrix where{T<:Real}
    rowdim,coldim = size(eyeboxpoly)
    rays = [Ray(SVector{rowdim}(point),opticalcenter(lens)-SVector{rowdim}(point)) for point in eachcol(eyeboxpoly)]

    points = collect([point(closestintersection(surfaceintersection(displayplane,ray),false)) for ray in rays])

    SMatrix{rowdim,coldim}(reinterpret(Float64,points)...)
end

"""System parameters for a typical HMD."""
function systemparameters() 
    return (eye_box = (10.0mm,9.0mm),fov = (90.0°,60.0°),eye_relief = 20.0mm, pupil_diameter = 4.0mm, display_sphere_radius = 40.0mm,min_fnumber = 1.6, pixel_pitch = .7μm) #these numbers give a reasonable system, albeit with 2000 lenslets. However displays are quite small (184x252)μm
end

"""System parameters which will yield fewer lenslets which will draw much faster"""
function smallsystemparameters()
    return (eye_box = (10.0mm,9.0mm),fov = (20.0°,20.0°),eye_relief = 20.0mm, pupil_diameter = 2mm, display_sphere_radius = 40.0mm,min_fnumber = 1.6, pixel_pitch = .7μm) #these numbers give a reasonable system, albeit with 2000 lenslets. However displays are quite small (184x252)μm
end
export smallsystemparameters


"""Coordinate frames for the eye/display system. The origin of this frame is at the geometric center of the eyeball. Positive Z axis is assumed to be the forward direction, the default direction of the eye's optical axis when looking directly ahead."""
function setup_coordinate_frames()
    eyeballframe = Transform()
    corneavertex = OpticSim.Data.cornea_to_eyecenter()
    eyeboxtransform = eyeballframe*OpticSim.translation(0.0,0.0,ustrip(mm,corneavertex))  #unfortunately can't use unitful values in transforms because the rotation and translation components would have different types which is not allowed in a Matrix.

    return (eyeball_frame = eyeballframe,eye_box_frame = eyeboxtransform)
end

"""repeat vec1 enough times to fill a vector of length(vec2)."""
function extend(vec1,vec2)
    numrepeats = Int64(ceil(length(vec1)/length(vec2)))
    remaining = mod(length(vec1),length(vec2))
    subdivs= similar(vec2,numrepeats*length(vec2)+remaining)
    empty!(subdivs)

    for i in 1:numrepeats-1   
        append!(subdivs,vec2)
    end
    append!(subdivs,vec2[1:remaining])

    return subdivs
end

"""Create lenslet system that will cover the eyebox and fov. If the numbers in the fov tuple do not have units the system assumes they represent an angle in radians. Is is better to use unitful angles though, either rad or °:
```
setup_system(eye_box,(1rad,1.2rad),...)
setup_system(eye_box,(1°,1.2°),...)
```
If you do this
```
setup_system(eye_box,(1,1.2),...)
```
the system will assume the angle is in radians."""
function setup_system(eye_box,fov,eye_relief,pupil_diameter,display_sphere_radius,min_fnumber,pixel_pitch)
    #All coordinates are ultimately transformed into the eyeball_frame coordinate systems
    (eyeball_frame,eye_box_frame) = setup_coordinate_frames()
        
    @info "Input parameters\n"

    parameters = (eye_box = eye_box,fov = fov,eye_relief = eye_relief,pupil_diameter = pupil_diameter,display_sphere_radius = display_sphere_radius,min_fnumber = min_fnumber,pixel_pitch = pixel_pitch)
    
    for key in keys(parameters) @info "$key = $(parameters[key])" end
    println("\n\n")

    eyeboxz = (eye_box_frame*SVector(0.0,0.0,0.0))[3]
    eyebox_plane = Plane([0.0,0.0,1.0],[0.0,0.0,eyeboxz])
    #get system properties
    props = system_properties(eye_relief,eye_box,fov,pupil_diameter,.2,11,pixelpitch = pixel_pitch, minfnumber = min_fnumber)

    subdivisions = props[:subdivisions] #tuple representing how the eyebox can be subdivided given the cluster used for the lenslets
    
    clusterdata = props[:cluster_data] 
    cluster = clusterdata[:cluster] #cluster that is repeated across the display to ensure continuous coverage of the eyebox and fov.
    focallength = ustrip(mm,props[:focal_length]) #strip units off because these don't work well with Transform
    
    #compute lenslets based on system properties. lattice_coordinates are the (i,j) integer lattice coordinates of the hexagonal lattice making up the display. These coordinates are used to properly assign color and subdivided eyebox to the lenslets.
    lenses,lattice_coordinates = spherelenslets(eyebox_plane,eye_relief,focallength,[0.0,0.0,1.0],display_sphere_radius,fov[1],fov[2],elementbasis(cluster))

    @info "$(length(lenses)) lenses required for this design"
    temp = display_plane.(lenses)
    displayplanes = [x[1] for x in temp]
    display_plane_centers = [x[2] for x in temp]

    lensletcolors = pointcolor.(lattice_coordinates,Ref(cluster))

    #compute subdivided eyebox polygons and assign to appropriate lenslets
    eyeboxpoly::SMatrix{3,4}  =  mm * (eye_box_frame * eyeboxpolygon(ustrip.(mm,eye_box)...)) #four corners of the eyebox frame which is assumed centered around the positive Z axis. Transformed to the eyeballframe. Have to switch back and forth between Unitful and unitless quantities because Transform doesn't work with Unitful values.

    @info "Subdividing eyebox polygon into $subdivisions sub boxes"
    subdivided_eyeboxpolys::Vector{SMatrix{3,4}} = compute_lenslet_eyebox_data(eye_box_frame,eyeboxpoly,subdivisions)
    test_compute_lenslet_eyebox_data(eye_box_frame,eyeboxpoly,subdivisions)

    #assign each lenslet one of the subdivided eyebox polygons
    lenslet_eye_boxes = eyebox_assignment.(lattice_coordinates,Ref(cluster),Ref(subdivided_eyeboxpolys))

    #assign each lenslet the index into subdivided_eyeboxpolys that corresponds to the subdivided eyebox polygon the lenslet is supposed to cover.
    lenslet_eyebox_numbers = eyebox_number.(lattice_coordinates,Ref(cluster),Ref(length(subdivided_eyeboxpolys)))

    #compute display rectangles that will cover the assigned eyebox polygons
    lensleteyeboxcenters = [Statistics.mean(eachcol(x)) for x in lenslet_eye_boxes] 

    
   
    # testresults = test_replace_optical_center.(lensleteyeboxcenters,lenslet_geometric_centers,display_plane_centers,lenses)
    # repropoints = [round.(x,digits = 2) for x in testresults]

    # @info "reprojected points $(repropoints)"

    #project eyebox into lenslet display plane and compute bounding box. This is the size of the display for this lenslet
    subdivs = extend(lenses,subdivided_eyeboxpolys)
    
    projected_eyeboxes = project_eyebox_to_display_plane.(subdivs,lenses,displayplanes) #repeate subdivided_eyeboxpolys enough times to cover all lenses
    eyebox_rectangle = Rectangle(ustrip(mm,eye_box[1]/2),ustrip(mm,eye_box[2]/2),[0.0,0.0,1.0],[0.0,0.0,eyeboxz], interface = opaqueinterface())

    projected_eyebox_centers = [Statistics.mean(eachcol(x)) for x in projected_eyeboxes]
    lenses = replace_optical_center.(lensleteyeboxcenters,projected_eyebox_centers,lenses)  #make new lenses with optical centers that will cause the centroid of the eyebox assigned to the lens to project to the center of the display plane.

    return LensletSystem{Float64}(eyebox_rectangle,subdivisions,lenses,lattice_coordinates,lensletcolors,projected_eyeboxes,displayplanes,lenslet_eye_boxes,lenslet_eyebox_numbers,lensleteyeboxcenters,props,subdivided_eyeboxpolys)

end
export setup_system







