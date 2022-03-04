# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.


""" Typical properties for near eye VR display """
nominal_system_inputs() = (eye_relief = 20mm, fov = (5°,5°),eyebox = (10mm,8mm),display_radius = 125.0mm, pupil_diameter = 3.5mm,pixel_pitch = .9μm, minfnumber = 2.0, mtf = .2, cycles_per_degree = 11, max_display_size = 250μm, )
export nominal_system_inputs

xcoords(a::SMatrix{3,4}) = a[1,:]
ycoords(a::SMatrix{3,4}) = a[2,:]
zcoords(a::SMatrix{3,4}) = a[3,:]

function matnormal(a::SMatrix{3,4}) 
    side1 = a[:,4] - a[:,1]
    side2 = a[:,2] - a[:,1]
    return normalize(cross(side1,side2))
end

matcentroid(a::SMatrix{3,4}) = Statistics.mean(eachcol(a))

function matrix2rectangle(a::SMatrix{3,4,T}) where{T}
    center = Statistics.mean(eachcol(a))
    xmin,xmax = extrema(xcoords(a))
    halfsizeu = (xmax-xmin)/2.0
    ymin,ymax = extrema(ycoords(a))
    halfsizev = (ymax-ymin)/2.0
    Rectangle(halfsizeu,halfsizev,matnormal(a),center, interface = opaqueinterface())
end

function test_paraxial_lens()
    lens = ParaxialLensRect(10.0,1.0,1.0,[0.0,0.0,1.0],[0.0,0.0,0.0])
    displaypoint = [0.0,0.0,-2.5]
    direction = [0.0,0.0,1.0]
    r1 = Ray(displaypoint,direction)
    intsct1 = OpticSim.surfaceintersection(lens,r1)

    # compute the refracted ray and project this backwards to its intersection with the optical axis. This should be the same position as returned by virtualpoint()
    refrac1,_,_ = OpticSim.processintersection(OpticSim.interface(lens),OpticSim.point(intsct1),OpticSim.normal(lens),OpticalRay(r1,1.0,.55),OpticSim.TEMP_REF,OpticSim.PRESSURE_REF,false)
    r2 = Ray(-displaypoint,-direction)

    intsct2= OpticSim.surfaceintersection(lens,r2)

    # compute the refracted ray and project this backwards to its intersection with the optical axis. This should be the same position as returned by virtualpoint()
    refrac2,_,_ = OpticSim.processintersection(OpticSim.interface(lens),OpticSim.point(intsct2),OpticSim.normal(lens),OpticalRay(r2,1.0,.55),OpticSim.TEMP_REF,OpticSim.PRESSURE_REF,false)

    # isapprox(-refrac2[1],refrac1[1]), refrac2[1], refrac1[1]
    Vis.draw(lens)
    Vis.draw!([Ray(OpticSim.point(intsct1)+ [.5,.5,0.0],refrac1),r1]) #offset first ray so it can be seen. The other ray will be written over it and obscure it.
    Vis.draw!([Ray(OpticSim.point(intsct2),refrac2),r2], color = "green") #visualization should have two rays in opposite directions coming out of lens plane but only see one. Error in ParaxialLens.
end
export test_paraxial_lens

"""Example that shows how to call setup_system with typical values"""
function setup_nominal_system(;eyebox_subdivisions = nothing)::LensletSystem
    (;eye_relief,fov,eyebox,display_radius,pupil_diameter,minfnumber,pixel_pitch) = nominal_system_inputs()
    setup_system(eyebox,fov,eye_relief,pupil_diameter,display_radius,minfnumber,pixel_pitch,eyebox_subdivisions = eyebox_subdivisions)
end
export setup_nominal_system

function compute_eyebox_rays(sys = setup_nominal_system())
    (;lenses,projected_eyeboxes,eyebox_rectangle) = sys
    centers = [x for x in centroid.(shape.(lenses))] #use the geometric center of the lenses not the optical center_point
    display_centers = matcentroid.(projected_eyeboxes)
    rays = Ray.(display_centers, centers.-display_centers)
    ray_traces = Vector{LensTrace{Float64,3}}(undef,0)

    tracked_rays = [trace(CSGOpticalSystem(LensAssembly(lens),eyebox_rectangle),OpticalRay(ray,1.0,.530),trackrays = ray_traces) for (lens,ray) in zip(lenses,rays)]

    return [x for x in ray_traces]
end
export compute_eyebox_rays

"""Example call of system_properties function"""
function systemproperties_call()
    (;eye_relief,fov,eyebox,pupil_diameter,minfnumber,pixel_pitch) = nominal_system_inputs()
    system_properties(eye_relief,eyebox,fov,pupil_diameter,.2,20,minfnumber = minfnumber,pixelpitch = pixel_pitch)
end
export systemproperties_call

"""Computes the lenslets on the display sphere and draws them."""
function drawspherelenslets()
    (;fov, eye_relief,eyebox,display_radius,pupil_diameter,pixel_pitch,minfnumber,mtf,cycles_per_degree) = nominal_system_inputs()
    
    computedprops = system_properties(eye_relief,eyebox,fov,pupil_diameter,mtf,cycles_per_degree,minfnumber = minfnumber,maxdisplaysize = 250μm,pixelpitch = pixel_pitch)
    focallength = ustrip(mm,computedprops[:focal_length])
    lenses,_ = spherelenslets(Plane(0.0,0.0,1.0,0.0,0.0,18.0),eye_relief,focallength,[0.0,0.0,1.0],display_radius,fov[1],fov[2],HexBasis1())
    Vis.draw() #clear screen
    Vis.draw!.(lenses)
    return nothing
end
export drawspherelenslets

"""returns the hex polygons on the sphere that will eventually be turned into hexagonal lenslets"""
function testspherepolygons()
    eyebox = Plane(0.0,0.0,1.0,0.0,0.0,12.0)
    dir = [0.0,0.0,1.0]
    fovθ = 10.0°
    fovϕ = 5.0°
    lattice = HexBasis1()
    displayradius = 125.0mm
    eyerelief = 20mm

    return spherepolygons(eyebox,eyerelief,displayradius,dir,fovθ,fovϕ,lattice)
end
export testspherepolygons


function testprojectonplane()
    verts = tilevertices(HexBasis1())
    verts = vcat(verts,[0 for _ in 1:6]')
    projectonbestfitplane(SMatrix{size(verts)...}(verts),[0.0,0.0,1.0])
end
export testprojectonplane

function testproject()
    norml = SVector(0.0, 0, 1)
    hex = HexBasis1()
    verts =vcat(tilevertices((0, 0), hex), [0.0 0 0 0 0 0])
    surf = Plane(norml, SVector(0.0, 0, 10))

    project(verts, norml, surf)
end
export testproject

function testprojection()
    pts = spherepoints(1.0,-.2,-.2,1.0,1.1)
    surf = Plane(0.0,0.0,-1.0,0.0,0.0,0.0)
    dir = [0.0,0.0,-1.0]
    project(pts,dir,surf)
end
export testprojection

"""returns hexagonal lenslets on the display sphere"""
function testspherelenslets()
    (;fov, eye_relief,eyebox,display_radius,pupil_diameter,pixel_pitch,minfnumber,mtf,cycles_per_degree,max_display_size) = nominal_system_inputs()
    
    fov = (10°,10°)
    computedprops = system_properties(eye_relief,eyebox,fov,pupil_diameter,mtf,cycles_per_degree,minfnumber = minfnumber,maxdisplaysize = max_display_size,pixelpitch = pixel_pitch)
    focallength = ustrip(mm,computedprops[:focal_length])
    spherelenslets(Plane(0.0,0.0,1.0,0.0,0.0,18.0),eye_relief,focallength,[0.0,0.0,1.0],display_radius,fov[1],fov[2],HexBasis1())
end
export testspherelenslets

function draw_projected_corners(system = setup_nominal_system())
    (;lenslet_eye_boxes, subdivisions_of_eyebox,lenses,displayplanes,projected_eyeboxes) = system
    colors = distinguishable_colors(reduce(*,subdivisions_of_eyebox))

    for (lenslet_eyebox,lens,display_plane,projected_eyebox) in zip(lenslet_eye_boxes,lenses,displayplanes,projected_eyeboxes)
        center = opticalcenter(lens)
        for (eyeboxpt,projected_point) in zip(eachcol(lenslet_eyebox),eachcol(projected_eyebox))
            r = Ray(eyeboxpt,center-eyeboxpt)
            intsct = surfaceintersection(display_plane,r)
            closest = closestintersection(intsct)
            display_intersection = point(closest)
            @assert isapprox(display_intersection,projected_point) "error display_intersection $display_intersection projected_point $projected_point"
            lenstrace = LensTrace(OpticalRay(r,1.0,.5),closest)
            Vis.draw!(lenstrace)
        end
    end
end
export draw_projected_corners

"""assigns each lenslet/display subsystem a rectangular sub part of the eyebox"""
function draw_eyebox_assignment(system = setup_nominal_system(),clear_screen = true;draw_eyebox = true)
    (;eyebox_rectangle,
    lenses,
    lenslet_colors,
    lattice_coordinates,
    subdivisions_of_eyebox,
    projected_eyeboxes,
    displayplanes,
    lenslet_eyebox_numbers) = system
    if clear_screen
        Vis.draw() #clear screen
    end

    if draw_eyebox
        Vis.draw!(eyebox_rectangle)
    end

    for (lens,lattice_coordinates,color) in zip(lenses,lattice_coordinates,lenslet_colors)
        Vis.draw!(lens,color = color)
    end

    num_distinct_eyeboxes = reduce(*,subdivisions_of_eyebox) 
    colors = distinguishable_colors(num_distinct_eyeboxes)

    #THIS is almost certainly wrong
    # for (eyebox,plane,eyeboxnum,lens) in zip(projected_eyeboxes,displayplanes,lenslet_eyebox_numbers,lenses)
    #     localframe = OpticSim.translation(-OpticSim.normal(lens))*OpticSim.localframe(shape(lens)) #the local frame of the lens is, unfortunately, stored only in the ConvexPoly shape. No other lens type has a local frame which is terrible design. Should be fixed.
    #     transformedpoints = (inv(localframe)*eyebox)[1:2,:]
    #     pts = [SVector{2}(x...) for x in eachcol(transformedpoints)]
    #     # Vis.draw!(plane,color = colors[eyeboxnum])
    #     Vis.draw!(ConvexPolygon(localframe,pts),color = colors[eyeboxnum])
    # end
end
export draw_eyebox_assignment

#compute the ray from each lenslet's eyebox center to the geometric center of the lenslet. Trace this ray through the lens and draw the refracted ray. The refracted ray should line up exactly with -normal(lenslet).
function draw_eyebox_rays(system = setup_nominal_system())
    (;lenses, lenslet_eyebox_centers,) = system

    for (lens,eyebox_center) in zip(lenses,lenslet_eyebox_centers)
        r = Ray(centroid(lens),centroid(lens)-eyebox_center)
        refrac_ray = ray(trace(LensAssembly(lens),OpticalRay(r,1.0,.5)))
        Vis.draw!(refrac_ray,color = "yellow")
    end
end

function draw_eyebox_polygons(system = setup_nominal_system())
    (;projected_eyebox_polygons) = system
    for poly in projected_eyebox_polygons
        Vis.draw!(poly)
    end
end
export draw_eyebox_polygons

function draw_system(system = setup_nominal_system())
    draw_eyebox_assignment(system,draw_eyebox=false)
    draw_subdivided_eyeboxes(system,false)
    # Vis.draw!(compute_eyebox_rays(system))
    draw_eyebox_rays(system)
    draw_projected_corners(system)
    draw_eyebox_polygons(system)
end
export draw_system

"""draw subdivided eyeboxes to make sure they are in the right place"""
function draw_subdivided_eyeboxes(system = setup_nominal_system(), clear_screen = true)
    (;subdivided_eyebox_polys) = system

    colors = distinguishable_colors(length(subdivided_eyebox_polys))

    if clear_screen
        Vis.draw()
    end
    
    for (eyebox,color) in zip(subdivided_eyebox_polys,colors)
        center = matcentroid(eyebox)
        Vis.draw!(matrix2rectangle(eyebox),color = color)
        r = Ray(center,matnormal(eyebox))
        Vis.draw!(r)
        Vis.draw!(leaf(Sphere(.1),Geometry.translation(center)),color = "white")
    end
end
export draw_subdivided_eyeboxes


function test_project_eyebox_to_display_plane()
    boxpoly = SMatrix{3,4}(
        -1.0,1,0.0,
        1.0,1.0,0.0,
        1.0,-1.0,0.0,
        -1.0,-1.0,0.0
        )
    correct_answer = SMatrix{3,4}(
        1.3,1.0, 11.5, 
        1.0,1.0, 11.5, 
        1.0,1.3, 11.5, 
        1.3,1.3, 11.5)
    fl = 1.5
    lenscenter = [1.0,1.0,10.0]
    lens = ParaxialLensRect(fl,.5,.5,[0.0,0.0,1.0],lenscenter)
    displayplane = Plane([0.0,0.0,1.0],[0.0,0.0,lenscenter[3]+fl])

    answer,polygons = project_eyebox_to_display_plane(boxpoly,lens,displayplane)
    @assert isapprox(answer,correct_answer) "computed points $answer"
    @info "eyebox points $answer"
    @info "eyebox polygons $polygons"
end
export test_project_eyebox_to_display_plane

"""For each display computes a ray from the display center to the geometric center of the lenslet. Then it computes the refracted ray and intersects it with the subdivided eyebox associated with the display. Draws all these rays. The rays from the displays should hit the center points of the subdivided eyeboxes"""
function test_ray_cast_from_display_center(clear_display = true)
    system = setup_nominal_system()
    (;eyebox_rectangle,lenses,projected_eyeboxes) = system

    if clear_display
        Vis.draw()
    end

    # draw_system(system)
    draw_subdivided_eyeboxes(system,false)

    for (lens,display_eyebox) in zip(lenses,projected_eyeboxes)
        lens_geometric_center = centroid(lens)

        display_center = matcentroid(display_eyebox)
        diff = lens_geometric_center-display_center

        r = OpticalRay(Ray(display_center,diff),1.0,.53)
       
        # Vis.draw!(eyebox_rectangle)
        system = CSGOpticalSystem(LensAssembly(lens),eyebox_rectangle)
        Vis.draw!(r,color = "black")
        tr = trace(system,r)

        if tr !== nothing
            Vis.draw!(tr,color = "black")
        end
    end
end
export test_ray_cast_from_display_center


"""refracts the ray from the center of the display passing through the geometric center of the lens and then computes the closes point on this ray to the eyeboxcentroid. If the lens optical center has been calculated correctly the ray should pass through the eyeboxcentroid"""
function test_replace_optical_center(eyeboxcentroid,lensgeometriccenter,displaycenter,lens)
    rdir = lensgeometriccenter-displaycenter
    r = Ray(displaycenter,rdir)

    assy = LensAssembly(lens)
    refracted_trace = trace(assy,OpticalRay(r,1.0,.530))

    r2 = Ray(eyeboxcentroid,lensgeometriccenter-eyeboxcentroid)

    #this doesn't work because of a bug in processintersection for paraxial interfaces. Ray points in the wrong direction.
    reverse_trace = trace(assy,OpticalRay(r2,1.0,.530))
    reverse_point = closestpointonray(ray(reverse_trace),displaycenter)

    return closestpointonray(ray(refracted_trace),eyeboxcentroid)
end
export test_replace_optical_center

function center_data()
    a = HexBasis1()
    verts = 5*tilevertices(a)
    vecofmat = collect(reinterpret(reshape,SVector{2,eltype(verts)},verts))
     
    lens = ParaxialLensConvexPoly(
        1.0,
        Geometry.translation(1.0,0.0,1.0)*Geometry.rotation(0.0,Float64(π),0.0),
        vecofmat,
        SVector(0.0,0.0)
    )
    display_center = [1.0,0.0,2.0]
    eyebox_center = [0.0,0.0,0.0]
    return lens,eyebox_center,display_center
end
export center_data

"""verify that paraxial lens properly transorms rays into local coordinate frame to compute intersections with convex poly shape"""
function test_paraxial_convex_poly()
    lens,eyebox_center,display_center = center_data()

    @info "normal of lens $(OpticSim.normal(lens))"
    @info "center of lens $(centroid(lens)) center of lens geometry $(centroid(shape(lens)))"
 
    display_ray = OpticalRay(Ray([1.0,0.0,2.0],[0.0,0.0,-1.0]),1.0,.53)
    
    ltrace = trace(LensAssembly(lens),display_ray)
    @info "trace with optical center aligned with geometric center $ltrace"
    new_center = compute_optical_center(eyebox_center,display_center,lens)
    @info "offset center $new_center"
    displaced_lens = replace_optical_center(eyebox_center,display_center,lens)
    @info "optical center of displacedlens $(opticalcenter(displaced_lens)) center of geometry $(centroid(shape(displaced_lens)))"

    dtrace = trace(LensAssembly(displaced_lens),display_ray)
    @info "trace with optical center displaced from geometric center  $dtrace"
    refracted_ray = ray(dtrace)
    intsct = surfaceintersection(Plane([0.0,0.0,1.0],eyebox_center),refracted_ray)
    
    @info "intsct should be [0,0,0] $(point(OpticSim.lower(intsct)))" #ray was computed to refract from the geometric lens center to the point [0,0,0] on the eyebox plane
    return intsct

end
export test_paraxial_convex_poly


