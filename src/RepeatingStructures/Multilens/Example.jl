""" Typical properties for near eye VR display """
nominal_system_inputs() = (eye_relief = 20mm, fov = (5°,5°),eyebox = (10mm,8mm),display_radius = 125.0mm, pupil_diameter = 3.5mm,pixel_pitch = .9μm, minfnumber = 2.0, mtf = .2, cycles_per_degree = 11, max_display_size = 250μm, )
export nominal_system_inputs

xcoords(a::SMatrix{3,4}) = a[1,:]
ycoords(a::SMatrix{3,4}) = a[2,:]
zcoords(a::SMatrix{3,4}) = a[3,:]

function normal(a::SMatrix{3,4}) 
    side1 = a[:,4] - a[:,1]
    side2 = a[:,2] - a[:,1]
    return normalize(cross(side1,side2))
end

centroid(a::SMatrix{3,4}) = Statistics.mean(eachcol(a))

function matrix2rectangle(a::SMatrix{3,4,T}) where{T}
    center = Statistics.mean(eachcol(a))
    xmin,xmax = extrema(xcoords(a))
    halfsizeu = (xmax-xmin)/2.0
    ymin,ymax = extrema(ycoords(a))
    halfsizev = (ymax-ymin)/2.0
    Rectangle(halfsizeu,halfsizev,normal(a),center, interface = opaqueinterface())
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
function setup_nominal_system()
    (;eye_relief,fov,eyebox,display_radius,pupil_diameter,minfnumber,pixel_pitch) = nominal_system_inputs()
    setup_system(eyebox,fov,eye_relief,pupil_diameter,display_radius,minfnumber,pixel_pitch)
end
export setup_nominal_system

function test_projected_eyeboxes(sys = setup_nominal_system())
    lenses = sys.lenses
    display_eyeboxes = sys.projected_eyeboxes
    centers = [x for x in centroid.(shape.(lenses))] #use the geometric center of the lenses not the optical center_point
     display_centers = [Statistics.mean(eachcol(x)) for x in display_eyeboxes]
    rays = Ray.(display_centers, centers.-display_centers)
    ray_traces = Vector{LensTrace{Float64,3}}(undef,0)

    tracked_rays = [trace(CSGOpticalSystem(LensAssembly(lens),Rectangle(10.0,9.0,[0.0,0.0,1.0],[0.0,0.0,13.5],interface = opaqueinterface())),OpticalRay(ray,1.0,.530),trackrays = ray_traces) for (lens,ray) in zip(lenses,rays)]
    # all_rays = vcat(tracked_rays,ray_traces)
    return [x for x in ray_traces]
end
export test_projected_eyeboxes

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
    normal = SVector(0.0, 0, 1)
    hex = HexBasis1()
    verts =vcat(tilevertices((0, 0), hex), [0.0 0 0 0 0 0])
    surf = Plane(normal, SVector(0.0, 0, 10))

    project(verts, normal, surf)
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

"""assigns each lenslet/display subsystem a rectangular sub part of the eyebox"""
function test_eyebox_assignment(system = setup_nominal_system())
    (;eyeboxrect,
    lenses,
    lenslet_colors,
    lattice_coordinates,
    subdivisions_of_eyebox,
    projected_eyeboxes,
    displayplanes,
    lenslet_eyebox_numbers) = system
    Vis.draw() #clear screen
    Vis.draw!(eyeboxrect)
    for (lens,lattice_coordinates,color) in zip(lenses,lattice_coordinates,lenslet_colors)
        Vis.draw!(lens,color = color)
    end

    num_distinct_eyeboxes = reduce(*,subdivisions_of_eyebox) 
    colors = distinguishable_colors(num_distinct_eyeboxes)

    for (eyebox,plane,eyeboxnum) in zip(projected_eyeboxes,displayplanes,lenslet_eyebox_numbers)
        localframe = Transform(pointonplane(plane),normal(plane))
        transformedpoints = (inv(localframe)*eyebox)[1:2,:]
        pts = [SVector{2}(x...) for x in eachcol(transformedpoints)]
        # Vis.draw!(plane,color = colors[eyeboxnum])
        Vis.draw!(ConvexPolygon(localframe,pts),color = colors[eyeboxnum])
    end
end
export test_eyebox_assignment

function draw_system(system = setup_nominal_system())
    test_eyebox_assignment(system)
    Vis.draw!(test_projected_eyeboxes(system))
end
export draw_system

"""draw subdivided eyeboxes to make sure they are in the right place"""
function test_eyebox_subdivision(clear_screen = true)
    (;subdivided_eyebox_polys) = setup_nominal_system()
    @info "eyebox subdivisions $subdivided_eyebox_polys"
    for subdiv in subdivided_eyebox_polys
        @info "subdivided rectangle $(matrix2rectangle(subdiv))"
    end
    colors = distinguishable_colors(length(subdivided_eyebox_polys))

    if clear_screen
        Vis.draw()
    end
    
    for (eyebox,color) in zip(subdivided_eyebox_polys,colors)
        center = centroid(eyebox)
        Vis.draw!(matrix2rectangle(eyebox),color = color)
        Vis.draw!(leaf(Sphere(.1),Geometry.translation(center)),color = "white")
        r = Ray(center,normal(eyebox))
        Vis.draw!(r)
    end
    Vis.draw!(Sphere(.001)) #need to draw another 3d object not on the plane or numbers on axes are impossible to read
end
export test_eyebox_subdivision

"""For each display computes a ray from the display center to the geometric center of the lenslet. Then it computes the refracted ray and intersects it with the subdivided eyebox associated with the display. Draws all these rays."""
function test_ray_cast_from_display_center(clear_display = true)
    system = setup_nominal_system()
    (;eyeboxrect,lenses,projected_eyeboxes) = system

    if clear_display
        Vis.draw()
    end

    # draw_system(system)

    for (lens,display_eyebox) in zip(lenses,projected_eyeboxes)
        lens_geometric_center = centroid(shape(lens))
        display_center = centroid(display_eyebox)
        diff = lens_geometric_center-display_center

        r = OpticalRay(Ray(display_center,diff),1.0,.53)
       
        Vis.draw!(eyeboxrect)
        system = CSGOpticalSystem(LensAssembly(lens),eyeboxrect)
        Vis.draw!(r,color = "black")
        tr = trace(system,r)
        @info "trace of system = $tr"
        if tr !== nothing
            Vis.draw!(tr,color = "black")
        end
    end
end
export test_ray_cast_from_display_center


