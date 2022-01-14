""" Typical properties for near eye VR display """
nominal_system_inputs() = (eye_relief = 20mm, fov = (90°,60°),eyebox = (10mm,8mm),display_radius = 125.0mm, pupil_diameter = 3.5mm,pixel_pitch = .9μm, minfnumber = 2.0, mtf = .2, cycles_per_degree = 11, max_display_size = 250μm, )
export nominal_system_inputs


"""Example that shows how to call setup_system with typical values"""
function setupsystem_call()
    (;eye_relief,eyebox,display_radius,pupil_diameter,minfnumber,pixel_pitch) = nominal_system_inputs()
    fov = (10°,10°)
    setup_system(eyebox,fov,eye_relief,pupil_diameter,display_radius,minfnumber,pixel_pitch)
end
export setupsystem_call

function test_projected_eyeboxes()
    sys = setupsystem_call()
    lenses = sys.lenses
    display_eyeboxes = sys.projected_eyeboxes
    centers = centroid.(lenses)
    display_centers = [Statistics.mean(eachcol(x)) for x in display_eyeboxes]
    rays = Ray.(display_centers, centers.-display_centers)
    refracted_rays = [trace(LensAssembly(lens),OpticalRay(ray,1.0,.530)) for (lens,ray) in zip(lenses,rays)]
    ray_traces = Vector{LensTrace{Float64,3}}(undef,0)

    tracked_rays = [trace(CSGOpticalSystem(LensAssembly(lens),Rectangle(10.0,9.0,[0.0,0.0,1.0],[0.0,0.0,13.0],interface = opaqueinterface())),OpticalRay(ray,1.0,.530),trackrays = ray_traces) for (lens,ray) in zip(lenses,rays)]
    # all_rays = vcat(tracked_rays,ray_traces)
    return [x for x in ray_traces]
end
export test_projected_eyeboxes

"""Example call of system_properties function"""
function systemproperties_call()
    (;eye_relief,eyebox,display_radius,pupil_diameter,minfnumber,pixel_pitch) = nominal_system_inputs()
    fov = (10°,10°)
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
function test_eyebox_assignment()
    (;fov, eye_relief,eyebox,display_radius,pupil_diameter,minfnumber,pixel_pitch) = nominal_system_inputs()

    fov = (15°,15°)
    
    (;eyeboxrect,
    lenses,
    lenslet_colors,
    lattice_coordinates,
    subdivisions_of_eyebox,
    projected_eyeboxes,
    displayplanes,
    lenslet_eyebox_numbers) = setup_system(eyebox,fov,eye_relief,pupil_diameter,display_radius,minfnumber,pixel_pitch)
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

