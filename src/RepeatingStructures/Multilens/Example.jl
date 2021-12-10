""" Typical properties for near eye VR display """
nominal_system_properties() = (eye_relief = 20mm, fov = (90°,60°),eyebox = (10mm,8mm),display_radius = 125.0mm, pupil_diameter = 3.5mm,pixel_pitch = .9μm, minfnumber = 2.0, mtf = .2, cycles_per_degree = 11, max_display_size = 250μm, )
export nominal_system_properties


function drawspherelenslets()
    (;fov, eye_relief,eyebox,display_radius,pupil_diameter,pixel_pitch,minfnumber,mtf,cycles_per_degree) = nominal_system_properties()
    
    computedprops = systemproperties(eye_relief,eyebox,fov,pupil_diameter,mtf,cycles_per_degree,minfnumber = minfnumber,maxdisplaysize = max_display_size,pixelpitch = pixel_pitch)
    focallength = ustrip(mm,computedprops[:focal_length])
    lenses,_ = spherelenslets(Plane(0.0,0.0,1.0,0.0,0.0,18.0),focallength,[0.0,0.0,-1.0],ustrip(mm,display_radius),fov[1],fov[2],HexBasis1())
    Vis.draw!.(lenses)
    return nothing
end
export drawspherelenslets

function testspherepolygons()
    eyebox = Plane(0.0,0.0,1.0,0.0,0.0,12.0)
    dir = [0.0,0.0,-1.0]
    fovθ = 10.0°
    fovϕ = 5.0°
    lattice = HexBasis1()
    displayradius = 125.0mm
    eyerelief = 20mm

    return spherepolygons(eyebox,eyerelief,displayradius,dir,fovθ,fovϕ,lattice)
end
export testspherepolygons

