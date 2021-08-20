# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

using ..OpticSim.ParaxialAnalysis.HeadEye

function example1()

    t = identitytransform()
    t = Transform(OpticSim.Vec4(t[:,1]), OpticSim.Vec4(t[:,2]), OpticSim.Vec4(t[:,3])*-1.0, OpticSim.Vec4(t[:,4]))
    head = HeadEye.Head(t)

    # HeadEye.lookAt!(HeadEye.eye(h, :left), SVector(-0.5, 0.0, 1.0))

    # -----------------------------------------------------------
    # mla stuff
    # -----------------------------------------------------------


    # compute center of mla - in front of left eye - oposite the tirection of the head
    left_eye_tr = Transform(origin(HeadEye.tr(head, :left_eye)), forward(HeadEye.tr(head, :head)))
    pos = local2world(left_eye_tr) * SVector(0.0, 0.0, 30.0) 
    mla_rt = Transform(pos, forward(HeadEye.tr(head, :head)) * -1.0)
    # eli = Ellipse(7.0, 7.0, unitZ3(), pos)
    # Vis.draw!(eli, transparency=true, color=RGBA(1.0, 0.2, 0.2, 0.4))

    # csg = HeadEye.csg_sphere(radius=20.0)
    # csg = HeadEye.csg_cylinder(radius = 20.0, added_rotation = rotationX(π/2.0))
     csg = HeadEye.csg_plane()

    # shapes_2d = HeadEye.get_shapes(HexBasis1(), resolution=(8,5), size=1.0)

     shapes_2d = HeadEye.get_shapes(RectangularBasis(), resolution=(15,15), size=1.0)

    shapes_3d = HeadEye.project(shapes_2d, csg)

    #create pupil lens

    # create the planar polygons and the lenses
    paraxial_lenses = Vector(undef, 0)
    emitters = Vector{Emitters.Sources.Source}(undef, 0)
    eyelens = HeadEye.leftpupil(head)

    for shape in shapes_3d

        # sanity test to change the lens optical center according to the shape's coordinates (relative to the center)
        # this one is for rectangles - need better algorithm for hexagons
        coord = HeadEye.coordinates(shape)

        # center_point = SVector(-3.0*coord[1], 0.0)
        # center_point = SVector(-3.0*coord[1], -3.0*coord[2])
        # center_point = SVector(0.0, -3.0*coord[2])
        center_point = SVector(0.0, 0.0)

        lens = HeadEye.build_paraxial_lens(
            shape, 
            local_center_point = center_point,
            parent_transform = mla_rt, 
            scale=0.99, 
            focal_length=1.0
        )
        push!(paraxial_lenses, OpticSim.LensAssembly(lens,eyelens))

        
        emitter = HeadEye.build_emitter(lens; distance=1.0, size=(.5, .5))
        push!(emitters, emitter)
    end

    detector = HeadEye.shape(HeadEye.retina(HeadEye.lefteye(head)))
    sys = OpticSim.CSGOpticalSystem.(paraxial_lenses, Ref(detector),500,500) #Julia idiomatic way of preventing broadcasting on an argument is to use Ref(arg)

    return head,Repeat.ArrayStructures.LensArray(sys,emitters)
end
export example1


function drawheadsystem()
   head, lensarray = example1()
    # -----------------------------------------------------------
    # rendering
    # -----------------------------------------------------------
    Vis.set_current_mode(:normal)


    Vis.draw(head, draw_head=true)
    emitters,lenses = (emitters(lensarray),lenses(lensarray))

    emitter = emitters[60]
    for onesys in sys
    # for (onesys,emitter) in zip(sys,emitters)
        Vis.drawtracerays!(onesys; raygenerator=emitter, trackallrays = true, colorbynhits = true, test = true, numdivisions = 100, drawsys = true, rayfilter = nothing)
        Vis.draw!(emitter, debug=false)
    end
end
export drawheadsystem
