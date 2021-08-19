# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

module Examples

using ..OpticSim
using ..OpticSim.Geometry
using ..OpticSim.Repeat
using ..OpticSim.Emitters
using ..OpticSim.ParaxialAnalysis
using ..OpticSim.ParaxialAnalysis.HeadEye

using StaticArrays
using Colors
import Makie


emitterarray(emitters, transform = identitytransform()) = Emitters.Sources.CompositeSource(transform, emitters)

function example1()

    t = identitytransform()
    t = Transform(OpticSim.Vec4(t[:,1]), OpticSim.Vec4(t[:,2]), OpticSim.Vec4(t[:,3])*-1.0, OpticSim.Vec4(t[:,4]))
    h = HeadEye.Head(t)

    # HeadEye.lookAt!(HeadEye.eye(h, :left), SVector(-0.5, 0.0, 1.0))

    # -----------------------------------------------------------
    # mla stuff
    # -----------------------------------------------------------


    # compute center of mla - in front of left eye - oposite the tirection of the head
    left_eye_tr = Transform(origin(HeadEye.tr(h, :left_eye)), forward(HeadEye.tr(h, :head)))
    pos = local2world(left_eye_tr) * SVector(0.0, 0.0, 30.0) 
    mla_rt = Transform(pos, forward(HeadEye.tr(h, :head)) * -1.0)
    # eli = Ellipse(7.0, 7.0, unitZ3(), pos)
    # Vis.draw!(eli, transparency=true, color=RGBA(1.0, 0.2, 0.2, 0.4))

    csg = HeadEye.csg_sphere(radius=20.0)
    # csg = HeadEye.csg_cylinder(radius = 20.0, added_rotation = rotationX(π/2.0))
    # csg = HeadEye.csg_plane()

    # shapes_2d = HeadEye.get_shapes(HexBasis1(), resolution=(8,5), size=1.0)

     shapes_2d = HeadEye.get_shapes(RectangularBasis(), resolution=(5,5), size=1.0)

    shapes_3d = HeadEye.project(shapes_2d, csg)

    # create the planar polygons and the lenses
    paraxial_lenses = Vector{ParaxialLens}(undef, 0)
    emitters = Vector{Emitters.Sources.Source}(undef, 0)
    for shape in shapes_3d

        # sanity test to change the lens optical center according to the shape's coordinates (relative to the center)
        # this one is for rectangles - need better algorithm for hexagons
        coord = HeadEye.coordinates(shape)
        center_point = SVector(-3.0*coord[1], 0.0)
        center_point = SVector(-3.0*coord[1], -3.0*coord[2])
        # center_point = SVector(0.0, -3.0*coord[2])
        # center_point = SVector(0.0, 0.0)

        lens = HeadEye.build_paraxial_lens(
            shape, 
            local_center_point = center_point,
            parent_transform = mla_rt, 
            scale=0.99, 
            focal_length=30.0
        )
        push!(paraxial_lenses, lens)

        emitter = HeadEye.build_emitter(lens; distance=1.0, size=(1.0, 1.0))
        push!(emitters, emitter)
    end

    mla = LensAssembly(paraxial_lenses...)
    display = emitterarray(emitters)

    pupil = HeadEye.pupil(HeadEye.eye(h, :left))
    pupil_tr = HeadEye.tr(h, :left_pupil)
    detector = Circle(HeadEye.pupilsize(pupil) / 2.0, forward(pupil_tr), origin(pupil_tr), interface = opaqueinterface())
    sys = CSGOpticalSystem(mla, detector)



    # -----------------------------------------------------------
    # rendering
    # -----------------------------------------------------------
    Vis.set_current_mode(:normal)

    if (true)
        Vis.drawtracerays(sys; raygenerator=display, trackallrays = true, colorbynhits = true, test = true, numdivisions = 100, drawgen = false)
        # Vis.draw!(sys)

        Vis.draw!(h, draw_head=false)

        Vis.draw!(display, debug=false)
    end
end
export example1




end # module Test