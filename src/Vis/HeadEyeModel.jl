



function OpticSim.Vis.draw!(scene::Makie.LScene, h::ParaxialAnalysis.HeadEye.Head; draw_head=true, kwargs...) where {T<:Real}
    if (draw_head)
        scale_factor = 60.0;
        face = ParaxialAnalysis.HeadEye.transform_mesh(ParaxialAnalysis.HeadEye.get_face_model(), ParaxialAnalysis.HeadEye.tr(h)*Geometry.scale(scale_factor, scale_factor, scale_factor))
        Makie.mesh!(scene, face, color=RGBA(1.0,0.2,0.2, 0.1), transparency=true)
    end

    ParaxialAnalysis.HeadEye.draw_local_frame(scene, ParaxialAnalysis.HeadEye.tr(h), "Head Transform")

    for e in ParaxialAnalysis.HeadEye.eyes(h)
        Vis.draw!(scene, e, parent_transform = ParaxialAnalysis.HeadEye.tr(h))
    end
end

function OpticSim.Vis.draw!(scene::Makie.LScene, e::ParaxialAnalysis.HeadEye.Eye; parent_transform::Geometry.Transform = Transform(), kwargs...) where {T<:Real}

    ParaxialAnalysis.HeadEye.draw_local_frame(scene, parent_transform*ParaxialAnalysis.HeadEye.tr(e), ParaxialAnalysis.HeadEye.text(e), scale=2)
    Vis.draw!(scene, ParaxialAnalysis.HeadEye.pupil(e), parent_transform = parent_transform*ParaxialAnalysis.HeadEye.tr(e))
end

function OpticSim.Vis.draw!(scene::Makie.LScene, p::ParaxialAnalysis.HeadEye.Pupil; parent_transform::Geometry.Transform = Transform(), kwargs...) where {T<:Real}

    ParaxialAnalysis.HeadEye.draw_local_frame(scene, parent_transform*ParaxialAnalysis.HeadEye.tr(p), "pupil", scale=1.0)
    
    # half_size = size(p) / 2.0
    # eli = Ellipse(half_size, half_size, rotation(parent_transform)*unitZ3(), parent_transform*origin(tr(p)))
    # Vis.draw!(scene, eli, transparency=true, color=RGBA(0.2, 0.2, 0.2, 0.4))
end
