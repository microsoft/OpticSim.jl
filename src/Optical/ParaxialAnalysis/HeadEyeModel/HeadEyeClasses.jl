# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

# --------------------------------------------------------
# Pupil
# --------------------------------------------------------
mutable struct Pupil
    _tr::Transform
    _size::Float64

    function Pupil(tr::Transform, size::Float64)
        return new(tr, size)
    end
end

tr(p::Pupil) = p._tr
text(p::Pupil) = p._size
size(p::Pupil) = p._size
shape(p::Pupil) = Circle(size(p) / 2.0, forward(tr(p)), origin(tr(p)), interface = opaqueinterface())


# --------------------------------------------------------
# Eye
# --------------------------------------------------------
mutable struct Eye
    _tr::Transform
    _text::String
    _pupil::Pupil

    function Eye(tr::Transform, text::String)
        pupil = Pupil(Transform(SVector(0.0, 0.0, 8)), 4.0)

        return new(tr, text, pupil)
    end
end

tr(e::Eye) =  e._tr
text(e::Eye) =  e._text
pupil(e::Eye) =  e._pupil
function tr!(e::Eye, t::Transform) 
    e._tr = t
end
function lookAt!(e::Eye, dir::Vec3) 
    tr!(e, Transform(origin(tr(e)), dir))
end


# --------------------------------------------------------
# Head
# --------------------------------------------------------
mutable struct Head 
    _tr::Transform
    _eyes::Vector{Eye}

    function Head(tr::Transform)
        x_dist = 30.0
        y_dist = 3.0
        z_dist = 4.0
        left_eye = Eye(Transform(Vec3(-x_dist, y_dist, z_dist)), "Left Eye")
        right_eye = Eye(Transform(Vec3(+x_dist, y_dist, z_dist)), "Right Eye")
        return new(tr, [left_eye, right_eye])
    end
end

tr(h::Head) =  h._tr
eyes(h::Head) =  h._eyes
eye(h::Head, lr::Symbol) = h._eyes[eye_index(lr)]
eye_index(lr::Symbol) = lr == :left ? 1 : 2


function tr(h::Head, part::Symbol) 
    if (part == :head)
        return tr(h)
    elseif (part == :left_eye)
        return tr(h, :head) * tr(eyes(h)[1])
    elseif (part == :right_eye)
        return tr(h, :head) * tr(eyes(h)[2])
    elseif (part == :left_pupil)
        return return tr(h, :left_eye) * tr(pupil(eyes(h)[1]))
    elseif (part == :right_pupil)
        return tr(h, :right_eye) * tr(pupil(eyes(h)[2]))
    else
        @error "unrecognized part in Head::tr [$part]"
    end
end

function shape(h::Head, part::Symbol)
end

# -------------------------------------------------------


__face_model = nothing
function get_face_model()
    global __face_model
    if (__face_model === nothing)
        fn = joinpath(dirname(@__FILE__), "face_model.obj")
        __face_model = FileIO.load(fn)
    end
    return __face_model
end


function draw_arrow(scene::Makie.LScene, from::Vec3, to::Vec3; scale=1.0, color=Makie.MakieCore.Automatic())
    scale = scale * 0.5
    to2 = to - from
    Makie.arrows!(scene, 
        [from[1]], [from[2]], [from[3]], 
        [to2[1]], [to2[2]], [to2[3]], 
        linewidth=scale*0.25, 
        arrowsize=scale, 
        arrowcolor=color,
        linecolor=color, 
        transparency=true,
    )
end

function quat(tr::Transform)
    rot, _, _  =  decomposeRTS(tr)
	w = sqrt(1.0 + tr[1,1] + tr[2,2] + tr[3,3]) / 2.0
	w4 = (4.0 * w)
	x = (tr[3,2] - tr[2,3]) / w4
	y = (tr[1,3] - tr[3,1]) / w4
	z = (tr[2,1] - tr[1,2]) / w4
    return Makie.Quaternionf0(x, y, z, w)
end


function draw_local_frame(scene::Makie.LScene, tr::Transform, caption::String = "???"; scale=5.0)
    o = zero(Vec3)
    len = 6
    ux = unitX3() * len
    uy = unitY3() * len
    uz = unitZ3() * len
    l2w = local2world(tr)

    wo = l2w * o
    wux = l2w * ux
    wuy = l2w * uy
    wuz = l2w * uz

    Makie.scatter!(scene, [wo[1]], [wo[2]], [wo[3]], color=RGBA(0.9, 0.9, 0.1, 0.2), markersize=20, transparency=true)
    draw_arrow(scene, wo, wux, scale=scale, color=RGBA(0.9, 0.1, 0.1, 0.2))
    draw_arrow(scene, wo, wuy, scale=scale, color=color=RGBA(0.1, 0.9, 0.1, 0.2))
    draw_arrow(scene, wo, wuz, scale=scale, color=color=RGBA(0.1, 0.1, 0.9, 0.2))

    
    q = quat(tr)
    text_pos = l2w * (o + Vec3(0.0, 0.01, 0.0))
    Makie.text!( scene, 
        caption,
        rotation = q,
        position = text_pos, 
        color = RGBA(0.1, 0.1, 0.1, 0.2),
        align = (:left, :baseline),
        textsize = scale *2.0/5.0,
        space = :data
    )

    
end


# transform a mesh coordinates
function transform_mesh(mesh::GeometryBasics.Mesh, tr::Transform)
    tr = convert(Transform{Float32}, tr)

    vertices = GeometryBasics.coordinates(mesh)
    faces = GeometryBasics.faces(mesh)
    vertices = [GeometryBasics.Point3f0(tr * Vec3(p...)) for p in vertices]
    return GeometryBasics.Mesh(vertices, faces)
end




