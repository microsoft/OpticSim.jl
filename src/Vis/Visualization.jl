# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

using ..OpticSim
using ..OpticSim: euclideancontrolpoints, evalcsg, vertex, makiemesh, detector, centroid, lower, upper, intervals, α
using ..OpticSim.Geometry

using Unitful
using ImageView
using Images
using ColorTypes
using ColorSchemes
using StaticArrays
using LinearAlgebra
import Makie
import GeometryBasics
import Plots
import Luxor
using FileIO

#############################################################################

function drawcurve!(canvassize::Int, curve::Spline{P,S,N,M}, numpoints::Int; linewidth = 0.5, curvecolor = RGB(1, 0, 1), controlpointcolor = RGB(0, 0, 0), controlpointsize = 5, controlpolygoncolor = RGB(0, 0, 0)) where {P,S,N,M}
    step = 1.0 / numpoints

    point1 = point(curve, 0.0) .* canvassize
    # println("point1 $point1")
    point1 = [point1[1], canvassize - point1[2]] # flip y because of the coordinate system SVG uses

    Luxor.setline(linewidth)
    Luxor.sethue(curvecolor)

    for u in step:step:1.0
        point2 = point(curve, u) .* canvassize
        point2 = [point2[1], canvassize - point2[2]]

        Luxor.line(Luxor.Point(point1...), Luxor.Point(point2...), :fillstroke)
        point1 = copy(point2)
    end

    Luxor.sethue(controlpointcolor)
    # draw control points
    controlpoints = euclideancontrolpoints(curve)
    for point in controlpoints
        pt = point .* canvassize
        pt = [pt[1], canvassize - pt[2]]
        Luxor.circle(Luxor.Point(pt...), controlpointsize, :fill)
    end

    Luxor.sethue(controlpolygoncolor)
    for i in 1:(size(controlpoints)[1] - 1)
        pt1 = controlpoints[i] .* canvassize
        pt1 = [pt1[1], canvassize - pt1[2]]
        pt2 = controlpoints[i + 1] .* canvassize
        pt2 = [pt2[1], canvassize - pt2[2]]
        Luxor.line(Luxor.Point(pt1...), Luxor.Point(pt2...), :fillstroke)
    end
end

function drawcurves(curves::Vararg{Spline{P,S,N,M}}; numpoints::Int = 200, canvassize::Int = 2000) where {P,S,N,M}
    canvas = Luxor.Drawing(canvassize, canvassize)
    Luxor.background("white")
    drawcurve!(canvassize, curves[1], numpoints, linewidth = 5, curvecolor = RGB(0, 1, 1))
    for curve in curves[2:end]
        drawcurve!(canvassize, curve, numpoints, linewidth = 1, controlpointcolor = RGB(1, 0, 0), controlpointsize = 3, controlpolygoncolor = RGB(0, 1, 0))
    end
    Luxor.finish()
    Luxor.preview()
end

#############################################################################

# all functions follow the pattern draw(obj) and draw!(scene, obj) where the first case draws the object in a blank scene and displays it, and
# the second case draws the object in an existing scene, draw!(obj) can also be used to draw the object in the current scene

global current_main_scene = nothing
global current_layout_scene = nothing
global current_3d_scene = nothing
global current_mode = nothing           # modes:    nothing, :default  -> Original Vis beheviour    
                                        #           :pluto             -> support pluto notebooks 
                                        #           :docs              -> support documenter figures 

# added the following 2 functions to allow us to hack the drawing mechanisim while in a pluto notebook
set_current_main_scene(scene) = (global current_main_scene = scene)
set_current_3d_scene(lscene) = (global current_3d_scene = lscene)

get_current_mode() = begin global current_mode; return current_mode end
set_current_mode(mode) = (global current_mode = mode)

show(image) = imshow(image)
display(scene = current_main_scene) = begin 
    global current_mode; 
    if (get_current_mode() == :pluto || get_current_mode() == :docs)
        return scene
    end
    Makie.display(scene)
end

"""
    scene(resolution = (1000, 1000))

Create a new Makie scene with the given resolution including control buttons.
"""
function scene(resolution = (1000, 1000))
    @assert resolution[1] > 0 && resolution[2] > 0

    scene, layout = Makie.layoutscene(resolution = resolution)
    global current_main_scene = scene
    global current_layout_scene = layout
    lscene = layout[1, 1] = Makie.LScene(scene, scenekw = (camera = Makie.cam3d_cad!, axis_type = Makie.axis3d!, raw = false))
    global current_3d_scene = lscene

    # in these modes we want to skip the creation of the utility buttons as these modes are not interactive
    if (get_current_mode() == :pluto || get_current_mode() == :docs)
        return scene, lscene
    end

    threedbutton = Makie.Button(scene, label = "3D", buttoncolor = RGB(0.8, 0.8, 0.8), height = 40, width = 80)
    twodxbutton = Makie.Button(scene, label = "2D-x", buttoncolor = RGB(0.8, 0.8, 0.8), height = 40, width = 80)
    twodybutton = Makie.Button(scene, label = "2D-y", buttoncolor = RGB(0.8, 0.8, 0.8), height = 40, width = 80)
    savebutton = Makie.Button(scene, label = "Screenshot", buttoncolor = RGB(0.8, 0.8, 0.8), height = 40, width = 160)

    Makie.on(threedbutton.clicks) do nclicks
        make3d(lscene)
        yield()
    end

    Makie.on(twodybutton.clicks) do nclicks
        make2dy(lscene)
        yield()
    end

    Makie.on(twodxbutton.clicks) do nclicks
        make2dx(lscene)
        yield()
    end

    Makie.on(savebutton.clicks) do nclicks
        Vis.save("screenshot.png")
        yield()
    end

    layout[2, 1] = Makie.grid!(hcat(threedbutton, twodxbutton, twodybutton, savebutton), tellwidth = false, tellheight = true)
    return scene, lscene
end

function make3d(scene::Makie.LScene = current_3d_scene)
    s = scene.scene
    # use 3d camera
    Makie.cam3d_cad!(s)
    # reset scene rotation
    s.transformation.rotation[] = Makie.Quaternion(0.0, 0.0, 0.0, 1.0)
    # show all the axis ticks
    s[Makie.OldAxis].attributes.showticks[] = (true, true, true)
    # reset tick and axis label rotation and position
    rot = (Makie.Quaternion(0.0, 0.0, -0.7071067811865476, -0.7071067811865475), Makie.Quaternion(0.0, 0.0, 1.0, 0.0), Makie.Quaternion(0.0, 0.7071067811865475, 0.7071067811865476, 0.0))
    s[Makie.OldAxis].attributes.ticks.rotation[] = rot
    s[Makie.OldAxis].attributes.names.rotation[] = rot
    s[Makie.OldAxis].attributes.ticks.align = ((:left, :center), (:right, :center), (:right, :center))
    s[Makie.OldAxis].attributes.names.align = ((:left, :center), (:right, :center), (:right, :center))
    # reset scene limits to automatic
    s.limits[] = Makie.Automatic()
    Makie.update!(s)
end

function make2dy(scene::Makie.LScene = current_3d_scene)
    s = scene.scene
    # use 2d camera
    Makie.cam2d!(s)

    scene_transform = Makie.qrotation(Makie.Vec3f0(0, 1, 0), 0.5pi)
    scene_transform_inv = Makie.qrotation(Makie.Vec3f0(0, 1, 0), -0.5pi)    # to use with the ticks and names

    # set rotation to look onto yz plane
    s.transformation.rotation[] = scene_transform
    # hide x ticks

    # there is a bug in Makie 0.14.2 which cause an exception setting the X showticks to false. 
    # we work wround it by making sure the labels we want to turn off are ortogonal to the view direction 
    # s[Makie.OldAxis].attributes.showticks[] = (false, true, true)
    s[Makie.OldAxis].attributes.showticks[] = (true, true, true)

    # set tick and axis label rotation and position
    s[Makie.OldAxis].attributes.ticks.rotation[] = (0.0, scene_transform_inv, scene_transform_inv)
    s[Makie.OldAxis].attributes.names.rotation[] = s[Makie.OldAxis].attributes.ticks.rotation[]
    s[Makie.OldAxis].attributes.ticks.align = ((:right, :center), (:right, :center), (:center, :right))
    s[Makie.OldAxis].attributes.names.align = ((:left, :center), (:left, :center), (:center, :left))
    # update the scene limits automatically to get true reference values
    s.limits[] = Makie.Automatic()
    Makie.update_limits!(s)
    # manually set the scene limits to draw the axes correctly
    o, w = Makie.origin(s.data_limits[]), Makie.widths(s.data_limits[])
    s.limits[] = Makie.FRect3D((1000.0f0, o[2], o[3]), (w[2], w[2], w[3]))
    # set the eye (i.e. light) position to behind the camera
    s.camera.eyeposition[] = (0, 0, -100)
    Makie.update!(s)
end

function make2dx(scene::Makie.LScene = current_3d_scene)
    s = scene.scene
    # use 2d camera
    Makie.cam2d!(s)

    scene_transform= Makie.qrotation(Makie.Vec3f0(0, 0, 1), 0.5pi) * Makie.qrotation(Makie.Vec3f0(1, 0, 0), 0.5pi)
    scene_transform_inv=Makie.qrotation(Makie.Vec3f0(1, 0, 0), -0.5pi) * Makie.qrotation(Makie.Vec3f0(0, 0, 1), -0.5pi) 

    # set rotation to look onto yz plane
    s.transformation.rotation[] = scene_transform
    # hide y ticks

    # there is a bug in Makie 0.14.2 which cause an exception setting the X showticks to false. 
    # we work wround it by making sure the labels we want to turn off are ortogonal to the view direction 
    s[Makie.OldAxis].attributes.showticks[] = (true, false, true)

    # set tick and axis label rotation and position
    s[Makie.OldAxis].attributes.ticks.rotation[] = (scene_transform_inv, 0.0, scene_transform_inv)
    s[Makie.OldAxis].attributes.names.rotation[] = s[Makie.OldAxis].attributes.ticks.rotation[]
    s[Makie.OldAxis].attributes.ticks.align = ((:right, :center), (:right, :center), (:center, :center))
    s[Makie.OldAxis].attributes.names.align = ((:left, :center), (:right, :center), (:center, :center))
    # update the scene limits automatically to get true reference values
    s.limits[] = Makie.Automatic()
    Makie.update_limits!(s)
    # manually set the scene limits to draw the axes correctly
    o, w = Makie.origin(s.data_limits[]), Makie.widths(s.data_limits[])
    s.limits[] = Makie.FRect3D((o[1], -1000.0f0, o[3]), (w[1], w[1], w[3]))
    # set the eye (i.e. light) position to behind the camera
    s.camera.eyeposition[] = (0, 0, -100)
    Makie.update!(s)
end

"""
    draw(ob; resolution = (1000, 1000), kwargs...)

Draw an object in a new scene.
`kwargs` depends on the object type.
"""
function draw(ob; resolution = (1000, 1000), kwargs...)
    scene, lscene = Vis.scene(resolution)
    draw!(lscene, ob; kwargs...)
    display(scene)

    if (get_current_mode() == :pluto || get_current_mode() == :docs)
        return scene
    end
end

"""
    draw!([scene = currentscene], ob; kwargs...)

Draw an object in an existing scene.
`kwargs` depends on the object type.
"""
function draw!(ob; kwargs...)
    if current_3d_scene === nothing
        scene, lscene = Vis.scene()
    else
        scene = current_main_scene
        lscene = current_3d_scene
    end
    draw!(lscene, ob; kwargs...)
    display(scene)

    if (get_current_mode() == :pluto || get_current_mode() == :docs)
        return scene
    end
end

"""
    save(path::String)

Save the current Makie scene to an image file.
"""
function save(path::String)
    # save closes the window so just display it again as a work-around
    # for some reason the size isn't maintained automatically so we just reset it manually
    size = Makie.size(current_main_scene)
    Makie.save(path, current_3d_scene.scene)
    Makie.resize!(current_main_scene, size)
    display(current_main_scene)
end
function save(::Nothing) end

#############################################################################

function λtoRGB(λ::T, gamma::T = 0.8) where {T<:Real}
    wavelength = λ * 1000 # λ is in um, need in nm
    if (wavelength >= 380 && wavelength <= 440)
        attenuation = 0.3 + 0.7 * (wavelength - 380) / (440 - 380)
        R = ((-(wavelength - 440) / (440 - 380)) * attenuation)^gamma
        G = 0.0
        B = (1.0 * attenuation)^gamma
    elseif (wavelength >= 440 && wavelength <= 490)
        R = 0.0
        G = ((wavelength - 440) / (490 - 440))^gamma
        B = 1.0
    elseif (wavelength >= 490 && wavelength <= 510)
        R = 0.0
        G = 1.0
        B = (-(wavelength - 510) / (510 - 490))^gamma
    elseif (wavelength >= 510 && wavelength <= 580)
        R = ((wavelength - 510) / (580 - 510))^gamma
        G = 1.0
        B = 0.0
    elseif (wavelength >= 580 && wavelength <= 645)
        R = 1.0
        G = (-(wavelength - 645) / (645 - 580))^gamma
        B = 0.0
    elseif (wavelength >= 645 && wavelength <= 750)
        attenuation = 0.3 + 0.7 * (750 - wavelength) / (750 - 645)
        R = (1.0 * attenuation)^gamma
        G = 0.0
        B = 0.0
    else
        R = 0.0
        G = 0.0
        B = 0.0
    end
    return RGB(R, G, B)
end

indexedcolor(i::Int) = ColorSchemes.hsv[rem(i / (2.1 * π), 1.0)]
indexedcolor2(i::Int) = ColorSchemes.hsv[1.0 - rem(i / (2.1 * π), 1.0)] .* 0.5

#############################################################################

## displaying imported meshes

Base.:*(a::Transform, p::GeometryBasics.PointMeta) = a * p.main
Base.:*(a::Real, p::GeometryBasics.PointMeta{N,S}) where {S<:Real,N} = GeometryBasics.Point{N,S}((a * SVector{N,S}(p))...)
Base.:*(a::Transform, p::GeometryBasics.Point{N,S}) where {S<:Real,N} = GeometryBasics.Point{N,S}((a.rotation * SVector{N,S}(p) + a.translation)...)

function draw!(scene::Makie.LScene, ob::AbstractString; color = :gray, linewidth = 3, shaded::Bool = true, wireframe::Bool = false, transform::Transform{Float64} = identitytransform(Float64), scale::Float64 = 1.0, kwargs...)
    if any(endswith(lowercase(ob), x) for x in [".obj", "ply", ".2dm", ".off", ".stl"])
        meshdata = FileIO.load(ob)
        if transform != identitytransform(Float64) || scale != 1.0
            coords = [transform * (scale * p) for p in GeometryBasics.coordinates(meshdata)]
            meshdata = GeometryBasics.Mesh(coords, GeometryBasics.faces(meshdata))
        end
        Makie.mesh!(GeometryBasics.normal_mesh(meshdata); kwargs..., color = color, shading = shaded, visible = shaded)
        if wireframe
            if shaded
                Makie.wireframe!(scene[end][1], color = (:black, 0.1), linewidth = linewidth)
            else
                Makie.wireframe!(scene[end][1], color = color, linewidth = linewidth)
            end
        end
    else
        @error "Unsupported file type"
    end
end

## GEOMETRY

"""
    draw!(scene::Makie.LScene, surf::Surface{T}; numdivisions = 20, normals = false, normalcolor = :blue, kwargs...)

Transforms `surf` into a mesh using [`makemesh`](@ref) and draws the result.
`normals` of the surface can be drawn at evenly sampled points with provided `normalcolor`.
`numdivisions` determines the resolution with which the mesh is triangulated.
`kwargs` is passed on to the [`TriangleMesh`](@ref) drawing function.
"""
function draw!(scene::Makie.LScene, surf::Surface{T}; numdivisions::Int = 30, normals::Bool = false, normalcolor = :blue, kwargs...) where {T<:Real}
    mesh = makemesh(surf, numdivisions)
    if nothing === mesh
        return
    end
    draw!(scene, mesh; kwargs..., normals = false)
    if normals
        ndirs = Makie.Point3f0.(samplesurface(surf, normal, numdivisions ÷ 10))
        norigins = Makie.Point3f0.(samplesurface(surf, point, numdivisions ÷ 10))
        Makie.arrows!(scene, norigins, ndirs, arrowsize = 0.2, arrowcolor = normalcolor, linecolor = normalcolor, linewidth = 2)
    end
end

"""
    draw!(scene::Makie.LScene, tmesh::TriangleMesh{T}; linewidth = 3, shaded = true, wireframe = false, color = :orange, normals = false, normalcolor = :blue, transparency = false, kwargs...)

Draw a [`TriangleMesh`](@ref), optionially with a visible `wireframe`. `kwargs` are passed on to [`Makie.mesh`](http://makie.juliaplots.org/stable/plotting_functions.html#mesh).
"""
function draw!(scene::Makie.LScene, tmesh::TriangleMesh{T}; linewidth = 3, shaded::Bool = true, wireframe::Bool = false, color = :orange, normals::Bool = false, normalcolor = :blue, transparency::Bool = false, kwargs...) where {T<:Real}
    points, indices = makiemesh(tmesh)
    if length(points) > 0 && length(indices) > 0
        Makie.mesh!(scene, points, indices; kwargs..., color = color, shading = shaded, transparency = transparency, visible = shaded)
        if wireframe
            mesh = scene.scene[end][1]
            if shaded
                Makie.wireframe!(scene, mesh, color = (:black, 0.1), linewidth = linewidth)
            else
                Makie.wireframe!(scene, mesh, color = color, linewidth = linewidth)
            end
        end
    end
    if normals
        @warn "Normals being drawn from triangulated mesh, precision may be low"
        norigins = [Makie.Point3f0(centroid(t)) for t in tmesh.triangles[1:10:end]]
        ndirs = [Makie.Point3f0(normal(t)) for t in tmesh.triangles[1:10:end]]
        if length(norigins) > 0
            Makie.arrows!(scene, norigins, ndirs, arrowsize = 0.2, arrowcolor = normalcolor, linecolor = normalcolor, linewidth = 2)
        end
    end
end

"""
    draw!(scene::Makie.LScene, meshes::Vararg{S}; colors::Bool = false, kwargs...) where {T<:Real,S<:Union{TriangleMesh{T},Surface{T}}}

Draw a series of [`TriangleMesh`](@ref) or [`Surface`](@ref) objects, if `colors` is true then each mesh will be colored automatically with a diverse series of colors.
`kwargs` are is passed on to the drawing function for each element.
"""
function draw!(scene::Makie.LScene, meshes::Vararg{S}; colors::Bool = false, kwargs...) where {T<:Real,S<:Union{TriangleMesh{T},Surface{T}}}
    for i in 1:length(meshes)
        if colors
            col = indexedcolor2(i)
        else
            col = :orange
        end
        draw!(scene, meshes[i]; kwargs..., color = col)
    end
end

function draw(meshes::Vararg{S}; kwargs...) where {T<:Real,S<:Union{TriangleMesh{T},Surface{T}}}
    scene, lscene = Vis.scene()
    draw!(lscene, meshes...; kwargs...)
    Makie.display(scene)
    if (get_current_mode() == :pluto || get_current_mode() == :docs)
        return scene
    end
end

"""
    draw!(scene::Makie.LScene, csg::Union{CSGTree,CSGGenerator}; numdivisions::Int = 20, kwargs...)

Convert a CSG object ([`CSGTree`](@ref) or [`CSGGenerator`](@ref)) to a mesh using [`makemesh`](@ref) with resolution set by `numdivisions` and draw the resulting [`TriangleMesh`](@ref).
"""
draw!(scene::Makie.LScene, csg::CSGTree{T}; numdivisions::Int = 30, kwargs...) where {T<:Real} = draw!(scene, makemesh(csg, numdivisions); kwargs...)
draw!(scene::Makie.LScene, csg::CSGGenerator{T}; kwargs...) where {T<:Real} = draw!(scene, csg(); kwargs...)

"""
    draw!(scene::Makie.LScene, bbox::BoundingBox{T}; kwargs...)

Draw a [`BoundingBox`](@ref) as a wireframe, ie series of lines.
"""
function draw!(scene::Makie.LScene, bbox::BoundingBox{T}; kwargs...) where {T<:Real}
    p1 = SVector{3,T}(bbox.xmin, bbox.ymin, bbox.zmin)
    p2 = SVector{3,T}(bbox.xmin, bbox.ymax, bbox.zmin)
    p3 = SVector{3,T}(bbox.xmin, bbox.ymax, bbox.zmax)
    p4 = SVector{3,T}(bbox.xmin, bbox.ymin, bbox.zmax)
    p5 = SVector{3,T}(bbox.xmax, bbox.ymin, bbox.zmin)
    p6 = SVector{3,T}(bbox.xmax, bbox.ymax, bbox.zmin)
    p7 = SVector{3,T}(bbox.xmax, bbox.ymax, bbox.zmax)
    p8 = SVector{3,T}(bbox.xmax, bbox.ymin, bbox.zmax)
    Makie.linesegments!(scene, [p1, p2, p2, p3, p3, p4, p4, p1, p1, p5, p2, p6, p3, p7, p4, p8, p5, p6, p6, p7, p7, p8, p8, p5]; kwargs...)
end

## OPTICS

"""
    draw!(scene::Makie.LScene, ass::LensAssembly; kwargs...)

Draw each element in a [`LensAssembly`](@ref), with each element automatically colored differently.
"""
function draw!(scene::Makie.LScene, ass::LensAssembly{T}; kwargs...) where {T<:Real}
    for (i, e) in enumerate(elements(ass))
        draw!(scene, e; kwargs..., color = indexedcolor2(i))
    end
end

"""
    draw!(scene::Makie.LScene, sys::AbstractOpticalSystem; kwargs...)

Draw each element in the lens assembly of an [`AbstractOpticalSystem`](@ref), with each element automatically colored differently, as well as the detector of the system.
"""
function draw!(scene::Makie.LScene, sys::CSGOpticalSystem{T}; kwargs...) where {T<:Real}
    draw!(scene, sys.assembly; kwargs...)
    draw!(scene, sys.detector; kwargs...)
end

draw!(scene::Makie.LScene, sys::AxisymmetricOpticalSystem{T}; kwargs...) where {T<:Real} = draw!(scene, sys.system; kwargs...)

onlydetectorrays(system::Q, tracevalue::LensTrace{T,3}) where {T<:Real,Q<:AbstractOpticalSystem{T}} = onsurface(detector(system), point(tracevalue))

## RAY GEN
"""
    drawtracerays(system::Q; raygenerator::S = Source(transform = Transform.translation(0.0,0.0,10.0), origins = Origins.RectGrid(10.0,10.0,25,25),directions = Constant(0.0,0.0,-1.0)), test::Bool = false, trackallrays::Bool = false, colorbysourcenum::Bool = false, colorbynhits::Bool = false, rayfilter::Union{Nothing,Function} = onlydetectorrays, kwargs...)

Displays a model of the optical system.
`raygenerator` is an iterator that generates rays.
If `trackallrays` is true then ray paths from the emitter will be displayed otherwise just the final rays that intersect the image detector will be shown, not the entire ray path.
`colorbysourcenum` and `colorbynhits` will color rays accordingly, otherwise rays will be colored according to their wavelength.

By default only ray paths that eventually intersect the detector surface are displayed. If you want to display all ray paths set `rayfilter = nothing`.

Also `drawtracerays!` to add to an existing scene, with `drawsys` and `drawgen` to specify whether `system` and `raygenerator` should be drawn respectively.
"""
function drawtracerays(system::Q; raygenerator::S = Source(transform = translation(0.0,0.0,10.0), origins = Origins.RectGrid(10.0,10.0,25,25),directions = Constant(0.0,0.0,-1.0)), test::Bool = false, trackallrays::Bool = false, colorbysourcenum::Bool = false, colorbynhits::Bool = false, rayfilter::Union{Nothing,Function} = onlydetectorrays, verbose::Bool = false, resolution::Tuple{Int,Int} = (1000, 1000), kwargs...) where {T<:Real,Q<:AbstractOpticalSystem{T},S<:AbstractRayGenerator{T}}
    verbose && println("Drawing System...")
    s, ls = Vis.scene(resolution)

    drawtracerays!(ls, system, raygenerator = raygenerator, test = test, colorbysourcenum = colorbysourcenum, colorbynhits = colorbynhits, rayfilter = rayfilter, trackallrays = trackallrays, verbose = verbose, drawsys = true, drawgen = true; kwargs...)

    display(s)
    if (get_current_mode() == :pluto || get_current_mode() == :docs)
        return s
    end
end

drawtracerays!(system::Q; kwargs...) where {T<:Real,Q<:AbstractOpticalSystem{T}} = drawtracerays!(current_3d_scene, system; kwargs...)

function drawtracerays!(scene::Makie.LScene, system::Q; raygenerator::S = Source(transform = translation(0.0,0.0,10.0), origins = Origins.RectGrid(10.0,10.0,25,25),directions = Constant(0.0,0.0,-1.0)), test::Bool = false, trackallrays::Bool = false, colorbysourcenum::Bool = false, colorbynhits::Bool = false, rayfilter::Union{Nothing,Function} = onlydetectorrays, verbose::Bool = false, drawsys::Bool = false, drawgen::Bool = false, kwargs...) where {T<:Real,Q<:AbstractOpticalSystem{T},S<:AbstractRayGenerator{T}}
    raylines = Vector{LensTrace{T,3}}(undef, 0)

    drawgen && draw!(scene, raygenerator, norays = true; kwargs...)
    drawsys && draw!(scene, system; kwargs...)

    verbose && println("Tracing...")
    for (i, r) in enumerate(raygenerator)
        if i % 1000 == 0 && verbose
            print("\r $i / $(length(raygenerator))")
        end
        allrays = Vector{LensTrace{T,3}}(undef, 0)
        if trackallrays
            res = trace(system, r, trackrays = allrays, test = test)
        else
            res = trace(system, r, test = test)
        end

        if trackallrays && !isempty(allrays)
            if rayfilter === nothing || rayfilter(system, allrays[end])
                # filter on the trace that is hitting the detector
                for r in allrays
                    push!(raylines, r)
                end
            end
        elseif res !== nothing
            if rayfilter === nothing || rayfilter(system, res)
                push!(raylines, res)
            end
        end
    end
    verbose && print("\r")

    verbose && println("Drawing Rays...")
    draw!(scene, raylines, colorbysourcenum = colorbysourcenum, colorbynhits = colorbynhits; kwargs...)
end

"""
    drawtraceimage(system::Q; raygenerator::S = Source(transform = translation(0.0,0.0,10.0), origins = Origins.RectGrid(10.0,10.0,25,25),directions = Constant(0.0,0.0,-1.0)), test::Bool = false)

Traces rays from `raygenerator` through `system` and shows and returns the detector image. `verbose` will print progress updates.
"""
function drawtraceimage(system::Q; raygenerator::S = Source(transform = translation(0.0,0.0,10.0), origins = Origins.RectGrid(10.0,10.0,25,25),directions = Constant(0.0,0.0,-1.0)), test::Bool = false, verbose::Bool = false) where {T<:Real,Q<:AbstractOpticalSystem{T},S<:AbstractRayGenerator{T}}
    resetdetector!(system)
    start_time = time()
    for (i, r) in enumerate(raygenerator)
        if i % 1000 == 0 && verbose
            dif = round(time() - start_time, digits = 1)
            left = round((time() - start_time) / i * (length(raygenerator) - i), digits = 1)
            print("\rTraced: $i / $(length(raygenerator))     Elapsed: $(dif)s     Left: $(left)s")
        end
        res = trace(system, r, test = test)
    end
    verbose && print("\r")
    tot = round(time() - start_time, digits = 1)
    verbose && println("Completed in $(tot)s")
    show(detectorimage(system))
    return detectorimage(system)
end

## RAYS, LINES AND POINTS

"""
    draw!(scene::Makie.LScene, rays::AbstractVector{<:AbstractRay{T,N}}; kwargs...)

Draw a vector of [`Ray`](@ref) or [`OpticalRay`](@ref) objects.
"""
function draw!(scene::Makie.LScene, rays::AbstractVector{<:AbstractRay{T,N}}; kwargs...) where {T<:Real,N}
    for r in rays
        draw!(scene, r; kwargs...)
    end
end

"""
    draw!(scene::Makie.LScene, traces::AbstractVector{LensTrace{T,N}}; kwargs...)

Draw a vector of [`LensTrace`](@ref) objects.
"""
function draw!(scene::Makie.LScene, traces::AbstractVector{LensTrace{T,N}}; kwargs...) where {T<:Real,N}
    for t in traces
        draw!(scene, t; kwargs...)
    end
end

"""
    draw!(scene::Makie.LScene, trace::LensTrace{T,N}; colorbysourcenum::Bool = false, colorbynhits::Bool = false, kwargs...)

Draw a [`LensTrace`](@ref) as a line which can be colored automatically by its `sourcenum` or `nhits` attributes.
The alpha is determined by the `power` attribute of `trace`.
"""
function draw!(scene::Makie.LScene, trace::LensTrace{T,N}; colorbysourcenum::Bool = false, colorbynhits::Bool = false, kwargs...) where {T<:Real,N}
    if colorbysourcenum
        color = indexedcolor(sourcenum(trace))
    elseif colorbynhits
        color = indexedcolor(nhits(trace))
    else
        color = λtoRGB(wavelength(trace))
    end
    draw!(scene, (origin(ray(trace)), point(intersection(trace))); kwargs..., color = RGBA(color.r, color.g, color.b, sqrt(power(trace))), transparency = true)
end

"""
    draw!(scene::Makie.LScene, ray::OpticalRay{T,N}; colorbysourcenum::Bool = false, colorbynhits::Bool = false, kwargs...)

Draw an [`OpticalRay`](@ref) which can be colored automatically by its `sourcenum` or `nhits` attributes.
The alpha of the ray is determined by the `power` attribute of `ray`.
`kwargs` are passed to `draw!(scene, ray::Ray)`.
"""
function draw!(scene::Makie.LScene, r::OpticalRay{T,N}; colorbysourcenum::Bool = false, colorbynhits::Bool = false, kwargs...) where {T<:Real,N}
    if colorbysourcenum
        color = indexedcolor(sourcenum(r))
    elseif colorbynhits
        color = indexedcolor(nhits(r))
    else
        color = λtoRGB(wavelength(r))
    end
    draw!(scene, ray(r); kwargs..., color = RGBA(color.r, color.g, color.b, sqrt(power(r))), transparency = true, rayscale = power(r))
end

"""
    draw!(scene::Makie.LScene, ray::Ray{T,N}; color = :yellow, rayscale = 1.0, kwargs...)

Draw a [`Ray`](@ref) in a given `color` optionally scaling the size using `rayscale`.
`kwargs` are passed to [`Makie.arrows`](http://makie.juliaplots.org/stable/plotting_functions.html#arrows).
"""
function draw!(scene::Makie.LScene, ray::AbstractRay{T,N}; color = :yellow, rayscale = 1.0, kwargs...) where {T<:Real,N}
    arrow_size = min(0.05, rayscale * 0.05)
    Makie.arrows!(scene, [Makie.Point3f0(origin(ray))], [Makie.Point3f0(rayscale * direction(ray))]; kwargs..., arrowsize = arrow_size, arrowcolor = color, linecolor = color, linewidth=arrow_size * 0.5)
end

"""
    draw!(scene::Makie.LScene, du::DisjointUnion{T}; kwargs...)

Draw each [`Interval`](@ref) in a [`DisjointUnion`](@ref).
"""
function draw!(scene::Makie.LScene, du::DisjointUnion{T}; kwargs...) where {T<:Real}
    draw!(scene, intervals(du); kwargs...)
end

"""
    draw!(scene::Makie.LScene, intervals::AbstractVector{Interval{T}}; kwargs...)

Draw a vector of [`Interval`](@ref)s.
"""
function draw!(scene::Makie.LScene, intervals::AbstractVector{Interval{T}}; kwargs...) where {T<:Real}
    for i in intervals
        draw!(scene, i; kwargs...)
    end
end

"""
    draw!(scene::Makie.LScene, interval::Interval{T}; kwargs...)

Draw an [`Interval`](@ref) as a line with circles at each [`Intersection`](@ref) point.
"""
function draw!(scene::Makie.LScene, interval::Interval{T}; kwargs...) where {T<:Real}
    if !(interval isa EmptyInterval)
        l = lower(interval)
        u = upper(interval)
        if !(l isa RayOrigin)
            draw!(scene, l)
        else
            @warn "Negative half space, can't draw ray origin"
        end
        if !(u isa Infinity)
            draw!(scene, u)
        else
            @warn "Positive half space, can't draw end point"
        end
        if !(l isa RayOrigin) && !(u isa Infinity)
            draw!(scene, (point(l), point(u)); kwargs...)
        end
    end
end

"""
    draw!(scene::Makie.LScene, intersection::Intersection; normal::Bool = false, kwargs...)

Draw an [`Intersection`](@ref) as a circle, optionally showing the surface normal at the point.
"""
function draw!(scene::Makie.LScene, intersection::Intersection; normal::Bool = false, kwargs...)
    draw!(scene, point(intersection); kwargs...)
    if normal
        draw!(scene, Ray(point(intersection), normal(intersection)); kwargs...)
    end
end

"""
    draw!(scene::Makie.LScene, lines::AbstractVector{Tuple{AbstractVector{T},AbstractVector{T}}}; kwargs...)

Draw a vector of lines.
"""
function draw!(scene::Makie.LScene, lines::AbstractVector{Tuple{P,P}}; kwargs...) where {T<:Real,P<:AbstractVector{T}}
    for l in lines
        draw!(scene, l; kwargs...)
    end
end

"""
    draw!(scene::Makie.LScene, line::Tuple{AbstractVector{T},AbstractVector{T}}; color = :yellow, kwargs...)

Draw a line between two points, `kwargs` are passed to [`Makie.linesegments`](http://makie.juliaplots.org/stable/plotting_functions.html#linesegments).
"""
function draw!(scene::Makie.LScene, line::Tuple{P,P}; color = :yellow, kwargs...) where {T<:Real,P<:AbstractVector{T}}
    Makie.linesegments!(scene, [line[1], line[2]]; kwargs..., color = color)
end

"""
    draw!(s::Makie.LScene, point::AbstractVector{T}; kwargs...)

Draw a single point, `kwargs` are passed to `draw!(scene, points::AbstractVector{AbstractVector{T}})`.
"""
function draw!(s::Makie.LScene, point::AbstractVector{T}; kwargs...) where {T<:Real}
    draw!(s, [point]; kwargs...)
end

"""
    draw!(scene::Makie.LScene, points::AbstractVector{AbstractVector{T}}; markersize = 20, color = :black, kwargs...)

Draw a vector of points.
`kwargs` are passed to [`Makie.scatter`](http://makie.juliaplots.org/stable/plotting_functions.html#scatter).
"""
function draw!(scene::Makie.LScene, points::AbstractVector{P}; markersize = 20, color = :black, kwargs...) where {T<:Real,P<:AbstractVector{T}}
    Makie.scatter!(scene, points, markersize = markersize, color = color, strokewidth = 0; kwargs...)
end

#######################################################

function plotOPD!(sys::AxisymmetricOpticalSystem{T}; label = nothing, color = nothing, collimated::Bool = true, axis::Int = 1, samples::Int = 5, wavelength::T = 0.55, sourcepos::SVector{3,T} = SVector{3,T}(0.0, 0.0, 10.0), kwargs...) where {T<:Real}
    sysdiam = semidiameter(sys)
    rays = Vector{OpticalRay{T,3}}(undef, 0)
    positions = Vector{T}(undef, 0)
    if axis == 2
        if collimated
            dir = -sourcepos
            for v in 0:samples
                s = (2 * (v / samples) - 1.0) * sysdiam
                p = SVector{3,T}(0, s, 0)
                push!(positions, s)
                push!(rays, OpticalRay(p - dir, dir, 1.0, wavelength))
            end
        else
            for v in 0:samples
                s = (2 * (v / samples) - 1.0) * sysdiam
                p = SVector{3,T}(0, s, 0)
                dir = p - sourcepos
                push!(positions, s)
                push!(rays, OpticalRay(sourcepos, dir, 1.0, wavelength))
            end
        end
    elseif axis == 1
        if collimated
            dir = -sourcepos
            for v in 0:samples
                s = (2 * (v / samples) - 1.0) * sysdiam
                p = SVector{3,T}(s, 0, 0)
                push!(positions, s)
                push!(rays, OpticalRay(p - dir, dir, 1.0, wavelength))
            end
        else
            for v in 0:samples
                s = (2 * (v / samples) - 1.0) * sysdiam
                p = SVector{3,T}(s, 0, 0)
                dir = p - sourcepos
                push!(positions, s)
                push!(rays, OpticalRay(sourcepos, dir, 1.0, wavelength))
            end
        end
    end
    raygenerator = RayListSource(rays)
    chiefray = OpticalRay(sourcepos, -sourcepos, 1.0, wavelength)
    chiefres = trace(sys, chiefray, test = true)
    xs = Vector{T}()
    opls = Vector{T}()
    for (i, r) in enumerate(raygenerator)
        lt = trace(sys, r, test = true)
        if !(nothing === lt)
            push!(xs, positions[i + 1])
            push!(opls, (pathlength(lt) - pathlength(chiefres)) / (wavelength / 1000)) # wavelength is um while OPD is mm
        end
    end
    p = Plots.plot!(xs, opls, color = nothing === color ? λtoRGB(wavelength) : color, label = label)
    Plots.xlabel!("Relative position on entrance pupil")
    Plots.ylabel!("OPD (waves)")
    return p
end

"""
    spotdiag(sys::CSGOpticalSystem{T}, raygenerator::OpticalRayGenerator{T}; size = (500, 500), kwargs...)

Plot a spot diagram for an arbitrary [`CSGOpticalSystem`](@ref) and [`OpticalRayGenerator`](@ref).
All rays from `raygenerator` will be traced through `sys` and their intersection location on the detector plotted.

Also `spotdiag!` of the same arguments to add to an existing plot.
"""
function spotdiag(sys::CSGOpticalSystem{T}, raygenerator::OpticalRayGenerator{T}; kwargs...) where {T<:Real}
    Plots.plot()
    return spotdiag!(sys, raygenerator; kwargs...)
end

"""
    spotdiag(sys::AxisymmetricOpticalSystem{T}; size = (500, 500), hexapolar::Bool = true, collimated::Bool = true, samples::Int = 5, wavelength::T = 0.55, sourceangle::T = zero(T), sourcepos::SVector{3,T} = SVector{3,T}(0.0, 0.0, 10.0), kwargs...)

Plot a spot diagram for an [`AxisymmetricOpticalSystem`](@ref), rays are distributed across the entrance pupil of the system either in a hexapolar of rectangular grid pattern depending on `hexapolar`.

The input rays can be `collimated`, in which case the `sourceangle` parameter determines their direction.
Otherwise rays are treated as coming from a point source at `sourcepos`.

Also `spotdiag!` of the same arguments to add to an existing plot.
"""
function spotdiag(sys::AxisymmetricOpticalSystem{T}; kwargs...) where {T<:Real}
    Plots.plot()
    return spotdiag!(sys; kwargs...)
end

function spotdiag!(sys::AxisymmetricOpticalSystem{T}; color = nothing, hexapolar::Bool = true, collimated::Bool = true, samples::Int = 5, wavelength::T = 0.55, sourceangle::T = zero(T), sourcepos::SVector{3,T} = SVector{3,T}(0.0, 0.0, 10.0), kwargs...) where {T<:Real}
    if hexapolar
        source = HexapolarField(sys, samples = samples, collimated = collimated, wavelength = wavelength, sourcepos = sourcepos, sourceangle = sourceangle)
    else
        source = GridField(sys, samples = samples, collimated = collimated, wavelength = wavelength, sourcepos = sourcepos, sourceangle = sourceangle)
    end
    return spotdiag!(sys.system, source, color = color; kwargs...)
end

function spotdiag!(sys::CSGOpticalSystem{T}, raygenerator::OpticalRayGenerator{T}; color = nothing, label = nothing, size = (500, 500), kwargs...) where {T<:Real}
    xs = Vector{T}()
    ys = Vector{T}()
    det = detector(sys)
    for (i, r) in enumerate(raygenerator)
        lt = trace(sys, r, test = true)
        if !(nothing === lt)
            push!(xs, dot(point(lt) - centroid(det), det.uvec))
            push!(ys, dot(point(lt) - centroid(det), det.vvec))
        end
    end
    # relative to centroid
    cen = sum(xs) / length(xs), sum(ys) / length(ys)
    xs .-= cen[1]
    ys .-= cen[2]
    p = Plots.scatter!(xs, ys, color = (color === nothing) ? λtoRGB(wavelength(raygenerator[0])) : color, label = label, marker = :+, markerstrokewidth = 0, markerstrokealpha = 0, size = size, aspect_ratio = :equal)
    Plots.xlabel!(p, "x (mm)")
    Plots.ylabel!(p, "y (mm)")
    return p
end

spotdiaggrid(sys::AxisymmetricOpticalSystem{T}, raygenerators::Vector{<:OpticalRayGenerator{T}}; kwargs...) where {T<:Real} = spotdiaggrid(sys.system, raygenerators; kwargs...)

function spotdiaggrid(sys::CSGOpticalSystem{T}, raygenerators::Vector{<:OpticalRayGenerator{T}}; colorbynum::Bool = false, kwargs...) where {T<:Real}
    plots = []
    x = Int(floor(sqrt(length(raygenerators))))
    y = length(raygenerators) ÷ x
    size = 2000 / max(x, y)
    xmin = 0
    xmax = 0
    ymin = 0
    ymax = 0
    for (i, r) in enumerate(raygenerators)
        if colorbynum
            p = spotdiag(sys, r; color = indexedcolor2(i), kwargs..., size = (size, size))
        else
            p = spotdiag(sys, r; kwargs..., size = (size, size))
        end
        xmin = min(min(Plots.xlims(p)...), xmin)
        xmax = max(max(Plots.xlims(p)...), xmax)
        ymin = min(min(Plots.ylims(p)...), ymin)
        ymax = max(max(Plots.ylims(p)...), ymax)
        push!(plots, p)
    end
    # set all plots to the same range
    for p in plots
        Plots.xlims!(p, xmin, xmax)
        Plots.ylims!(p, ymin, ymax)
    end
    return Plots.plot(plots..., layout = Plots.grid(x, y), size = (y * size, x * size))
end

"""
    surfacesag(object::Union{CSGTree{T},Surface{T}}, resolution::Tuple{Int,Int}, halfsizes::Tuple{T,T}; offset::T = T(10), position::SVector{3,T} = SVector{3,T}(0.0, 0.0, 10.0), direction::SVector{3,T} = SVector{3,T}(0.0, 0.0, -1.0), rotationvec::SVector{3,T} = SVector{3,T}(0.0, 1.0, 0.0))

Calculates and displays the surface sag of an arbitrary [`Surface`](@ref) or [`CSGTree`](@ref).

Rays are shot in a grid of size defined by `resolution` across a arectangular area defined by `halfsizes`.
This rectangle is centred at `postion` with normal along `direction` and rotation defined by `rotationvec`.
`offset` is subtracted from the sag measurements to provide values relative to the appropriate zero level.
"""
function surfacesag(object::Union{CSGTree{T},Surface{T}}, resolution::Tuple{Int,Int}, halfsizes::Tuple{T,T}; offset::T = T(10), position::SVector{3,T} = SVector{3,T}(0.0, 0.0, 10.0), direction::SVector{3,T} = SVector{3,T}(0.0, 0.0, -1.0), rotationvec::SVector{3,T} = SVector{3,T}(0.0, 1.0, 0.0)) where {T<:Real}
    direction = normalize(direction)
    if abs(dot(rotationvec, direction)) == one(T)
        rotationvec = SVector{3,T}(1.0, 0.0, 0.0)
    end
    uvec = normalize(cross(normalize(rotationvec), direction))
    vvec = normalize(cross(direction, uvec))
    image = Array{T,2}(undef, resolution)
    image .= NaN
    Threads.@threads for x in 0:(resolution[1] - 1)
        for y in 0:(resolution[2] - 1)
            u = 2 * (x / (resolution[1] - 1)) - 1.0
            v = 2 * (y / (resolution[2] - 1)) - 1.0
            pos = position + (halfsizes[1] * u * uvec) + (halfsizes[2] * v * vvec)
            r = Ray(pos, direction)
            int = closestintersection(surfaceintersection(object, r), false)
            if int !== nothing
                int = int::Intersection{T,3}
                sag = α(int) - offset
                image[y + 1, x + 1] = -sag
            end
        end
    end
    p = Plots.heatmap((-halfsizes[1]):(2 * halfsizes[1] / resolution[1]):halfsizes[1], (-halfsizes[2]):(2 * halfsizes[2] / resolution[2]):halfsizes[2], image, c = :jet1, aspect_ratio = :equal, size = (1000, 1000))
    Plots.xlims!(p, -halfsizes[1], halfsizes[1])
    Plots.ylims!(p, -halfsizes[2], halfsizes[2])
    return p
end

"""
    eyebox_eval_eye(assembly::LensAssembly{T}, raygen::OpticalRayGenerator{T}, eye_rotation_x::T, eye_rotation_y::T, sample_points_x::Int, sample_points_y::Int; pupil_radius::T = T(2.0), resolution::Int = 512, eye_transform::Transform{T} = identitytransform(T))

Visualise the images formed when tracing `assembly` with a human eye for an evenly sampled `sample_points_x × sample_points_y` grid in the angular range of eyeball rotations `-eye_rotation_x:eye_rotation_x` and `-eye_rotation_y:eye_rotation_y` in each dimension respectively.
`resolution` is the size of the detector image (necessarily square).

The eye must be positioned appropriately relative to the system using `eye_transform`, this should transform the eye to the correct position and orientation when at 0 rotation in both dimensions.
By default the eye is directed along the positive z-axis with the vertex of the cornea at the origin.

The result is displayed as a 4D image - the image seen by the eye is shown in 2D as normal with sliders to vary eye rotation in x and y.
The idea being that the whole image should be visible for all rotations in the range.
"""
function eyebox_eval_eye(assembly::LensAssembly{T}, raygen::OpticalRayGenerator{T}, eye_rotation_x::T, eye_rotation_y::T, sample_points_x::Int, sample_points_y::Int; pupil_radius::T = T(2.0), resolution::Int = 512, eye_transform::Transform{T} = identitytransform(T)) where {T<:Real}
    out_image = zeros(Float32, resolution, resolution, sample_points_y, sample_points_x)
    xrange = sample_points_x > 1 ? range(-eye_rotation_x, stop = eye_rotation_x, length = sample_points_x) : [0.0]
    yrange = sample_points_y > 1 ? range(-eye_rotation_y, stop = eye_rotation_y, length = sample_points_y) : [0.0]

    for (xi, erx) in enumerate(xrange)
        for (yi, ery) in enumerate(yrange)
            print("Position: ($xi, $yi) / ($sample_points_x, $sample_points_y)\r")
            r = OpticSim.rotmatd(T, ery, erx, 0.0)
            ov = OpticSim.rotate(eye_transform, SVector{3,T}(0.0, 0.0, 13.0))
            t = r * ov - ov
            sys = ModelEye(assembly, pupil_radius = pupil_radius, detpixels = resolution, transform = eye_transform * Transform(r, t))
            # Vis.draw(sys)
            # Vis.draw!((eye_transform.translation - ov - SVector(0.0, 20.0, 0.0), eye_transform.translation - ov + SVector(0.0, 20.0, 0.0)))
            # Vis.draw!((eye_transform.translation - ov - SVector(20.0, 0.0, 0.0), eye_transform.translation - ov + SVector(20.0, 0.0, 0.0)))
            # Vis.draw!(eye_transform.translation, markersize = 500.0)
            out_image[:, :, yi, xi] = traceMT(sys, raygen, printprog = false)
        end
    end
    # this is a 4D image - ImageView makes sliders to visualise the dimensions >2!
    Vis.show(out_image)
end

"""
    eyebox_eval_planar(assembly::LensAssembly{T}, raygen::OpticalRayGenerator{T}, eyebox::Rectangle{T}, sample_points_x::Int, sample_points_y::Int, vsize::T; pupil_radius::T = T(2.0), resolution::Int = 512)

Visualise the images formed when tracing `assembly` for multiple pupil positions within a planar `eyebox`.
Any angles which are present in a circle radius `pupil_radius` around each sampling point on an even `sample_points_x × sample_points_y` grid on the `eyebox` are added to the sub-image at that grid point.

A paraxial lens focal length 1mm is placed at the eyebox and a detector of size `vsize × vsize`mm placed 1mm behind it.
The normal of the detector rectangle should point towards the system (and away from the fake detector).
Any rays which miss the detector are ignored.

The pupil is always fully contained in the eyebox, i.e., the extreme sample position in u would be `eyebox.halfsizeu - pupil_radius`, for example.

The result is displayed as a 4D image - each sub-image is shown as normal with sliders to vary eye rotation in x and y.
The idea being that the whole FoV should be visible for all rotations in the range.
"""
function eyebox_eval_planar(assembly::LensAssembly{T}, raygen::OpticalRayGenerator{T}, eyebox::Rectangle{T}, sample_points_x::Int, sample_points_y::Int, vsize::T; pupil_radius::T = T(2.0), resolution::Int = 512) where {T<:Real}
    sys = CSGOpticalSystem(assembly, eyebox, 1, 1)
    hits = tracehitsMT(sys, raygen)
    println("Calculating results...")
    det = Rectangle(vsize, vsize, normal(eyebox), centroid(eyebox) - normal(eyebox))
    if length(hits) > 0
        out_image = zeros(Float32, resolution, resolution, sample_points_y, sample_points_x)
        for res in hits
            dir = direction(ray(res))
            u, v = uv(res)
            x = u * eyebox.halfsizeu
            y = v * eyebox.halfsizev

            raydirection = normalize((centroid(eyebox) - point(res)) + direction(ray(res)) / dot(-normal(eyebox), direction(ray(res))))
            focal_plane_hit = surfaceintersection(det, Ray(point(res), raydirection))
            if focal_plane_hit isa EmptyInterval
                continue
            end
            px, py = OpticSim.uvtopix(det, uv(OpticSim.halfspaceintersection(focal_plane_hit)), (resolution, resolution))

            hsu = eyebox.halfsizeu - pupil_radius
            hsv = eyebox.halfsizev - pupil_radius
            xrange = sample_points_x > 1 ? range(-hsu, stop = hsu, length = sample_points_x) : [0.0]
            yrange = sample_points_y > 1 ? range(-hsv, stop = hsv, length = sample_points_y) : [0.0]

            for (xi, posx) in enumerate(xrange)
                for (yi, posy) in enumerate(yrange)
                    if (x - posx)^2 + (y - posy)^2 < pupil_radius^2
                        out_image[py, px, yi, xi] += sourcepower(ray(res))
                    end
                end
            end
        end
        # this is a 4D image - ImageView makes sliders to visualise the dimensions >2!
        Vis.show(out_image)
    end
end
