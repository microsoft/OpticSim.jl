# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

using ..OpticSim, .Geometry
using LinearAlgebra
using Distributions
using StaticArrays

import Makie

using .Emitters
using .Emitters.Spectrum
using .Emitters.Directions
using .Emitters.Origins
using .Emitters.AngularPower
using .Emitters.Sources

const ARRROW_LENGTH = 0.5
const ARRROW_SIZE = 0.01
const MARKER_SIZE = 1


#-------------------------------------
# draw debug information - local axes and positions
#-------------------------------------
function maybe_draw_debug_info(scene::Makie.LScene, o::Origins.AbstractOriginDistribution; transform::Geometry.Transform = Transform(), debug::Bool=false, kwargs...) where {T<:Real}

    dir = forward(transform)
    uv = SVector{3}(right(transform))
    vv = SVector{3}(up(transform))
    pos = origin(transform)

    if (debug)
        # this is a stupid hack to force makie to render in 3d - for some scenes, makie decide with no apparent reason to show in 2d instead of 3d
        Makie.scatter!(scene, [pos[1], pos[1]+0.1], [pos[2], pos[2]+0.1], [pos[3], pos[3]+0.1], color=:red, markersize=0)

        # draw the origin and normal of the surface
        Makie.scatter!(scene, pos, color=:blue, markersize = MARKER_SIZE * visual_size(o))

        # normal
        arrow_size = ARRROW_SIZE * visual_size(o)
        arrow_start = pos
        arrow_end = dir * ARRROW_LENGTH * visual_size(o) 
        Makie.arrows!(scene.scene, [Makie.Point3f(arrow_start)], [Makie.Point3f(arrow_end)], arrowsize=arrow_size, linewidth=arrow_size * 0.5, linecolor=:blue, arrowcolor=:blue)
        arrow_end = uv * 0.5 * ARRROW_LENGTH * visual_size(o) 
        Makie.arrows!(scene.scene, [Makie.Point3f(arrow_start)], [Makie.Point3f(arrow_end)], arrowsize= 0.5 * arrow_size, linewidth=arrow_size * 0.5, linecolor=:red, arrowcolor=:red)
        arrow_end = vv * 0.5 * ARRROW_LENGTH * visual_size(o) 
        Makie.arrows!(scene.scene, [Makie.Point3f(arrow_start)], [Makie.Point3f(arrow_end)], arrowsize= 0.5 * arrow_size, linewidth=arrow_size * 0.5, linecolor=:green, arrowcolor=:green)

        # draw all the samples origins
        positions = map(x -> transform*x, collect(o))
        positions = collect(Makie.Point3f, positions)
        Makie.scatter!(scene, positions, color=:green, markersize = MARKER_SIZE * visual_size(o))

        # positions = collect(Makie.Point3f, o)
        # Makie.scatter!(scene, positions, color=:green, markersize = MARKER_SIZE * visual_size(o))
    end

end


#-------------------------------------
# draw point origin
#-------------------------------------
function OpticSim.Vis.draw!(scene::Makie.LScene, o::Origins.Point; transform::Geometry.Transform = Transform(), kwargs...) where {T<:Real}
        maybe_draw_debug_info(scene, o; transform=transform, kwargs...)
end

#-------------------------------------
# draw RectGrid and RectUniform origins
#-------------------------------------
function OpticSim.Vis.draw!(scene::Makie.LScene, o::Union{Origins.RectGrid, Origins.RectUniform}; transform::Geometry.Transform = Transform(), kwargs...) where {T<:Real}
    dir = forward(transform)
    uv = SVector{3}(right(transform))
    vv = SVector{3}(up(transform))
    pos = origin(transform)

    # @info "RECT: transform $(pos)"

    plane = OpticSim.Plane(dir, pos)
    rect = OpticSim.Rectangle(plane, o.width / 2, o.height / 2, uv, vv)
    
    OpticSim.Vis.draw!(scene, rect;  kwargs...)

    maybe_draw_debug_info(scene, o; transform=transform, kwargs...)
end


#-------------------------------------
# draw hexapolar origin
#-------------------------------------
function OpticSim.Vis.draw!(scene::Makie.LScene, o::Origins.Hexapolar; transform::Geometry.Transform = Transform(), kwargs...) where {T<:Real}
    dir = forward(transform)
    uv = SVector{3}(right(transform))
    vv = SVector{3}(up(transform))
    pos = origin(transform)

    plane = OpticSim.Plane(dir, pos)
    ellipse = OpticSim.Ellipse(plane, o.halfsizeu, o.halfsizev, uv, vv)
    
    OpticSim.Vis.draw!(scene, ellipse;  kwargs...)

    maybe_draw_debug_info(scene, o; transform=transform, kwargs...)
end

#-------------------------------------
# draw source
#-------------------------------------
function OpticSim.Vis.draw!(scene::Makie.LScene, s::S; parent_transform::Geometry.Transform = Transform(), debug::Bool=false, kwargs...) where {T<:Real,S<:Sources.AbstractSource{T}}
   
    OpticSim.Vis.draw!(scene, Emitters.Sources.origins(s);  transform=parent_transform * Emitters.Sources.transform(s), debug=debug, kwargs...)

    if (debug)
        m = zeros(T, length(s), 7)
        for (index, optical_ray) in enumerate(s)
            ray = OpticSim.ray(optical_ray)
            ray = parent_transform * ray
            m[index, 1:7] = [ray.origin... ray.direction... OpticSim.power(optical_ray)]
        end
        
        m[:, 4:6] .*= m[:, 7] * ARRROW_LENGTH * visual_size(Emitters.Sources.origins(s))  

        # Makie.arrows!(scene, [Makie.Point3f(origin(ray))], [Makie.Point3f(rayscale * direction(ray))]; kwargs..., arrowsize = min(0.05, rayscale * 0.05), arrowcolor = color, linecolor = color, linewidth = 2)
        color = :yellow
        arrow_size = ARRROW_SIZE * visual_size(Emitters.Sources.origins(s))
        Makie.arrows!(scene, m[:,1], m[:,2], m[:,3], m[:,4], m[:,5], m[:,6]; kwargs...,  arrowcolor=color, linecolor=color, arrowsize=arrow_size, linewidth=arrow_size*0.5)
    end

    # for ray in o
    #     OpticSim.Vis.draw!(scene, ray)
    # end
end

#-------------------------------------
# draw optical rays
#-------------------------------------
function OpticSim.Vis.draw!(scene::Makie.LScene, rays::AbstractVector{OpticSim.OpticalRay{T, 3}}; kwargs...) where {T<:Real}
    m = zeros(T, length(rays)*2, 3)
    for (index, optical_ray) in enumerate(rays)
        ray = OpticSim.ray(optical_ray)
        m[(index-1)*2+1, 1:3] = [origin(optical_ray)...]
        m[(index-1)*2+2, 1:3] = [(OpticSim.origin(optical_ray) + OpticSim.direction(optical_ray) * 1 * OpticSim.power(optical_ray))... ]
    end
    
    color = :green
    Makie.linesegments!(scene, m[:,1], m[:,2], m[:,3]; kwargs...,  color = color, linewidth = 2, )
end

#-------------------------------------
# draw composite source
#-------------------------------------
function OpticSim.Vis.draw!(scene::Makie.LScene, s::Sources.CompositeSource{T}; parent_transform::Geometry.Transform = Transform(), kwargs...) where {T<:Real}
    for source in s.sources
        OpticSim.Vis.draw!(scene, source; parent_transform=parent_transform*Emitters.Sources.transform(s), kwargs...)
    end
end
