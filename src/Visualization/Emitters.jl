import ..OpticSim
using ..OpticSim.Emitters
using ..OpticSim.Emitters.Geometry
using ..OpticSim.Emitters.Spectrum
using ..OpticSim.Emitters.Directions
using ..OpticSim.Emitters.Origins
using ..OpticSim.Emitters.AngularPower
using ..OpticSim.Emitters.Sources

using LinearAlgebra
using Distributions
using StaticArrays

import Makie
import Makie.AbstractPlotting
import Makie.AbstractPlotting.MakieLayout

const ARRROW_LENGTH = 0.5
const ARRROW_SIZE = 0.01
const MARKER_SIZE = 15


#-------------------------------------
# draw debug information - local axes and positions
#-------------------------------------
function maybe_draw_debug_info(scene::MakieLayout.LScene, o::Origins.AbstractOriginDistribution; transform::Geometry.Transform = Transform(), debug::Bool=false, kwargs...) where {T<:Real}

    dir = forward(transform)
    uv = SVector{3}(right(transform))
    vv = SVector{3}(up(transform))
    pos = OpticSim.origin(transform)

    if (debug)
        # this is a stupid hack to force makie to render in 3d - for some scenes, makie decide with no apperent reason to show in 2d instead of 3d
        AbstractPlotting.scatter!(scene, [pos[1], pos[1]+0.1], [pos[2], pos[2]+0.1], [pos[3], pos[3]+0.1], color=:red, markersize=0)

        # draw the origin and normal of the surface
        Makie.scatter!(scene, pos, color=:blue, markersize = MARKER_SIZE * visual_size(o))

        # normal
        arrow_start = pos
        arrow_end = dir * ARRROW_LENGTH * visual_size(o) 
        Makie.arrows!(scene.scene, [AbstractPlotting.Point3f0(arrow_start)], [AbstractPlotting.Point3f0(arrow_end)], arrowsize=ARRROW_SIZE * visual_size(o), arrowcolor=:blue)
        arrow_end = uv * 0.5 * ARRROW_LENGTH * visual_size(o) 
        Makie.arrows!(scene.scene, [AbstractPlotting.Point3f0(arrow_start)], [AbstractPlotting.Point3f0(arrow_end)], arrowsize= 0.5 * ARRROW_SIZE * visual_size(o), arrowcolor=:red)
        arrow_end = vv * 0.5 * ARRROW_LENGTH * visual_size(o) 
        Makie.arrows!(scene.scene, [AbstractPlotting.Point3f0(arrow_start)], [AbstractPlotting.Point3f0(arrow_end)], arrowsize= 0.5 * ARRROW_SIZE * visual_size(o), arrowcolor=:green)

        # draw all the samples origins
        positions = map(x -> local2world(transform, x), collect(o))
        positions = collect(AbstractPlotting.Point3f0, positions)
        Makie.scatter!(scene, positions, color=:green, markersize = MARKER_SIZE * visual_size(o))

        # positions = collect(AbstractPlotting.Point3f0, o)
        # Makie.scatter!(scene, positions, color=:green, markersize = MARKER_SIZE * visual_size(o))
    end

end


#-------------------------------------
# draw point origin
#-------------------------------------
# function OpticSim.Vis.draw!(scene::MakieLayout.LScene, o::Origins.Point; transform::Geometry.Transform = Transform(), kwargs...) where {T<:Real}
function OpticSim.Vis.draw!(scene::MakieLayout.LScene, o::Origins.Point; transform::Geometry.Transform = Transform(), kwargs...) where {T<:Real}
        maybe_draw_debug_info(scene, o; transform=transform, kwargs...)
end

#-------------------------------------
# draw RectGrid and RectUniform origins
#-------------------------------------
function OpticSim.Vis.draw!(scene::MakieLayout.LScene, o::Union{Origins.RectGrid, Origins.RectUniform}; transform::Geometry.Transform = Transform(), kwargs...) where {T<:Real}
    dir = forward(transform)
    uv = SVector{3}(right(transform))
    vv = SVector{3}(up(transform))
    pos = OpticSim.origin(transform)

    # @info "RECT: transform $(pos)"

    plane = OpticSim.Plane(dir, pos)
    rect = OpticSim.Rectangle(plane, o.width / 2, o.height / 2, uv, vv)
    
    OpticSim.Vis.draw!(scene, rect;  kwargs...)

    maybe_draw_debug_info(scene, o; transform=transform, kwargs...)
end


#-------------------------------------
# draw hexapolar origin
#-------------------------------------
function OpticSim.Vis.draw!(scene::MakieLayout.LScene, o::Origins.Hexapolar; transform::Geometry.Transform = Transform(), kwargs...) where {T<:Real}
    dir = forward(transform)
    uv = SVector{3}(right(transform))
    vv = SVector{3}(up(transform))
    pos = OpticSim.origin(transform)

    plane = OpticSim.Plane(dir, pos)
    ellipse = OpticSim.Ellipse(plane, o.halfsizeu, o.halfsizev, uv, vv)
    
    OpticSim.Vis.draw!(scene, ellipse;  kwargs...)

    maybe_draw_debug_info(scene, o; transform=transform, kwargs...)
end

#-------------------------------------
# draw source
#-------------------------------------
function OpticSim.Vis.draw!(scene::MakieLayout.LScene, s::Sources.Source{T}; parent_transform::Geometry.Transform = Transform(), debug::Bool=false, kwargs...) where {T<:Real}
   
    OpticSim.Vis.draw!(scene, s.origins;  transform=parent_transform * s.transform, debug=debug, kwargs...)

    if (debug)
        m = zeros(T, length(s), 7)
        for (index, optical_ray) in enumerate(s)
            ray = OpticSim.ray(optical_ray)
            ray = parent_transform * ray
            m[index, 1:7] = [ray.origin... ray.direction... OpticSim.power(optical_ray)]
        end
        
        m[:, 4:6] .*= m[:, 7] * ARRROW_LENGTH * visual_size(s.origins)  

        # Makie.arrows!(scene, [Makie.Point3f0(origin(ray))], [Makie.Point3f0(rayscale * direction(ray))]; kwargs..., arrowsize = min(0.05, rayscale * 0.05), arrowcolor = color, linecolor = color, linewidth = 2)
        color = :yellow
        Makie.arrows!(scene, m[:,1], m[:,2], m[:,3], m[:,4], m[:,5], m[:,6]; kwargs...,  arrowcolor = color, linecolor = color, linewidth = 2, arrowsize=ARRROW_SIZE * visual_size(s.origins))
    end

    # for ray in o
    #     OpticSim.Vis.draw!(scene, ray)
    # end
end

#-------------------------------------
# draw optical rays
#-------------------------------------
function OpticSim.Vis.draw!(scene::MakieLayout.LScene, rays::AbstractVector{OpticSim.OpticalRay{T, 3}}; kwargs...) where {T<:Real}
    m = zeros(T, length(rays)*2, 3)
    for (index, optical_ray) in enumerate(rays)
        ray = OpticSim.ray(optical_ray)
        m[(index-1)*2+1, 1:3] = [OpticSim.origin(optical_ray)...]
        m[(index-1)*2+2, 1:3] = [(OpticSim.origin(optical_ray) + OpticSim.direction(optical_ray) * 1 * OpticSim.power(optical_ray))... ]
    end
    
    color = :green
    Makie.linesegments!(scene, m[:,1], m[:,2], m[:,3]; kwargs...,  color = color, linewidth = 2, )
end

#-------------------------------------
# draw composite source
#-------------------------------------
function OpticSim.Vis.draw!(scene::MakieLayout.LScene, s::Sources.CompositeSource{T}; parent_transform::Geometry.Transform = Transform(), kwargs...) where {T<:Real}

    # @info "Composite"
    for source in s.sources
        # @info "Composite: local transform $(translation(s.transform))"
        # @info "Composite: child transform $(translation(source.transform))"
        OpticSim.Vis.draw!(scene, source; parent_transform=parent_transform*s.transform, kwargs...)
    end
    # OpticSim.Vis.draw!(scene, s.origins;  transform=s.transform, kwargs...)

    # m = zeros(T, length(s), 7)
    # for (index, optical_ray) in enumerate(s)
    #     ray = OpticSim.ray(optical_ray)
    #     m[index, 1:7] = [ray.origin... ray.direction... OpticSim.power(optical_ray)]
    # end
    
    # m[:, 4:6] .*= m[:, 7] * ARRROW_LENGTH * visual_size(s.origins)  

    # # Makie.arrows!(scene, [Makie.Point3f0(origin(ray))], [Makie.Point3f0(rayscale * direction(ray))]; kwargs..., arrowsize = min(0.05, rayscale * 0.05), arrowcolor = color, linecolor = color, linewidth = 2)
    # color = :yellow
    # Makie.arrows!(scene, m[:,1], m[:,2], m[:,3], m[:,4], m[:,5], m[:,6]; kwargs...,  arrowcolor = color, linecolor = color, linewidth = 2, arrowsize=ARRROW_SIZE * visual_size(s.origins))

    # # for ray in o
    # #     OpticSim.Vis.draw!(scene, ray)
    # # end
end
