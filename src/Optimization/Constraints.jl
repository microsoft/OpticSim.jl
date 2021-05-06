# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

"""Interfaces, and some implementations, for common constraints widely used in optical optimization"""

#Notes: we probably need a OpticalElement type, which represents a physical lens object with a defined shape. One can compute thickness of an element, i.e., the distance through an object with a homogeneous material when pierced by a ray with a specific origin and direction. LensAssembly is too general to serve this purpose.
# Need a generalization of AbstractOpticalSystem with more than one detector.

# 

#maybe something like this:

abstract type OpticalElement end
function thickness(a::OpticalElement) end

struct Lens<:OpticalElement
    samplingdirection::Ray
    samplingPattern # note sure yet
end

function intersectionpoints() end
function centroid(a::AbstractVector) end

#standard measurements

function RMS_spot_size(a::AbstractVector{T}, b::AxisymmetricOpticalSystem{T}, samples::Int = 3) where {T}
    # RMSE spot size
    lens = Optimization.updateoptimizationvariables(b, a)
    field = Sources.Source(transform = Geometry.translation(0.0,0.0,10.0), origins = Origins.Hexapolar(5,5.0,5.0),directions = Directions.Constant(0.0,0.0,-1.0))
    error = zero(T)
    hits = 0
    for r in field
        traceres = OpticSim.trace(lens, r, test = true)
        if traceres !== nothing
            hitpoint = point(traceres)
            if abs(hitpoint[1]) > eps(T) && abs(hitpoint[2]) > eps(T)
                dist_to_axis = hitpoint[1]^2 + hitpoint[2]^2
                error += dist_to_axis
            end
            hits += 1
        end
    end
    if hits > 0
        error = sqrt(error / hits)
    end
    return error
end

function RMS_spot_size(system::OpticalSystem, rays::Emitters.Sources.Source, samples::Int = 3)
    error = zero(T)
    hits = 0
    #compute centroid of rays
    
    for r in rays
        traceres = OpticSim.trace(system, r, test = true)
        if traceres !== nothing
            temp = centroid-point(traceres)
            error += temp⋅temp
            hits += 1
        end
    end
    if hits > 0
        return error / hits
    else
        return nothing
    end
end

# function RMS_spot_size(lens, x::T...) where {T}
#     lens = Optimization.updateoptimizationvariables(lens, collect(x))
#     field = HexapolarField(lens, collimated = true, samples = 10)
#     error = zero(T)
#     hits = 0
#     for r in field
#         traceres = OpticSim.trace(lens, r, test = true)

#         if !(nothing === traceres)
#             hitpoint = point(traceres)
#             if abs(hitpoint[1]) > eps(T) && abs(hitpoint[2]) > eps(T)
#                 dist_to_axis = sqrt(hitpoint[1]^2 + hitpoint[2]^2)
#                 error += dist_to_axis
#             end
#             hits += 1
#         end
#     end

#     error /= hits

#     return error
# end

