# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

using Colors

# These utility functions provide a simple interface to the Emitters package to create emitters that are commonly used in optical systems. Many optical systems by convention have their optical axis parallel to the z axis. These emitters direct the rays in the negative z direction, toward the entrance of the optical system.

"""
    pointemitter(origin::AbstractVector{T}, coneangle; λ::Length = 500nm, numrays = 100) where {T<:Real}

Creates a point source with Lambertian emission power and cone distribution of rays, emitting in the -z direction. λ is
a unitful Length quantity, e.g., 550nm.
"""
function pointemitter(origin::AbstractVector{T}, coneangle; λ::Length = 500nm, numrays = 100) where {T<:Real}
    @assert length(origin) == 3

    return Sources.Source(
        Geometry.identitytransform(T), 
        Spectrum.DeltaFunction(λ),
        Origins.Point(Geometry.Vec3(origin)),
        Directions.UniformCone(-Geometry.unitZ3(T), coneangle, numrays),
        AngularPower.Lambertian()
    )
end

"""
    collimatedemitter(origin::AbstractVector{T}, halfsquaresize; λ::Length = 500nm, numrays = 100) where {T<:Real}

Creates a square collimated emitter, emitting rays in the -z direction. Rays are emitted on a square grid with
sqrt(numrays) on a side. λ can be a unitful quantity, e.g., 550nm, or a number. In the latter case the units are
implicitly microns.
"""
function collimatedemitter(origin::AbstractVector{T}, halfsquaresize; λ::Length = 500nm, numrays = 100) where {T<:Real}
    samples = Int(round(sqrt(numrays)))

    return Sources.Source(
        Geometry.translation(origin...),
        Spectrum.DeltaFunction(λ),
        Origins.RectGrid(2 * halfsquaresize, 2 * halfsquaresize, samples, samples),
        Directions.Constant(-Geometry.unitZ3(T)),
        AngularPower.Lambertian()
    )
end

"""
    imageemitter()


"""
function imageemitter(
    image::AbstractArray{<:AbstractGray},
    pixel_size::Tuple{T,T},
    pixel_pitch::Tuple{T,T},
    pixel_pos::Tuple{T,T};
    transform::Transform{T} = identitytransform(T)
) where {T<:Real}
    @assert length(size(image)) == 2

    pixel_source = Sources.Source(;
        transform = Geometry.translation(pixel_pos[1], pixel_pos[2], zero(T)),
        origins = Origins.RectGrid(pixel_size[1], pixel_size[2], 1, 1),
        directions = Directions.HexapolarCone(deg2rad(10), 2)
    )

    sources::Vector{Sources.CompositeSource{T}} = []
    for i = 1:size(image, 1)
        for j = 1:size(image, 2)
            if image[i, j] < 0.5
                continue
            end

            push!(
                sources,
                Sources.CompositeSource(
                    Geometry.translation((i - 1) * pixel_pitch[1], (j - 1) * pixel_pitch[2], zero(T)),
                    [pixel_source]
                )
            )
        end
    end

    return Sources.CompositeSource(transform, sources)
end
