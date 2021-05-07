# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

# These utility functions provide a simple interface to the Emitters package to create emitters that are commonly used in optical systems. Many optical systems by convention have their optical axis parallel to the z axis. These emitters direct the rays in the negative z direction, toward the entrance of the optical system.

"""
Creates a point source with Lambertian emission power and cone distribution of rays, emitting in the -z direction. λ is a unitful Length quantity, e.g., 550nm.
"""
function pointemitter(origin::AbstractVector{T},coneangle;λ::Length = 500nm,numrays = 100) where{T<:Real}
    @assert length(origin) == 3

    return Sources.Source(
        Geometry.identitytransform(T), 
        Spectrum.DeltaFunction(λ),
        Origins.Point(Geometry.Vec3(origin)),
        Directions.UniformCone(-Geometry.unitZ3(T),coneangle,numrays),
        AngularPower.Lambertian()
        )
end

"""
Creates a square collimated emitter, emitting rays in the -z direction. Rays are emitted on a square grid with sqrt(numrays) on a side. λ can be a unitful quantity, e.g., 550nm, or a number. In the latter case the units are implicitly microns.
"""
function collimatedemitter(origin::AbstractVector{T},halfsquaresize;λ::Length = 500nm,numrays = 100) where{T<:Real}
    samples = Int(round(sqrt(numrays)))

    return Sources.Source(
        Geometry.translation(origin...),
        Spectrum.DeltaFunction(λ),
        Origins.RectGrid(2*halfsquaresize,2*halfsquaresize,samples,samples),
        Directions.Constant(-Geometry.unitZ3(T)),
        AngularPower.Lambertian())
end
