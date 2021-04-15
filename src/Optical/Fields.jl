# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

"""
    HexapolarField(sys::AxisymmetricOpticalSystem; collimated = true, samples = 8, wavelength = 0.55, sourcepos = (0.0, 0.0, 3.0), sourceangle = 0.0, sourcenum = 0)

Distributes rays over the entrance pupil of the system in a hexapolar pattern.
"""
function HexapolarField(sys::AxisymmetricOpticalSystem{T}; kwargs...) where {T<:Real}
    return HexapolarField(semidiameter(sys) * 0.95, zero(SVector{3,T}); kwargs...)
end

"""
    HexapolarField(semidiameter, pupilpos; collimated = true, samples = 8, wavelength = 0.55, sourcepos = (0.0, 0.0, 3.0), sourceangle = 0.0, sourcenum = 0)

Distributes rays over a circular pupil with half-diameter defined by `semidiameter`, centred at `pupilpos` in a hexapolar pattern.
`samples` is the number of rings in the hexapolar pattern, so the number of rays in total is `samples * (samples + 1) / 2) * 6 + 1`.
"""
function HexapolarField(semidiameter::T, pupilpos::SVector{3,T}; collimated::Bool = true, samples::Int = 8, wavelength::T = T(0.55), sourcepos::SVector{3,T} = SVector{3,T}(0.0, 0.0, 3.0), sourceangle::T = zero(T), sourcenum::Int = 0) where {T<:Real}
    if sourceangle != zero(T)
        sourcepos = norm(sourcepos) .* SVector{3,T}(0.0, tan(sourceangle), 1.0) + pupilpos
    end
    pointsonpupil = Vector{SVector{3,T}}(undef, 0)
    for n in 0:samples
        if n == 0
            push!(pointsonpupil, pupilpos)
        else
            ρ = n / samples
            ninring = n * 6
            for i in 1:ninring
                ϕ = (i / ninring) * 2π
                u = cos(ϕ) * semidiameter
                v = sin(ϕ) * semidiameter
                push!(pointsonpupil, SVector{3,T}(ρ * u, ρ * v, 0) + pupilpos)
            end
        end
    end
    N = length(pointsonpupil)
    rays = Vector{OpticalRay{T,3}}(undef, N)
    if collimated
        direction = pupilpos - sourcepos
        for i in 1:N
            rays[i] = OpticalRay(pointsonpupil[i] - direction, direction, one(T), wavelength, sourcenum = sourcenum)
        end
        raygenerator = RayListSource(rays)
    else
        for i in 1:N
            rays[i] = OpticalRay(sourcepos, pointsonpupil[i] - sourcepos, one(T), wavelength, sourcenum = sourcenum)
        end
        raygenerator = RayListSource(rays)
    end
    return raygenerator
end
export HexapolarField


"""
    GridField(sys::AxisymmetricOpticalSystem; collimated = true, samples = 20, wavelength = 0.55, sourcepos = (0.0, 0.0, 3.0), sourceangle = 0.0, sourcenum = 0)

Distributes rays over the entrance pupil of the system in a rectangular grid pattern.
"""
function GridField(sys::AxisymmetricOpticalSystem{T}; kwargs...) where {T<:Real}
    return GridField(semidiameter(sys) * 0.95, zero(SVector{3,T}); kwargs...)
end

"""
    GridField(semidiameter, pupilpos; collimated = true, samples = 20, wavelength = 0.55, sourcepos = (0.0, 0.0, 3.0), sourceangle = 0.0, sourcenum = 0)

Distributes rays over a circular pupil with half-diameter defined by `semidiameter`, centred at `pupilpos` in a rectangular grid pattern.
`samples` is the number of rays on each side of the grid, so there are `samples×samples` rays in total.
"""
function GridField(semidiameter::T, pupilpos::SVector{3,T}; collimated::Bool = true, samples::Int = 20, wavelength::T = T(0.55), sourcepos::SVector{3,T} = SVector{3,T}(0.0, 0.0, 3.0), sourceangle::T = zero(T), sourcenum::Int = 0) where {T<:Real}
    if sourceangle != zero(T)
        sourcepos = norm(sourcepos) .* SVector{3,T}(0.0, tan(sourceangle), 1.0) + pupilpos
    end
    pointsonpupil = Vector{SVector{3,T}}(undef, 0)
    for x in 0:samples
        for y in 0:samples
            u = (2.0 * (x / samples) - one(T)) * semidiameter
            v = (2.0 * (y / samples) - one(T)) * semidiameter
            if (u == zero(T) && v == zero(T)) || sqrt(u^2 + v^2) < semidiameter
                push!(pointsonpupil, SVector{3,T}(u, v, 0) + pupilpos)
            end
        end
    end

    N = length(pointsonpupil)
    rays = Vector{OpticalRay{T,3}}(undef, N)
    if collimated
        direction = pupilpos - sourcepos
        for i in 1:N
            rays[i] = OpticalRay(pointsonpupil[i] - direction, direction, one(T), wavelength, sourcenum = sourcenum)
        end
        raygenerator = RayListSource(rays)
    else
        for i in 1:N
            rays[i] = OpticalRay(sourcepos, pointsonpupil[i] - sourcepos, one(T), wavelength, sourcenum = sourcenum)
        end
        raygenerator = RayListSource(rays)
    end
    return raygenerator
end
export GridField
