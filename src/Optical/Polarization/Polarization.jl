module Polarization

using StaticArrays
using LinearAlgebra
using ..OpticSim
import Unitful
import Unitful.DefaultSymbols

abstract type AbstractPolarization{T<:Real} end

struct Chipman{T<:Real} <: AbstractPolarization{T}
    P::SMatrix{3,3,Complex{T},9} #3D polarization matrix
    
    Chipman{T}() where{T<:Real}= new{T}(SMatrix{3,3,Complex{T},9}(0.0 + 0.0im,0.0 + 0.0im,0.0 + 0.0im,0.0 + 0.0im,0.0 + 0.0im,0.0 + 0.0im,0.0 + 0.0im,0.0 + 0.0im,0.0 + 0.0im))
end

"""Computes the P matrix (Polarized Light and Optical Systems by Chimpan et al, chapter 9) which represents the s-p coordinate frame for an interface. In OpticSim surface normals are defined positive pointing outward from a solid. In Chipman the surface normal is defined positive on the exiting side of an interface. Need to take this account when passing `surfacenormal` argument."""
function localtoworld(surfacenormal::AbstractVector{T},incidentvector::AbstractVector{T}) where{T<:Real}
    k̂ = normalize(incidentvector)
    ŝ = normalize(k̂ × surfacenormal)
    p̂ = incidentvector × surfacenormal
    return SMatrix{3,3,T,9}(ŝ[1],ŝ[2],ŝ[3],p̂[1],p̂[2],p̂[3],k̂[1],k̂[2],k̂[3])
end

worldtolocal(loctoworld::SMatrix{3,3,T,9}) where{T<:Real}= loctoworld'

"""For Fresnel reflection need to compute reflected and/or refracted P matrix. This looks like:
see if normal is on the """

# function composepmatrix(pinputmatrix::Chipman{T},interface::OpticalInterface{T}, ray::OpticalRay{T}) where{T<:Real}
#     #compute worldtolocal for incident ray
#     plocal = pinputmatrix*worldtolocal(localtoworld(normal(interface),direction(ray)))
#     #compute s,p components for reflected and transmitted values these are Jones matrix values, potentially complex.
#     #compute reflected,refracted transformation

# end


# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

mᵢandmₜ(matout, matin, surfacenormal::SVector{N,T}, r::AbstractRay{T,N}) where {T<:Real,N} = dot(surfacenormal, direction(r)) < zero(T) ? (matout, matin) : (matin, matout)

"""
    snell(surfacenormal::AbstractVector{T}, raydirection::AbstractVector{T}, nᵢ::T, nₜ::T) -> Tuple{T,T}

`nᵢ` is the index of refraction on the incidence side of the interface.
`nₜ` is the index of refraction on the transmission side.

Returns `sinθᵢ` and `sinθₜ` according to [Snell's law](https://en.wikipedia.org/wiki/Snell%27s_law).
"""
function snell(surfacenormal::S, raydirection::S, nᵢ::T, nₜ::T) where {T<:Real,S<:AbstractArray{T}}
    # snells law: nᵢsin(θᵢ) = nₜsin(θₜ)
    # sin(θ) = sin(180-θ) so don't have to worry about reversing sign of r

    # prevent div by 0
    if surfacenormal == raydirection || surfacenormal == -raydirection
        return zero(T), zero(T)
    end

    sinθᵢ = norm(cross(surfacenormal, raydirection))
    #only interested in magnitude not direction of the vector normal to surfacenormal and raydirection.
    # If wanted vector then should use -raydirection.

    @assert zero(T) <= sinθᵢ

    sinθₜ = nᵢ / nₜ * sinθᵢ
    @assert zero(T) <= sinθₜ

    return sinθᵢ, sinθₜ
end

function reflectedray(surfacenormal::S, raydirection::S) where {T<:Real,S<:AbstractArray{T}}
    r = raydirection
    nₛ = surfacenormal
    tdot = dot(r, nₛ)

    if tdot == zero(T)
        return -r
    end

    return r - 2 * nₛ * (tdot)
end

function refractedray(incidenceindex::T, transmittedindex::T, surfacenormal::S, raydirection::S) where {T<:Real,S<:AbstractArray{T,1}}

    sinθᵢ, sinθₜ = snell(surfacenormal, raydirection, incidenceindex, transmittedindex)
    # (nᵢ, nₜ) = nᵢandnₜ(incidenceindex, transmittedindex, surfacenormal, raydirection) #redundant computation here checking for nₛ⋅r < 0 in two places.
    # Optimize later if compiler doesn't optimize this away by inlining.

    if (sinθᵢ >= transmittedindex / incidenceindex) # 100% reflectance, zero transmission
        return nothing
    else
        nₛ = surfacenormal
        r = raydirection

        sgn = sign(dot(nₛ, r))
        pₜ = r - nₛ * dot(r, nₛ)
        if sinθᵢ == zero(T)
            refractedparallel = SVector{3,T}(0.0, 0.0, 0.0)
        else
            refractedparallel = pₜ * sinθₜ / sinθᵢ
        end
        refractedperpendicular = sgn * nₛ * sqrt(one(T) - sinθₜ^2)
        refracted = refractedparallel + refractedperpendicular

        return normalize(refracted)
    end
end

"""
    fresnel(nᵢ::T, nₜ::T, sinθᵢ::T, sinθₜ::T) -> Tuple{T,T}

Returns reflectance and tranmission power coefficients according to the [Fresnel equations](https://en.wikipedia.org/wiki/Fresnel_equations).

The power of the transmitted and refracted rays may not sum to one because of the area correction applied to the transmitted component.
The intensity per area can increase or decrease depending on the indices of refraction. This is recorded in the Tₐ variable.

`nᵢ` is the RI of the material which the incident ray travels in, `nₜ` is the RI of the material the transmitted ray travels in.
`sinθᵢ` and `sinθₜ` are the sin of the angles of incidence and transmission respectively.
"""
function fresnel(nᵢ::T, nₜ::T, sinθᵢ::T, sinθₜ::T) where {T<:Real}
    if (sinθᵢ >= nₜ / nᵢ) # 100% reflectance, zero transmission
        return (one(T), zero(T))
    end

    cosθₜ = sqrt(one(T) - sinθₜ^2)
    cosθᵢ = sqrt(one(T) - sinθᵢ^2)

    rₛ = (nᵢ * cosθₜ - nₜ * cosθᵢ) / (nᵢ * cosθₜ + nₜ * cosθᵢ)
    tₛ = 2nᵢ * cosθᵢ / (nᵢ * cosθₜ + nₜ * cosθᵢ)
    rₚ = (nᵢ * cosθᵢ - nₜ * cosθₜ) / (nᵢ * cosθᵢ + nₜ * cosθₜ)
    tₚ = 2nᵢ * cosθᵢ / (nᵢ * cosθᵢ + nₜ * cosθₜ)

    # transmitted amplitude scale factor
    Tₐ = ((nₜ * cosθₜ) / (nᵢ * cosθᵢ))

    return rₛ,tₛ,rₚ,tₚ,Tₐ
end

"""Fresnel equations for dielectric/metal interface. Metal will have complex index of refraction"""
function fresnel(nᵢ::T, nₜ::Complex{T}, sinθᵢ::T) where {T<:Real}
    n² = (nₜ/nᵢ)^2

    cosθᵢ = sqrt(one(T) - sinθᵢ^2)

    rₛ = (cosθᵢ - sqrt(n²-sinθᵢ^2))/(cosθᵢ + sqrt(n²-sinθᵢ^2))
    rₚ = (n²*cosθᵢ - sqrt(n²-sinθᵢ^2))/(n²*cosθᵢ + sqrt(n²-sinθᵢ^2))

    return rₛ,T(0),rₚ,T(0),T(1) #for metals transmission is assumed to be zero. Tₐ = 1 because reflected area = incident area
end

function aluminumfresnel()
    #index of refraction of aluminum at 633nm is 1.374 + 7.620im
    nₜ = 1.374 + 7.620im
    nᵢ = 1.0 #assume air
    s = Vector{Complex{Float64}}(undef,0)
    p = Vector{Complex{Float64}}(undef,0)

    for θ in 0.0:.01:π/2
        rₛ,_,rₚ,_,_ = fresnel(nᵢ,nₜ,sin(θ))
        push!(s,rₛ)
        push!(p,rₚ)
    end

    return s,p
end

############################################################################################################################

"""
    processintersection(opticalinterface::OpticalInterface{T}, point::SVector{N,T}, normal::SVector{N,T}, incidentray::OpticalRay{T,N}, temperature::T, pressure::T, ::Bool, firstray::Bool = false) -> Tuple{SVector{N,T}, T, T}

Processes an intersection of an [`OpticalRay`](@ref) with an [`OpticalInterface`](@ref), distinct behaviors must be implemented for each subclass of `OpticalInterface`.

`point` is the 3D intersection point in global space, `normal` is the surface normal at the intersection point.

If `test` is true then the behavior of the ray should be deterministic.
`firstray` indicates that this ray is the first segment of the trace and therefore the origin is not offset.

The values returned are the normalized direction of the ray after the intersection, the _instantaneous_ power of the ray after the intersection and the optical path length of the ray up to the intersection.

`nothing` is returned if the ray should stop here, in order to obtain the correct intensity on the detector through monte carlo integration `nothing` should be returned proportionally to create the correct power distribution.
i.e. If the interface should modulate power to 76% then 24% of calls to this function should return `nothing`.
"""
# function processintersection(opticalinterface::FresnelInterface{T}, point::SVector{N,T}, normal::SVector{N,T}, incidentray::OpticalRay{T,N}, temperature::T, pressure::T, test::Bool, firstray::Bool = false) where {T<:Real,N}
#     λ = wavelength(incidentray)
#     mᵢ, mₜ = mᵢandmₜ(outsidematerialid(opticalinterface), insidematerialid(opticalinterface), normal, incidentray)
#     nᵢ = one(T)
#     nₜ = one(T)
#     α = zero(T)
#     if !isair(mᵢ)
#         mat = glassforid(mᵢ)::OpticSim.GlassCat.Glass
#         nᵢ = index(mat, λ, temperature = temperature, pressure = pressure)::T
#         α = absorption(mat, λ, temperature = temperature, pressure = pressure)::T
#     end
#     if !isair(mₜ)
#         nₜ = index(glassforid(mₜ)::OpticSim.GlassCat.Glass, λ, temperature = temperature, pressure = pressure)::T
#     end
#     (sinθᵢ, sinθₜ) = snell(normal, direction(incidentray), nᵢ, nₜ)
#     (powᵣ, powₜ) = fresnel(nᵢ, nₜ, sinθᵢ, sinθₜ)

#     incident_pow = power(incidentray)

#     # optical distance from ray origin to point of intersection in mm. Compensate for the fact that the ray has been slightly shortened.
#     geometricpathlength = norm(point - origin(incidentray)) + (firstray ? zero(T) : T(RAY_OFFSET))
#     thisraypathlength = nᵢ * geometricpathlength
#     raypathlength = pathlength(incidentray) + thisraypathlength

#     # compute updated power based on absorption coefficient of material using Beer's law
#     internal_trans = one(T)
#     if α > zero(T)
#         internal_trans = exp(-α * geometricpathlength)
#     end

#     # TODO - this is an approximation (total hack) for now until we get better modeling of thin film reflectors, antireflection coatings, etc.
#     powᵣ = max(powᵣ, reflectance(opticalinterface)) * internal_trans
#     powₜ = powₜ * transmission(opticalinterface) * internal_trans

#     # generate new rays using Monte Carlo sampling proportional to power. For most optical surfaces the vast majority of rays will be refracted rays.
#     # So could leave this turned on all the time with little impact on performance and get approximate scattering effects for free.
#     r = !test * rand()
#     # assuming (powᵣ + powₜ) <= 1 (asserted in constructor)
#     if interfacemode(opticalinterface) == Transmit || (interfacemode(opticalinterface) == ReflectOrTransmit && r < powₜ)
#         # refraction
#         raydirection = refractedray(nᵢ, nₜ, normal, direction(incidentray))
#         raypower = powₜ * incident_pow
#     elseif interfacemode(opticalinterface) == Reflect || (interfacemode(opticalinterface) == ReflectOrTransmit && r < (powᵣ + powₜ))
#         # reflection
#         raypower = powᵣ * incident_pow
#         raydirection = reflectedray(normal, direction(incidentray))
#     else
#         return nothing
#     end

#     if raydirection === nothing
#         return nothing
#     else
#         return normalize(raydirection), raypower, raypathlength
#     end
# end


end #module
export Polarization