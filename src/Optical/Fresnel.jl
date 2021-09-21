# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

mᵢandmₜ(matout, matin, surfacenormal::SVector{N,T}, r::SVector{N,T}) where {T<:Real,N} = dot(surfacenormal, r) < zero(T) ? (matout, matin) : (matin, matout)

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

    if isa(nₜ,Complex)
        sinθₜ = zero(T) #this is a bit hacky. When the transmitted index is complex, which it is for metals, sinθₜ is not defined. Set it to zero. It will be ignored in the specialized version of fresnel which handles complex index of refraction.
    else
        sinθₜ = nᵢ / nₜ * sinθᵢ
    end
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

    return normalize(r - 2 * nₛ * (tdot))
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
For geometric ray tracing this coefficient can be used directly to compute intensity on the detector plane. For Huygens phase optics need to take the square root to compute the amplitude.
The power of the transmitted and refracted rays may not sum to one because of the area correction applied to the transmitted component.
The intensity per area can increase or decrease depending on the indices of refraction.

`nᵢ` is the RI of the material which the incident ray travels in, `nₜ` is the RI of the material the transmitted ray travels in.
`sinθᵢ` and `sinθₜ` are the sin of the angles of incidence and transmission respectively.
"""
function fresnel(nᵢ::T, nₜ::T, sinθᵢ::T, sinθₜ::T) where {T<:Real}
    if (sinθᵢ >= nₜ / nᵢ) # 100% reflectance, zero transmission
        return (one(T), zero(T))
    end

    cosθₜ = sqrt(one(T) - sinθₜ^2)
    cosθᵢ = sqrt(one(T) - sinθᵢ^2)

    # nᵢ*sin(θᵢ) = nₜ*sin(θₜ)
    rₛ = ((nᵢ * cosθᵢ - nₜ * cosθₜ) / (nᵢ * cosθᵢ + nₜ * cosθₜ))
    tₛ = (2nᵢ * cosθᵢ / (nᵢ * cosθᵢ + nₜ * cosθₜ))
    rₚ = ((nₜ * cosθᵢ - nᵢ * cosθₜ) / (nₜ * cosθᵢ + nᵢ * cosθₜ))
    tₚ = (2nᵢ * cosθᵢ / (nₜ * cosθᵢ + nᵢ * cosθₜ))

    # transmitted amplitude scale factor
    Tₐ = ((nₜ * cosθₜ) / (nᵢ * cosθᵢ))

    return rₛ,tₛ,rₚ,tₚ,Tₐ
end

"""Fresnel equations for dielectric/metal interface. Metal will have complex index of refraction. This takes much longer than the case for real index of refraction so use specialized function. The assumption is that the light not reflected from the metal surface is completely absorbed. This will need to be extended for partially silvered mirrors."""
function fresnelcomplex(nᵢ::T, nₜ::Complex{T}, sinθᵢ::T) where {T<:Real}
    n² = (nₜ/nᵢ)^2

    cosθᵢ = sqrt(one(T) - sinθᵢ^2)

    rₛ = (cosθᵢ - sqrt(n²-sinθᵢ^2))/(cosθᵢ + sqrt(n²-sinθᵢ^2))
    rₚ = (n²*cosθᵢ - sqrt(n²-sinθᵢ^2))/(n²*cosθᵢ + sqrt(n²-sinθᵢ^2))

    return rₛ,T(0),rₚ,T(0),T(1) #for metals transmission is assumed to be zero. Tₐ = 1 because reflected area = incident area
end


############################################################################################################################

"""assuming random polarization this computes the average power across all polarization angles for the reflected and transmitted paths"""
function averagepower(rₛ, rₚ, tₛ, tₚ, Tₐ)
    # assuming random polarization
    powᵣ = 0.5 * (rₛ + rₚ)
    powₜ = 0.5 * (tₛ + tₚ) * Tₐ  #compensate for change in intensity per area caused by the differing angles of incidence and refraction as well as differing indices.
    return powᵣ,powₜ
end

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
function processintersection(opticalinterface::FresnelInterface{T}, point::SVector{N,T}, normal::SVector{N,T}, incidentray::OpticalRay{T,N,OpticSim.Polarization.NoPolarization{T}}, temperature::T, pressure::T, test::Bool, firstray::Bool = false) where {T<:Real,N}
    λ = wavelength(incidentray)
    mᵢ, mₜ = mᵢandmₜ(outsidematerialid(opticalinterface), insidematerialid(opticalinterface), normal, direction(incidentray))
    nᵢ = one(T)
    nₜ = one(T)
    α = zero(T)
    if !isair(mᵢ)
        mat = glassforid(mᵢ)::OpticSim.GlassCat.Glass
        nᵢ = index(mat, λ, temperature = temperature, pressure = pressure)::T
        α = absorption(mat, λ, temperature = temperature, pressure = pressure)::T
    end
    if !isair(mₜ)
        nₜ = index(glassforid(mₜ)::OpticSim.GlassCat.Glass, λ, temperature = temperature, pressure = pressure)::T
    end
    sinθᵢ, sinθₜ = snell(normal, direction(incidentray), nᵢ, nₜ)
    rₛ,tₛ,rₚ,tₚ,Tₐ = fresnel(nᵢ, nₜ, sinθᵢ, sinθₜ)

    
    powᵣ,powₜ = averagepower(rₛ,tₛ,rₚ,tₚ,Tₐ)

    incident_pow = power(incidentray)

    # optical distance from ray origin to point of intersection in mm. Compensate for the fact that the ray has been slightly shortened.
    geometricpathlength = norm(point - origin(incidentray)) + (firstray ? zero(T) : T(RAY_OFFSET))
    thisraypathlength = nᵢ * geometricpathlength
    raypathlength = pathlength(incidentray) + thisraypathlength

    # compute updated power based on absorption coefficient of material using Beer's law
    internal_trans = one(T)
    if α > zero(T)
        internal_trans = exp(-α * geometricpathlength)
    end

    # TODO - this is an approximation (total hack) for now until we get better modeling of thin film reflectors, antireflection coatings, etc.
    powᵣ = max(powᵣ, reflectance(opticalinterface)) * internal_trans
    powₜ = powₜ * transmission(opticalinterface) * internal_trans

    # generate new rays using Monte Carlo sampling proportional to power. For most optical surfaces the vast majority of rays will be refracted rays.
    # So could leave this turned on all the time with little impact on performance and get approximate scattering effects for free.
    r = !test * rand()
    # assuming (powᵣ + powₜ) <= 1 (asserted in constructor)
    if interfacemode(opticalinterface) == Transmit || (interfacemode(opticalinterface) == ReflectOrTransmit && r < powₜ)
        # refraction
        raydirection = refractedray(nᵢ, nₜ, normal, direction(incidentray))
        raypower = powₜ * incident_pow
    elseif interfacemode(opticalinterface) == Reflect || (interfacemode(opticalinterface) == ReflectOrTransmit && r < (powᵣ + powₜ))
        # reflection
        raypower = powᵣ * incident_pow
        raydirection = reflectedray(normal, direction(incidentray))
    else
        return nothing
    end

    if raydirection === nothing
        return nothing
    else
        return raydirection, raypower, raypathlength, Polarization.NoPolarization{T}()
    end
end

"""extracts the indices of refraction from the optical interface in the correct order, n₁, nₜ, where nᵢ is the incident index of refraction and nₜ is the transmitted index of refraction"""
function refractiveindices(opticalinterface::FresnelInterface{T},normal::SVector{N,T},λ::T,incidentray::SVector{N,T},temperature,pressure) where{T,N}
    mᵢ, mₜ = mᵢandmₜ(outsidematerialid(opticalinterface), insidematerialid(opticalinterface), normal, incidentray)
    nᵢ = one(T)
    nₜ = one(T)
    α = zero(T)
    if !isair(mᵢ)
        mat = glassforid(mᵢ)::OpticSim.GlassCat.Glass
        nᵢ = index(mat, λ, temperature = temperature, pressure = pressure)::T
        α = absorption(mat, λ, temperature = temperature, pressure = pressure)::T
    end
    if !isair(mₜ)
        #need to add extensions to GlassCat to allow for complex indices of refraction for metals.
        nₜ = index(glassforid(mₜ)::OpticSim.GlassCat.Glass, λ, temperature = temperature, pressure = pressure)::T
    end

    return nᵢ, nₜ
end


"""
processintersection(opticalinterface::OpticalInterface{T}, point::SVector{N,T}, normal::SVector{N,T}, incidentray::OpticalRay{T,N}, temperature::T, pressure::T, ::Bool, firstray::Bool = false) -> Tuple{SVector{N,T}, T, T}

Specialized method for polarized Fresnel interfaces. This takes much more time to compute so don't want to make it the default.

Processes an intersection of an [`OpticalRay`](@ref) with an [`OpticalInterface`](@ref), distinct behaviors must be implemented for each subclass of `OpticalInterface`.

`point` is the 3D intersection point in global space, `normal` is the surface normal at the intersection point.

If `test` is true then the behavior of the ray should be deterministic.
`firstray` indicates that this ray is the first segment of the trace and therefore the origin is not offset.

The values returned are the normalized direction of the ray after the intersection, the _instantaneous_ power of the ray after the intersection and the optical path length of the ray up to the intersection.

`nothing` is returned if the ray should stop here, in order to obtain the correct intensity on the detector through monte carlo integration `nothing` should be returned proportionally to create the correct power distribution.
i.e. If the interface should modulate power to 76% then 24% of calls to this function should return `nothing`.
"""
function processintersection(opticalinterface::FresnelInterface{T}, point::SVector{N,T}, normal::SVector{N,T}, incidentray::OpticalRay{T,N,Polarization.Chipman{T}}, temperature::T, pressure::T, test::Bool, firstray::Bool = false) where {T<:Real,N}
    nᵢ,nₜ = refractiveindices(opticalinterface, normal, wavelength(incidentray), direction(incidentray),temperature,pressure)
    
    sinθᵢ, sinθₜ = snell(normal, direction(incidentray), nᵢ, nₜ)

    #TODO the choice of which fresnel function to use depending on whether nₜ is complex should be hidden better in fresnel function. Maybe put the call to snell in fresnel?
    if is(nₜ,Complex)
        rₛ,tₛ,rₚ,tₚ,Tₐ = fresnelcomplex(nᵢ, nₜ, sinθᵢ)
    else
        rₛ,tₛ,rₚ,tₚ,Tₐ = fresnel(nᵢ, nₜ, sinθᵢ, sinθₜ)
    end
   
    powᵣ,powₜ = averagepower(rₛ,tₛ,rₚ,tₚ,Tₐ)
    incident_pow = power(incidentray)

    # optical distance from ray origin to point of intersection in mm. Compensate for the fact that the ray has been slightly shortened.
    geometricpathlength = norm(point - origin(incidentray)) + (firstray ? zero(T) : T(RAY_OFFSET))
    thisraypathlength = nᵢ * geometricpathlength
    raypathlength = pathlength(incidentray) + thisraypathlength

    # compute updated power based on absorption coefficient of material using Beer's law
    internal_trans = one(T)

    #TODO fix this after getting polarization stuff working
    # if α > zero(T)
    #     internal_trans = exp(-α * geometricpathlength)
    # end

    r = !test * rand()


    # assuming (powᵣ + powₜ) <= 1 (asserted in constructor)
    if interfacemode(opticalinterface) == Transmit || (interfacemode(opticalinterface) == ReflectOrTransmit && r < powₜ)
        # refraction
        raydirection = refractedray(nᵢ, nₜ, normal, direction(incidentray))
        raypower = powₜ * incident_pow
        polarizationinfo = composepolarization(Complex(tₛ),Complex(tₚ),Tₐ,normal,direction(incidentray),raydirection,polarization(incidentray))
    elseif interfacemode(opticalinterface) == Reflect || (interfacemode(opticalinterface) == ReflectOrTransmit && r < (powᵣ + powₜ))
        # reflection
        raypower = powᵣ * incident_pow
        raydirection = reflectedray(normal, direction(incidentray))
        polarizationinfo = composepolarization(Complex(rₛ),Complex(rₚ),Tₐ,normal,direction(incidentray),raydirection,polarization(incidentray))
    else
        return nothing
    end

    if raydirection === nothing
        return nothing
    else

        return normalize(raydirection), raypower, raypathlength, polarizationinfo
    end

end
