# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

"""
    WrapperSurface{T,S<:Surface{T}} <: Surface{T}

A generic surface type which serves as a basis for extension of [`Surface`](@ref)s for custom [`OpticalInterface`](@ref) subclasses.
Essentially just forwards all `Surface` and `ParametricSurface` methods to a field of the `WrapperSurface` named `surface`.
Also provides a generic implementation of [`surfaceintersection`](@ref) which tests for an intersection with the underlying surface and returns either an [`EmptyInterval`](@ref) or a half space (never a closed interval).
"""
abstract type WrapperSurface{T,S<:Surface{T}} <: Surface{T} end
export WrapperSurface

# forward all methods, can be specialised if not applicable for subclass
interface(r::WrapperSurface{T}) where {T<:Real} = r.interface
centroid(r::WrapperSurface{T,S}) where {T<:Real,S<:Surface{T}} = centroid(r.surface::S)
normal(r::WrapperSurface{T,S}) where {T<:Real,S<:Surface{T}} = normal(r.surface::S)
normal(r::WrapperSurface{T}, ::T, ::T) where {T<:Real} = normal(r)
point(r::WrapperSurface{T,S}, u::T, v::T) where {T<:Real,S<:Surface{T}} = point(r.surface::S, u, v)
uv(r::WrapperSurface{T,S}, x::T, y::T, z::T) where {T<:Real,S<:Surface{T}} = uv(r, SVector{3,T}(x, y, z))
uv(r::WrapperSurface{T,S}, p::SVector{3,T}) where {T<:Real,S<:Surface{T}} = uv(r.surface::S, p)
uvrange(::Type{WrapperSurface{T,S}}) where {T<:Real,S<:Surface{T}} = uvrange(S)
uvrange(::WrapperSurface{T,S}) where {T<:Real,S<:Surface{T}} = uvrange(S)
onsurface(a::WrapperSurface{T,S}, x::T, y::T, z::T) where {T<:Real,S<:Surface{T}} = onsurface(a.surface::S, x, y, z)
onsurface(a::WrapperSurface{T,S}, p::SVector{3,T}) where {T<:Real,S<:Surface{T}} = onsurface(a.surface::S, p)
inside(a::WrapperSurface{T,S}, x::T, y::T, z::T) where {T<:Real,S<:Surface{T}} = inside(a.surface::S, x, y, z)
inside(a::WrapperSurface{T,S}, p::SVector{3,T}) where {T<:Real,S<:Surface{T}} = inside(a.surface::S, p)
partials(a::WrapperSurface{T,S}, u::T, v::T) where {T<:Real,S<:Surface{T}} = partials(a.surface::S, u, v)
makemesh(l::WrapperSurface{T,S}, subdivisions::Int = 20) where {T<:Real,S<:Surface{T}} = makemesh(l.surface::S, subdivisions)

function surfaceintersection(l::WrapperSurface{T,S}, r::AbstractRay{T,3}) where {T<:Real,S<:Surface{T}}
    itvl = surfaceintersection(l.surface::S, r)
    if itvl isa EmptyInterval{T}
        return EmptyInterval(T)
    elseif itvl isa DisjointUnion{T}
        throw(ErrorException("Grating/HOE must use a simple (planar) surface for now"))
    else
        intsct = halfspaceintersection(itvl)
        u, v = uv(intsct)
        intsct = Intersection(α(intsct), point(intsct), normal(intsct), u, v, interface(l), flippednormal = flippednormal(intsct))
        if dot(normal(intsct), direction(r)) > zero(T)
            return rayorigininterval(intsct)
        else
            return positivehalfspace(intsct)
        end
    end
end

######################################################################

"""
    ThinGratingSurface{T,S} <: WrapperSurface{T,S}

Surface type for use with [`ThinGratingInterface`](@ref).

```julia
ThinGratingSurface(surface::Surface{T}, interface::ThinGratingInterface{T})
```
"""
struct ThinGratingSurface{T,S} <: WrapperSurface{T,S}
    surface::S
    interface::ThinGratingInterface{T}
    ThinGratingSurface(surface::S, interface::ThinGratingInterface{T}) where {T<:Real,S<:Surface{T}} = new{T,S}(surface, interface)
end
export ThinGratingSurface

interface(r::ThinGratingSurface{T}) where {T<:Real} = r.interface

function processintersection(opticalinterface::ThinGratingInterface{T}, point::SVector{N,T}, normal::SVector{N,T}, incidentray::OpticalRay{T,N}, temperature::T, pressure::T, test::Bool, firstray::Bool = false) where {T<:Real,N}
    # we want the surface to work if hit from either side, so we just reverse the normal if it is hit on the back side
    if dot(direction(incidentray), normal) > zero(T)
        normal = -normal
    end

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

    incident_pow = power(incidentray)
    # compute updated power based on absorption coefficient of material using Beer's law
    internal_trans = one(T)
    if α > zero(T)
        internal_trans = exp(-α * geometricpathlength)
    end

    # optical distance from ray origin to point of intersection in mm. Compensate for the fact that the ray has been slightly shortened.
    geometricpathlength = norm(point - origin(incidentray)) + (firstray ? zero(T) : T(RAY_OFFSET))
    thisraypathlength = nᵢ * geometricpathlength
    raypathlength = pathlength(incidentray) + thisraypathlength

    λ = wavelength(incidentray)
    r = !test * rand()
    if opticalinterface.period < λ / nₜ
        # TODO!! what to do in this case?
        raypower = transmission(opticalinterface, 0) * internal_trans
        if r >= raypower
            return nothing
        end
        raydirection = direction(incidentray)
    else
        # do the diffraction
        order = rand((opticalinterface.minorder):(opticalinterface.maxorder)) # TODO not the most sensible way to sample this...
        refl = reflectance(opticalinterface, order) * internal_trans
        trans = transmission(opticalinterface, order) * internal_trans
        if r < refl
            raypower = refl * incident_pow
            normal = -normal
        elseif r < trans
            raypower = trans * incident_pow
        else
            return nothing
        end
        incident_mag = 2π * nᵢ / λ
        grating_dir = order * opticalinterface.vector * (2π / opticalinterface.period)
        u = cross(normal, ((direction(incidentray) * incident_mag) + grating_dir))
        output_mag = 2π * nₜ / λ
        q = output_mag^2 - sum(u .^ 2)
        if q < zero(T)
            @warn "Order $order not valid for this configuration, skipping this ray" maxlog = 1
            return nothing
        else
            raydirection = normalize(normal * sqrt(q) - cross(normal, u))
        end
    end
    return raydirection, raypower, raypathlength
end

######################################################################

"""
    HologramSurface{T,S} <: WrapperSurface{T,S}

Surface type for use with [`HologramInterface`](@ref).

```julia
HologramSurface(surface::Surface{T}, interface::HologramInterface{T})
```
"""
struct HologramSurface{T,S} <: WrapperSurface{T,S}
    surface::S
    interface::HologramInterface{T}
    HologramSurface(surface::S, interface::HologramInterface{T}) where {T<:Real,S<:Surface{T}} = new{T,S}(surface, interface)
end
export HologramSurface

"""
    MultiHologramSurface{T,S} <: WrapperSurface{T,S}

Surface type for use with [`MultiHologramInterface`](@ref).

```julia
MultiHologramSurface(surface::Surface{T}, interface::MultiHologramInterface{T})
```
"""
struct MultiHologramSurface{T,S} <: WrapperSurface{T,S}
    surface::S
    interface::MultiHologramInterface{T}
    MultiHologramSurface(surface::S, interface::MultiHologramInterface{T}) where {T<:Real,S<:Surface{T}} = new{T,S}(surface, interface)
end
export MultiHologramSurface

function processintersection(opticalinterface::HologramInterface{T}, point::SVector{N,T}, normal::SVector{N,T}, incidentray::OpticalRay{T,N}, temperature::T, pressure::T, test::Bool, firstray::Bool = false) where {T<:Real,N}
    hitback = dot(direction(incidentray), normal) > zero(T)
    # we want the surface to work if hit from either side, so we just reverse the normal and interfaces if it is hit on the back side
    if hitback
        frontinterface = FresnelInterface{T}(opticalinterface.substratematerial, opticalinterface.aftermaterial)
        backinterface = FresnelInterface{T}(opticalinterface.substratematerial, opticalinterface.beforematerial)
        normal = -normal
    else
        frontinterface = FresnelInterface{T}(opticalinterface.substratematerial, opticalinterface.beforematerial)
        backinterface = FresnelInterface{T}(opticalinterface.substratematerial, opticalinterface.aftermaterial)
    end

    # refract the incoming ray going from beforematerial to substrate material, can have Fresnel/TI reflections
    tmp = processintersection(frontinterface, point, normal, incidentray, temperature, pressure, test, firstray)
    if tmp === nothing
        return nothing
    end
    input_ray, input_power, input_opl = tmp
    if dot(input_ray, normal) > zero(T)
        # just reflected of the surface through fresenel reflection or TIR, so return without diffraction
        return input_ray, input_power, input_opl
    end

    λ = wavelength(incidentray)
    mat = glassforid(opticalinterface.substratematerial)::OpticSim.GlassCat.Glass
    # get the index of the playback ray in the substrate
    nₛ = index(mat, λ, temperature = temperature, pressure = pressure)::T
    # get the index of the recording ray in the substrate
    nₛrec = index(mat, opticalinterface.recordingλ, temperature = temperature, pressure = pressure)::T

    # get the index of the recording ray in the material of the signal beam
    if !isair(opticalinterface.signalrecordingmaterial)
        signalbeammaterial = glassforid(opticalinterface.signalrecordingmaterial)::OpticSim.GlassCat.Glass
        nsig = index(signalbeammaterial, opticalinterface.recordingλ, temperature = temperature, pressure = pressure)::T
    else
        nsig = one(T)
    end
    # get the signal magnitue in the substrate
    mag_Ksig = 2π * nₛ / opticalinterface.recordingλ
    # and direction
    if opticalinterface.signalbeamstate === CollimatedBeam
        Ksig_raw = opticalinterface.signalpointordir
    elseif opticalinterface.signalbeamstate === ConvergingBeam
        Ksig_raw = normalize(opticalinterface.signalpointordir - point)
    elseif opticalinterface.signalbeamstate === DivergingBeam
        Ksig_raw = normalize(point - opticalinterface.signalpointordir)
    else
        throw(ErrorException("Invalide beam state"))
    end
    # refract it for the substrate interface
    sigonside = dot(Ksig_raw, normal) < zero(T)
    refr = refractedray(nsig, nₛrec, sigonside ? -normal : normal, Ksig_raw) # TODO not sure normal flip is needed
    if refr === nothing
        return nothing
    end
    Ksig_refracted = mag_Ksig * refr

    # get the index of the recording ray in the material of the reference beam
    if !isair(opticalinterface.referencerecordingmaterial)
        referencebeammaterial = glassforid(opticalinterface.referencerecordingmaterial)::OpticSim.GlassCat.Glass
        nref = index(referencebeammaterial, opticalinterface.recordingλ, temperature = temperature, pressure = pressure)::T
    else
        nref = one(T)
    end
    # get the reference magnitue in the substrate
    mag_Kref = 2π * nₛ / opticalinterface.recordingλ
    if opticalinterface.referencebeamstate === CollimatedBeam
        Kref_raw = opticalinterface.referencepointordir
    elseif opticalinterface.referencebeamstate === ConvergingBeam
        Kref_raw = normalize(opticalinterface.referencepointordir - point)
    elseif opticalinterface.referencebeamstate === DivergingBeam
        Kref_raw = normalize(point - opticalinterface.referencepointordir)
    else
        throw(ErrorException("Invalide beam state"))
    end
    # refract it for the substrate interface
    refonside = dot(Kref_raw, normal) < zero(T)
    refr = refractedray(nref, nₛrec, refonside ? -normal : normal, Kref_raw) # TODO not sure normal flip is needed
    if refr === nothing
        return nothing
    end
    Kref_refracted = mag_Kref * refr

    # check whether we should be doing reflection or transmission
    isreflection = sigonside != refonside
    if isreflection
        # if the vectors are opposite in relation to the normal then the HOE should work in reflection
        dnormal = normal
    else
        # otherwise it is transmission
        dnormal = -normal
    end

    # calulcate the grating vector
    grating_vec = (Ksig_refracted - Kref_refracted)
    # choose random order for which Kogelnik is valid
    if !opticalinterface.include0order
        order = 1
    else
        order = rand(0:1)
    end

    # calculate the magnitudes of input and output
    incident_mag = 2π * nₛ / λ
    incident_ray = incident_mag * input_ray
    # do the diffraction
    u = cross(dnormal, (incident_ray + order * grating_vec))
    output_mag = 2π * nₛ / λ
    q = output_mag^2 - sum(u .^ 2)
    if q < zero(T)
        @warn "Order $order not valid for this configuration, skipping this ray" maxlog = 1
        return nothing
    else
        diffracted_out = dnormal * sqrt(q) - cross(dnormal, u)
    end
    diffracted_out = normalize(diffracted_out)

    # calculate the efficiency in this case
    d = opticalinterface.thickness # in microns as λ in microns
    Δn = opticalinterface.RImodulation
    σ = incident_ray + grating_vec
    cₒ = dot(σ, -normal) / incident_mag
    cᵣ = dot(incident_ray, -normal) / incident_mag
    υ = (norm(incident_ray)^2 - norm(σ)^2) / (2 * incident_mag) # Kick from Kogelnik eq. 17
    ξ = (υ * d) / (2 * cₒ) # Kogelnik eq. 42
    ν = (π * Δn * d) / (λ * sqrt(complex(cₒ * cᵣ))) # Kogelnik eq. 42
    refl = dot(diffracted_out, normal) > zero(T)
    if refl # reflected
        # modulate the power according to Kogelnik - Lossless Dielectric Grating (Reflection)
        ν = 1im * ν # Kogelnik eq. 55
        ξ = -ξ # Kogelnik eq. 55
        η = 1 / (1 + (1 - (ξ / ν)^2) / (sinh(sqrt(ν^2 - ξ^2))^2)) # Kogelnik eq. 57, Kick eq. 6
    else # refracted
        # modulate the power according to Kogelnik - Lossless Dielectric Grating (Transmission)
        η = sin(sqrt(ν^2 + ξ^2))^2 / (1 + (ξ / ν)^2) # Kogelnik eq. 43, Kick eq. 4
    end

    # kill the ray based on efficiency for monte carlo
    if !(imag(η) == zero(T) && real(η) > zero(T))
        @warn "Invalid η: $η"
        return nothing
    else
        if order === 0
            k = (1 - real(η))
            if !test * rand() >= k
                return nothing
            end
            diffracted_power = input_power * k
        else
            k = real(η)
            if !test * rand() >= k
                return nothing
            end
            diffracted_power = input_power * k
        end
    end
    outray = OpticalRay(point, diffracted_out, diffracted_power, λ)
    # refract the output ray for the substrate/aftermaterial interface
    # for refraction back interface is the other way around so reverse normal
    tmp = processintersection(refl ? frontinterface : backinterface, point, refl ? normal : -normal, outray, temperature, pressure, test, true)
    if tmp === nothing
        return nothing
    end
    outpur_dir, output_pow, _ = tmp
    # TODO nothing happening to OPL currently
    return outpur_dir, output_pow, input_opl
end

function processintersection(opticalinterface::MultiHologramInterface{T}, point::SVector{N,T}, normal::SVector{N,T}, incidentray::OpticalRay{T,N}, temperature::T, pressure::T, test::Bool, firstray::Bool = false) where {T<:Real,N}
    # minη = typemax(T)
    # r = rand(T)
    # Ση = zero(T)
    # for i in 1:(opticalinterface.numinterfaces)
    #     int = interface(opticalinterface, i)::HologramInterface{T}
    #     tmp = processintersection(int, point, normal, incidentray, temperature, pressure, scatter, firstray)
    #     if tmp !== nothing
    #         dir, pow, opl = tmp
    #         η = pow / power(incidentray)
    #         Ση += η / opticalinterface.numinterfaces
    #         # @assert Ση <= 1.0
    #         if r < Ση
    #             # don't modulate the power because we sample the efficiencies proportionally
    #             return dir, power(incidentray), opl
    #         end
    #     end
    # end
    # return nothing
    # when test = true behavior is totally invalid so that the test is deteministic...
    i = test ? 1 : rand(1:(opticalinterface.numinterfaces)) # TODO definitely not the best way to sample this...
    int = interface(opticalinterface, i)::HologramInterface{T}
    return processintersection(int, point, normal, incidentray, temperature, pressure, test, firstray)
end
