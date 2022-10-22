# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

opticinterface(::Type{S}, insidematerial, outsidematerial, reflectance = zero(S), interfacemode = ReflectOrTransmit) where {S<:Real} = FresnelInterface{S}(insidematerial, outsidematerial, reflectance = reflectance, transmission = one(S) - reflectance, interfacemode = interfacemode)

"""
    SphericalLens(insidematerial, frontvertex, frontradius, backradius, thickness, semidiameter;  lastmaterial = OpticSim.GlassCat.Air, nextmaterial = OpticSim.GlassCat.Air, frontsurfacereflectance = 0.0, backsurfacereflectance = 0.0, frontdecenter = (0, 0), backdecenter = (0, 0), interfacemode = ReflectOrTransmit)

Constructs a simple cylindrical lens with spherical front and back surfaces. The side walls of the lens are absorbing.
"""
function SphericalLens(insidematerial::T, frontvertex::S, frontradius::S, backradius::S, thickness::S, semidiameter::S; lastmaterial::Q = OpticSim.GlassCat.Air, nextmaterial::R = OpticSim.GlassCat.Air, frontsurfacereflectance::S = zero(S), backsurfacereflectance::S = zero(S), frontdecenter::Tuple{S,S} = (zero(S), zero(S)), backdecenter::Tuple{S,S} = (zero(S), zero(S)), interfacemode = ReflectOrTransmit) where {R<:OpticSim.GlassCat.AbstractGlass,Q<:OpticSim.GlassCat.AbstractGlass,T<:OpticSim.GlassCat.AbstractGlass,S<:Real}
    @assert !isnan(frontradius)
    @assert !isnan(backradius)

    return ConicLens(insidematerial, frontvertex, frontradius, zero(S), backradius, zero(S), thickness, semidiameter, lastmaterial = lastmaterial, nextmaterial = nextmaterial, frontsurfacereflectance = frontsurfacereflectance, backsurfacereflectance = backsurfacereflectance, frontdecenter = frontdecenter, backdecenter = backdecenter, interfacemode = interfacemode)
end
export SphericalLens

"""
    ConicLens(insidematerial, frontvertex, frontradius, frontconic, backradius, backconic, thickness, semidiameter;  lastmaterial = OpticSim.GlassCat.Air, nextmaterial = OpticSim.GlassCat.Air, frontsurfacereflectance = 0.0, backsurfacereflectance = 0.0, frontdecenter = (0, 0), backdecenter = (0, 0), interfacemode = ReflectOrTransmit)

Constructs a simple cylindrical lens with front and back surfaces with a radius and conic term. The side walls of the lens are absorbing.
"""
function ConicLens(insidematerial::T, frontvertex::S, frontradius::S, frontconic::S, backradius::S, backconic::S, thickness::S, semidiameter::S; lastmaterial::Q = OpticSim.GlassCat.Air, nextmaterial::R = OpticSim.GlassCat.Air, frontsurfacereflectance::S = zero(S), backsurfacereflectance::S = zero(S), nsamples::Int = 17, frontdecenter::Tuple{S,S} = (zero(S), zero(S)), backdecenter::Tuple{S,S} = (zero(S), zero(S)), interfacemode = ReflectOrTransmit) where {R<:OpticSim.GlassCat.AbstractGlass,Q<:OpticSim.GlassCat.AbstractGlass,T<:OpticSim.GlassCat.AbstractGlass,S<:Real}
    @assert !isnan(frontradius)
    @assert !isnan(frontconic)

    return AsphericLens(insidematerial, frontvertex, frontradius, frontconic, nothing, backradius, backconic, nothing, thickness, semidiameter, lastmaterial = lastmaterial, nextmaterial = nextmaterial, frontsurfacereflectance = frontsurfacereflectance, backsurfacereflectance = backsurfacereflectance, nsamples = nsamples, frontdecenter = frontdecenter, backdecenter = backdecenter, interfacemode = interfacemode)
end
export ConicLens

"""
    AsphericLens(insidematerial, frontvertex, frontradius, frontconic, frontaspherics, backradius, backconic, backaspherics, thickness, semidiameter;  lastmaterial = OpticSim.GlassCat.Air, nextmaterial = OpticSim.GlassCat.Air, frontsurfacereflectance = 0.0, backsurfacereflectance = 0.0, frontdecenter = (0, 0), backdecenter = (0, 0), interfacemode = ReflectOrTransmit)

Cosntructs a simple cylindrical lens with front and back surfaces with a radius, conic and apsheric terms. The side walls of the lens are absorbing.
"""
function AsphericLens(insidematerial::T, frontvertex::S, frontradius::S, frontconic::S, frontaspherics::Union{Nothing,Vector{Pair{Int,S}}}, backradius::S, backconic::S, backaspherics::Union{Nothing,Vector{Pair{Int,S}}}, thickness::S, semidiameter::S; lastmaterial::Q = OpticSim.GlassCat.Air, nextmaterial::R = OpticSim.GlassCat.Air, frontsurfacereflectance::S = zero(S), backsurfacereflectance::S = zero(S), nsamples::Int = 17, frontdecenter::Tuple{S,S} = (zero(S), zero(S)), backdecenter::Tuple{S,S} = (zero(S), zero(S)), interfacemode = ReflectOrTransmit) where {R<:OpticSim.GlassCat.AbstractGlass,Q<:OpticSim.GlassCat.AbstractGlass,T<:OpticSim.GlassCat.AbstractGlass,S<:Real}
    @assert semidiameter > zero(S)
    @assert !isnan(frontradius)

    frontdecenter_l = abs(frontdecenter[1]) > eps(S) && abs(frontdecenter[2]) > eps(S) ? sqrt(frontdecenter[1]^2 + frontdecenter[2]^2) : zero(S)
    backdecenter_l = abs(backdecenter[1]) > eps(S) && abs(backdecenter[2]) > eps(S) ? sqrt(backdecenter[1]^2 + backdecenter[2]^2) : zero(S)

    if isinf(frontradius) && (frontaspherics === nothing)
        #planar
        lens_front = leaf(Plane(SVector{3,S}(0, 0, 1), SVector{3,S}(0, 0, frontvertex), vishalfsizeu = semidiameter, vishalfsizev = semidiameter, interface = interface = opticinterface(S, insidematerial, lastmaterial, frontsurfacereflectance, interfacemode)))
    elseif frontconic == zero(S) && (frontaspherics === nothing)
        #spherical
        if frontradius < zero(S)
            ϕmax = NaNsafeasin(min(abs((min(semidiameter, abs(frontradius)) + frontdecenter_l) / frontradius), one(S))) + S(π / 50)
            lens_front = leaf(SphericalCap(abs(frontradius), ϕmax, SVector{3,S}(0, 0, 1), SVector{3,S}(frontdecenter[1], frontdecenter[2], frontvertex), interface = opticinterface(S, insidematerial, lastmaterial, frontsurfacereflectance, interfacemode)))
        else
            offset = translation(S, frontdecenter[1], frontdecenter[2], frontvertex - frontradius)
            surf = Sphere(abs(frontradius), interface = opticinterface(S, insidematerial, lastmaterial, frontsurfacereflectance, interfacemode))
            lens_front = leaf(surf, offset)
        end
        if abs(frontradius) - frontdecenter_l <= semidiameter
            # if optical surface is smaller than semidiameter then use a plane to fill the gap
            plane = leaf(Plane(SVector{3,S}(0, 0, 1), SVector{3,S}(0, 0, frontvertex - frontradius), vishalfsizeu = semidiameter, vishalfsizev = semidiameter, interface = interface = opticinterface(S, insidematerial, lastmaterial, frontsurfacereflectance, interfacemode)))
            if frontradius > zero(S)
                lens_front = lens_front ∪ plane
            else
                lens_front = lens_front ∩ plane
            end
        end
    else
        #conic or aspheric
        if frontaspherics !== nothing
            frontaspherics = [(i, -k) for (i, k) in frontaspherics]
        end
        surf = AcceleratedParametricSurface(AsphericSurface(semidiameter + frontdecenter_l + S(0.01), radius = -frontradius, conic = frontconic, aspherics = frontaspherics), nsamples, interface = opticinterface(S, insidematerial, lastmaterial, frontsurfacereflectance, interfacemode))
        lens_front = leaf(surf, translation(S, frontdecenter[1], frontdecenter[2], frontvertex))
    end
    if isinf(backradius) && (backaspherics === nothing)
        #planar
        lens_rear = leaf(Plane(SVector{3,S}(0, 0, -1), SVector{3,S}(0, 0, frontvertex - thickness), vishalfsizeu = semidiameter, vishalfsizev = semidiameter, interface = interface = opticinterface(S, insidematerial, nextmaterial, backsurfacereflectance, interfacemode)))
    elseif backconic == zero(S) && (backaspherics === nothing)
        #spherical
        if backradius < zero(S)
            offset = translation(S, backdecenter[1], backdecenter[2], (frontvertex - thickness) - backradius)
            surf = Sphere(abs(backradius), interface = opticinterface(S, insidematerial, nextmaterial, backsurfacereflectance, interfacemode))
            lens_rear = leaf(surf, offset)
            if backradius == semidiameter
                # in this case we need to clip the sphere
                plane = leaf(Plane(SVector{3,S}(0, 0, 1), SVector{3,S}(0, 0, frontvertex - thickness + backradius)))
                lens_rear = lens_rear ∩ plane
            end
        else
            ϕmax = NaNsafeasin(min(abs((min(semidiameter, abs(backradius)) + backdecenter_l) / backradius), one(S))) + S(π / 50)
            lens_rear = leaf(SphericalCap(abs(backradius), ϕmax, SVector{3,S}(0, 0, -1), SVector{3,S}(backdecenter[1], backdecenter[2], frontvertex - thickness), interface = opticinterface(S, insidematerial, nextmaterial, backsurfacereflectance, interfacemode)))
        end
        if abs(backradius) - backdecenter_l <= semidiameter
            plane = leaf(Plane(SVector{3,S}(0, 0, -1), SVector{3,S}(0, 0, frontvertex - thickness - backradius), vishalfsizeu = semidiameter, vishalfsizev = semidiameter, interface = interface = opticinterface(S, insidematerial, lastmaterial, backsurfacereflectance, interfacemode)))
            if backradius < zero(S)
                lens_rear = lens_rear ∪ plane
            else
                lens_rear = lens_rear ∩ plane
            end
        end
    else
        #conic or aspheric
        if backaspherics !== nothing
            backaspherics = Tuple{Int,S}.(backaspherics)
        end
        surf = AcceleratedParametricSurface(AsphericSurface(semidiameter + backdecenter_l + S(0.01), radius = backradius, conic = backconic, aspherics = backaspherics), nsamples, interface = opticinterface(S, insidematerial, nextmaterial, backsurfacereflectance, interfacemode))
        lens_rear = leaf(surf, Transform(zero(S), S(π), zero(S), backdecenter[1], backdecenter[2], frontvertex - thickness))
    end
    extra_front = frontradius >= zero(S) || isinf(frontradius) ? zero(S) : abs(frontradius) - sqrt(frontradius^2 - semidiameter^2)
    extra_back = backradius >= zero(S) || isinf(backradius) ? zero(S) : abs(backradius) - sqrt(backradius^2 - semidiameter^2)
    barrel_center = ((frontvertex + extra_front) + (frontvertex - thickness - extra_back)) / 2
    lens_barrel = leaf(Cylinder(semidiameter, (thickness + extra_back + extra_front) * 2, interface = FresnelInterface{S}(insidematerial, OpticSim.GlassCat.Air, reflectance = zero(S), transmission = zero(S))), translation(S, zero(S), zero(S), barrel_center))
    lens_csg = lens_front ∩ lens_rear ∩ lens_barrel
    return lens_csg
end
export AsphericLens

"""
    FresnelLens(insidematerial, frontvertex, radius, thickness, semidiameter, groovedepth; conic = 0.0, aspherics = nothing, outsidematerial = OpticSim.GlassCat.Air)

Create a Fresnel lens as a CSG object, can be concave or convex. Groove positions are found iteratively based on `groovedepth`. For negative radii the vertex on the central surface is at `frontvertex`, so the total thickness of the lens is `thickness` + `groovedepth`.
**Aspherics currently not supported**.
"""
function FresnelLens(insidematerial::G, frontvertex::T, radius::T, thickness::T, semidiameter::T, groovedepth::T; conic::T = 0.0, aspherics::Union{Nothing,Vector{Pair{Int,T}}} = nothing, outsidematerial::H = OpticSim.GlassCat.Air, reverse::Bool = false) where {T<:Real,G<:OpticSim.GlassCat.AbstractGlass,H<:OpticSim.GlassCat.AbstractGlass}
    @assert abs(radius) > semidiameter
    interface = FresnelInterface{T}(insidematerial, outsidematerial)

    if aspherics !== nothing
        aspherics = [(i, -k) for (i, k) in aspherics]
    end

    # make the fresnel
    if conic == zero(T)
        if radius == zero(T) || isinf(radius)
            @error "Invalid radius"
        end
        # spherical lens which makes things much easier, we can just calculate groove positions directly
        n = 1
        if radius > zero(T)
            sphere = Sphere(radius, interface = interface)
            fresnel = leaf(sphere, translation(T, zero(T), zero(T), frontvertex - radius))
        else
            ϕmax = NaNsafeasin(abs(semidiameter / radius)) + T(π / 50)
            sphere = SphericalCap(abs(radius), ϕmax, interface = interface)
            fresnel = leaf(sphere, translation(T, zero(T), zero(T), frontvertex))
        end
        while true
            offset = n * groovedepth
            cylrad = sqrt(offset * (2 * abs(radius) - offset))
            if cylrad >= semidiameter
                break
            end
            if radius > zero(T)
                cylinder = leaf(Cylinder(cylrad, groovedepth * 1.5, interface = interface), translation(T, zero(T), zero(T), frontvertex - groovedepth / 2))
                newsurf = leaf(sphere, translation(T, zero(T), zero(T), frontvertex - radius + offset)) - cylinder
                fresnel = newsurf ∪ fresnel
            else
                cylinder = leaf(Cylinder(cylrad, groovedepth * 1.5, interface = interface), translation(T, zero(T), zero(T), frontvertex + groovedepth / 2))
                newsurf = leaf(sphere, translation(T, zero(T), zero(T), frontvertex - offset)) ∪ cylinder
                fresnel = newsurf ∩ fresnel
            end
            n += 1
        end
    elseif (aspherics === nothing)
        if radius == zero(T) || isinf(radius)
            @error "Invalid radius"
        end
        n = 1
        surface = AcceleratedParametricSurface(ZernikeSurface(semidiameter + T(0.01), radius = -radius, conic = conic), interface = interface)
        fresnel = leaf(surface, translation(T, zero(T), zero(T), frontvertex))
        while true
            offset = n * groovedepth
            if radius < zero(T)
                offset *= -1
            end
            cylrad = sqrt(2 * radius * offset - (conic + one(T)) * offset^2)
            if cylrad >= semidiameter
                break
            end
            if radius > zero(T)
                cylinder = leaf(Cylinder(cylrad, groovedepth * 1.5, interface = interface), translation(T, zero(T), zero(T), frontvertex - groovedepth / 2))
                newsurf = leaf(surface, translation(T, zero(T), zero(T), frontvertex + offset)) - cylinder
                fresnel = newsurf ∪ fresnel
            else
                cylinder = leaf(Cylinder(cylrad, groovedepth * 1.5, interface = interface), translation(T, zero(T), zero(T), frontvertex + groovedepth / 2))
                newsurf = leaf(surface, translation(T, zero(T), zero(T), frontvertex + offset)) ∪ cylinder
                fresnel = newsurf ∩ fresnel
            end
            n += 1
        end
    else
        # TODO - CSG gets more complicated as the surface can be both convex and concave and we need to find the groove positions iteratively
        @error "Ashperics not supported for Fresnel lenses"
    end

    outer_barrel = leaf(Cylinder(semidiameter, thickness * 2, interface = interface), translation(T, zero(T), zero(T), frontvertex - thickness / 2))
    fresnel = outer_barrel ∩ fresnel
    backplane = leaf(Plane(zero(T), zero(T), -one(T), zero(T), zero(T), frontvertex - thickness, vishalfsizeu = semidiameter, vishalfsizev = semidiameter))
    fresnel = backplane ∩ fresnel
    if !reverse
        fresnel = leaf(fresnel, rotationd(T, zero(T), T(180.0), zero(T)))
    end
    return fresnel
end
export FresnelLens
