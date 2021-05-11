# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

using .OpticSim.GlassCat: AbstractGlass, Air

export SphericalLens, ConicLens, AsphericLens, FresnelLens

function opticinterface(
    ::Type{T},
    insidematerial::AbstractGlass,
    outsidematerial::AbstractGlass,
    reflectance::T = zero(T),
    interfacemode::InterfaceMode = ReflectOrTransmit
) where {T<:Real}
    transmission = one(T) - reflectance
    return FresnelInterface{T}(insidematerial, outsidematerial; reflectance, transmission, interfacemode)
end

"""
    SphericalLens(
        insidematerial,
        frontvertex,
        frontradius,
        backradius,
        thickness,
        semidiameter;
        lastmaterial = Air,
        nextmaterial = Air,
        frontsurfacereflectance = 0.0,
        backsurfacereflectance = 0.0,
        frontdecenter = (0, 0),
        backdecenter = (0, 0),
        interfacemode = ReflectOrTransmit
    )

Constructs a simple cylindrical lens with spherical front and back surfaces. The side walls of the lens are absorbing.
"""
function SphericalLens(
    insidematerial::AbstractGlass,
    frontvertex::T,
    frontradius::T,
    backradius::T,
    thickness::T,
    semidiameter::T;
    lastmaterial::AbstractGlass = Air,
    nextmaterial::AbstractGlass = Air,
    frontsurfacereflectance::T = zero(T),
    backsurfacereflectance::T = zero(T),
    frontdecenter::Tuple{T,T} = (zero(T), zero(T)),
    backdecenter::Tuple{T,T} = (zero(T), zero(T)),
    interfacemode::InterfaceMode = ReflectOrTransmit
) where {T<:Real}
    @assert !isnan(frontradius)
    @assert !isnan(backradius)

    return ConicLens(
        insidematerial,
        frontvertex,
        frontradius,
        zero(T),
        backradius,
        zero(T),
        thickness,
        semidiameter;
        lastmaterial,
        nextmaterial,
        frontsurfacereflectance,
        backsurfacereflectance,
        frontdecenter,
        backdecenter,
        interfacemode
    )
end

"""
    ConicLens(
        insidematerial,
        frontvertex,
        frontradius,
        frontconic,
        backradius,
        backconic,
        thickness,
        semidiameter;
        lastmaterial = Air,
        nextmaterial = Air,
        frontsurfacereflectance = 0.0,
        backsurfacereflectance = 0.0,
        frontdecenter = (0, 0),
        backdecenter = (0, 0),
        interfacemode = ReflectOrTransmit
    )

Constructs a simple cylindrical lens with front and back surfaces with a radius and conic term. The side walls of the
lens are absorbing.
"""
function ConicLens(
    insidematerial::AbstractGlass,
    frontvertex::T,
    frontradius::T,
    frontconic::T,
    backradius::T,
    backconic::T,
    thickness::T,
    semidiameter::T;
    lastmaterial::AbstractGlass = Air,
    nextmaterial::AbstractGlass = OpticSim.GlassCat.Air,
    frontsurfacereflectance::T = zero(T),
    backsurfacereflectance::T = zero(T),
    nsamples::Int = 17,
    frontdecenter::Tuple{T,T} = (zero(T), zero(T)),
    backdecenter::Tuple{T,T} = (zero(T), zero(T)),
    interfacemode::InterfaceMode = ReflectOrTransmit
) where {T<:Real}
    @assert !isnan(frontradius)
    @assert !isnan(frontconic)

    return AsphericLens(
        insidematerial,
        frontvertex,
        frontradius,
        frontconic,
        nothing,
        backradius,
        backconic,
        nothing,
        thickness,
        semidiameter;
        lastmaterial,
        nextmaterial,
        frontsurfacereflectance,
        backsurfacereflectance,
        nsamples,
        frontdecenter,
        backdecenter,
        interfacemode
    )
end

"""
    AsphericLens(
        insidematerial,
        frontvertex,
        frontradius,
        frontconic,
        frontaspherics,
        backradius,
        backconic,
        backaspherics,
        thickness,
        semidiameter;
        lastmaterial = Air,
        nextmaterial = Air,
        frontsurfacereflectance = 0.0,
        backsurfacereflectance = 0.0,
        frontdecenter = (0, 0),
        backdecenter = (0, 0),
        interfacemode = ReflectOrTransmit
    )

Cosntructs a simple cylindrical lens with front and back surfaces with a radius, conic and apsheric terms. The side
walls of the lens are absorbing.
"""
function AsphericLens(
    insidematerial::AbstractGlass,
    frontvertex::T,
    frontradius::T,
    frontconic::T,
    frontaspherics::Union{Nothing,Vector{Pair{Int,T}}},
    backradius::T,
    backconic::T,
    backaspherics::Union{Nothing,Vector{Pair{Int,T}}},
    thickness::T,
    semidiameter::T;
    lastmaterial::AbstractGlass = Air,
    nextmaterial::AbstractGlass = Air,
    frontsurfacereflectance::T = zero(T),
    backsurfacereflectance::T = zero(T),
    nsamples::Int = 17,
    frontdecenter::Tuple{T,T} = (zero(T), zero(T)),
    backdecenter::Tuple{T,T} = (zero(T), zero(T)),
    interfacemode::InterfaceMode = ReflectOrTransmit
) where {T<:Real}
    @assert semidiameter > zero(T)
    @assert !isnan(frontradius)

    frontdecenter_l = (
        abs(frontdecenter[1]) > eps(T) && abs(frontdecenter[2]) > eps(T) ?
        sqrt(frontdecenter[1]^2 + frontdecenter[2]^2) :
        zero(T)
    )
    backdecenter_l = (
        abs(backdecenter[1]) > eps(T) && abs(backdecenter[2]) > eps(T) ?
        sqrt(backdecenter[1]^2 + backdecenter[2]^2) :
        zero(T)
    )

    # front lens
    frontinterface = opticinterface(T, insidematerial, lastmaterial, frontsurfacereflectance, interfacemode)
    if isinf(frontradius) && frontaspherics === nothing
        # planar
        lens_front = leaf(Plane(
            SVector{3,T}(0, 0, 1),
            SVector{3,T}(0, 0, frontvertex);
            vishalfsizeu = semidiameter,
            vishalfsizev = semidiameter,
            interface = frontinterface
        ))
    elseif frontconic == zero(T) && frontaspherics === nothing
        # spherical
        if frontradius < zero(T)
            ϕmax = NaNsafeasin(min(
                abs((min(semidiameter, abs(frontradius)) + frontdecenter_l) / frontradius),
                one(T)
            )) + convert(T, π / 50)
            lens_front = leaf(SphericalCap(
                abs(frontradius),
                ϕmax,
                SVector{3,T}(0, 0, 1),
                SVector{3,T}(frontdecenter[1], frontdecenter[2], frontvertex),
                interface = frontinterface
            ))
        else
            offset = translation(T, frontdecenter[1], frontdecenter[2], frontvertex - frontradius)
            surf = Sphere(abs(frontradius), interface = frontinterface)
            lens_front = leaf(surf, offset)
        end
        if abs(frontradius) - frontdecenter_l <= semidiameter
            # if optical surface is smaller than semidiameter then use a plane to fill the gap
            plane = leaf(Plane(
                SVector{3,T}(0, 0, 1),
                SVector{3,T}(0, 0, frontvertex - frontradius),
                vishalfsizeu = semidiameter,
                vishalfsizev = semidiameter,
                interface = frontinterface
            ))
            if frontradius > zero(S)
                lens_front = lens_front ∪ plane
            else
                lens_front = lens_front ∩ plane
            end
            lens_front = frontradius > zero(T) ? csgunion(lens_front, plane) : csgintersection(lens_front, plane)
        end
    else
        # conic or aspheric
        surf = AcceleratedParametricSurface(
            ZernikeSurface(
                semidiameter + frontdecenter_l + convert(T, 0.01);
                radius = -frontradius,
                conic = frontconic,
                aspherics = frontaspherics === nothing ? nothing : [(i, -k) for (i, k) in frontaspherics]
            ),
            nsamples,
            interface = frontinterface
        )
        lens_front = leaf(surf, translation(T, frontdecenter[1], frontdecenter[2], frontvertex))
    end

    # back lens
    backinterface = opticinterface(T, insidematerial, nextmaterial, backsurfacereflectance, interfacemode)
    if isinf(backradius) && backaspherics === nothing
        # planar
        lens_rear = leaf(Plane(
            SVector{3,T}(0, 0, -1),
            SVector{3,T}(0, 0, frontvertex - thickness),
            vishalfsizeu = semidiameter,
            vishalfsizev = semidiameter,
            interface = backinterface
        ))
    elseif backconic == zero(T) && backaspherics === nothing
        # spherical
        if backradius < zero(T)
            offset = translation(T, backdecenter[1], backdecenter[2], (frontvertex - thickness) - backradius)
            surf = Sphere(abs(backradius), interface = backinterface)
            lens_rear = leaf(surf, offset)
            if backradius == semidiameter
                # in this case we need to clip the sphere
                plane = leaf(Plane(SVector{3,S}(0, 0, 1), SVector{3,S}(0, 0, frontvertex - thickness + backradius)))
                lens_rear = lens_rear ∩ plane
            end
        else
            ϕmax = NaNsafeasin(min(
                abs((min(semidiameter, abs(backradius)) + backdecenter_l) / backradius),
                one(T)
            )) + convert(T, π / 50)
            lens_rear = leaf(SphericalCap(
                abs(backradius),
                ϕmax,
                SVector{3,T}(0, 0, -1),
                SVector{3,T}(backdecenter[1], backdecenter[2], frontvertex - thickness),
                interface = backinterface
            ))
        end
        if abs(backradius) - backdecenter_l <= semidiameter
            plane = leaf(Plane(
                SVector{3,T}(0, 0, -1),
                SVector{3,T}(0, 0, frontvertex - thickness - backradius),
                vishalfsizeu = semidiameter,
                vishalfsizev = semidiameter,
                interface = backinterface
            ))
            if backradius < zero(S)
                lens_rear = lens_rear ∪ plane
            else
                lens_rear = lens_rear ∩ plane
            end
        end
    else
        #conic or aspheric
        surf = AcceleratedParametricSurface(
            ZernikeSurface(
                semidiameter + backdecenter_l + convert(T, 0.01);
                radius = backradius,
                conic = backconic,
                aspherics = backaspherics === nothing ? nothing : Tuple{Int,T}.(backaspherics)
            ),
            nsamples,
            interface = backinterface
        )
        lens_rear = leaf(
            surf,
            Transform{T}(zero(T), convert(T, π), zero(T), backdecenter[1], backdecenter[2], frontvertex - thickness)
        )
    end

    extra_front = (
        frontradius >= zero(T) || isinf(frontradius) ?
        zero(T) :
        abs(frontradius) - sqrt(frontradius^2 - semidiameter^2)
    )
    extra_back = (
        backradius >= zero(T) || isinf(backradius) ?
        zero(T) :
        abs(backradius) - sqrt(backradius^2 - semidiameter^2)
    )
    barrel_center = ((frontvertex + extra_front) + (frontvertex - thickness - extra_back)) / 2
    lens_barrel = leaf(
        Cylinder(
            semidiameter,
            (thickness + extra_back + extra_front) * 2,
            interface = FresnelInterface{T}(insidematerial, Air; reflectance = zero(T), transmission = zero(T))
        ),
        translation(T, zero(T), zero(T), barrel_center)
    )
    lens_csg = lens_front ∩ lens_rear ∩ lens_barrel
    return lens_csg
end

"""
    FresnelLens(
        insidematerial,
        frontvertex,
        radius,
        thickness,
        semidiameter,
        groovedepth;
        conic = 0.0,
        aspherics = nothing,
        outsidematerial = Air
    )

Create a Fresnel lens as a CSG object, can be concave or convex. Groove positions are found iteratively based on
`groovedepth`. For negative radii the vertex on the central surface is at `frontvertex`, so the total thickness of the
lens is `thickness` + `groovedepth`.

**Aspherics currently not supported**.
"""
function FresnelLens(
    insidematerial::AbstractGlass,
    frontvertex::T,
    radius::T,
    thickness::T,
    semidiameter::T,
    groovedepth::T;
    conic::T = 0.0,
    aspherics::Union{Nothing,Vector{Tuple{Int,T}}} = nothing,
    outsidematerial::AbstractGlass = Air,
    reverse::Bool = false
) where {T<:Real}
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
        # spherical lens which makes things much easier, we can just caluclate groove positions directly
        n = 1
        if radius > zero(T)
            sphere = Sphere(radius; interface)
            fresnel = leaf(sphere, translation(T, zero(T), zero(T), frontvertex - radius))
        else
            ϕmax = NaNsafeasin(abs(semidiameter / radius)) + convert(T, π / 50)
            sphere = SphericalCap(abs(radius), ϕmax; interface)
            fresnel = leaf(sphere, translation(T, zero(T), zero(T), frontvertex))
        end
        while true
            offset = n * groovedepth
            cylrad = sqrt(offset * (2 * abs(radius) - offset))
            if cylrad >= semidiameter
                break
            end
            if radius > zero(T)
                cylinder = leaf(
                    Cylinder(cylrad, groovedepth * 1.5; interface),
                    translation(T, zero(T), zero(T), frontvertex - groovedepth / 2)
                )
                newsurf = leaf(sphere, translation(T, zero(T), zero(T), frontvertex - radius + offset)) - cylinder
                fresnel = newsurf ∪ fresnel
            else
                cylinder = leaf(
                    Cylinder(cylrad, groovedepth * 1.5; interface),
                    translation(T, zero(T), zero(T), frontvertex + groovedepth / 2)
                )
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
        surface = AcceleratedParametricSurface(
            ZernikeSurface(semidiameter + convert(T, 0.01); radius = -radius, conic);
            interface
        )
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
                cylinder = leaf(
                    Cylinder(cylrad, groovedepth * 1.5; interface),
                    translation(T, zero(T), zero(T), frontvertex - groovedepth / 2)
                )
                newsurf = leaf(surface, translation(T, zero(T), zero(T), frontvertex + offset)) - cylinder
                fresnel = newsurf ∪ fresnel
            else
                cylinder = leaf(
                    Cylinder(cylrad, groovedepth * 1.5; interface),
                    translation(T, zero(T), zero(T), frontvertex + groovedepth / 2)
                )
                newsurf = leaf(surface, translation(T, zero(T), zero(T), frontvertex + offset)) ∪ cylinder
                fresnel = newsurf ∩ fresnel
            end
            n += 1
        end
    else
        # TODO - CSG gets more complicated as the surface can be both convex and concave and we need to find the groove
        # positions iteratively
        @error "Ashperics not supported for Fresnel lenses"
    end

    outer_barrel = leaf(
        Cylinder(semidiameter, thickness * 2; interface),
        translation(T, zero(T), zero(T), frontvertex - thickness / 2)
    )
    fresnel = outer_barrel ∩ fresnel
    backplane = leaf(Plane(
        SVector{3,T}(0, 0, -1),
        SVector{3,T}(0, 0, frontvertex - thickness);
        vishalfsizeu = semidiameter,
        vishalfsizev = semidiameter
    ))
    fresnel = reverse ? backplane ∩ fresnel : leaf(fresnel, rotationd(T, zero(T), convert(T, 180), zero(T)))
    return fresnel
end
