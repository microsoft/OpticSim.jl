# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.



#########################################################################################################



"""
    AsphericSurface{T,N,Q} <: ParametricSurface{T,N}

Surface incorporating an aspheric polynomial - radius, conic and aspherics are defined relative to absolute semi-diameter,.
`T` is the datatype, `N` is the dimensionality,  and `Q` is the number of aspheric terms. 

The surface is centered at the origin and treated as being the cap of an infinite cylinder, thus creating a true half-space.
Outside of `0 <= ρ <= 1` the height of the surface is not necessarily well defined, so NaN may be returned.


```julia
AsphericSurface(semidiameter; radius, conic, aspherics=nothing, normradius = semidiameter)
```

`aspherics` should be a vector containing the aspheric polynomial coefficients starting at A1

The sag is defined by the equation

```math
z(r,\\phi) = \\frac{cr^2}{1 + \\sqrt{1 - (1+k)c^2r^2}} + \\sum_{i}^{Q}\\alpha_ir^{i}
```

where ``\\rho = \\frac{r}{\\texttt{normradius}}``, ``c = \\frac{1}{\\texttt{radius}}``, ``k = \\texttt{conic}`` .
"""
struct AsphericSurface{T,N,Q} <: ParametricSurface{T,N}
    semidiameter::T
    aspherics::SVector{Q,T}
    curvature::T
    conic::T
    normradius::T
    boundingcylinder::Cylinder{T,N}

    function AsphericSurface(semidiameter::T; radius::T = typemax(T), conic::T= zero(T), aspherics::Union{Nothing,Vector{Tuple{Int,T}}} = nothing, normradius::T = semidiameter) where {T<:Real}
        @assert semidiameter > 0
        @assert !isnan(semidiameter) && !isnan(radius) && !isnan(conic)
        @assert one(T) - (1 / radius)^2 * (conic + one(T)) * semidiameter^2 > 0 "Invalid surface (conic/radius combination: $radius, $conic)"

        acs = []
        if aspherics !== nothing
            asphericTerms = [i for (i, ) in aspherics]
            minAsphericTerm = minimum(asphericTerms)
            maxAsphericTerm = maximum(asphericTerms)
            @assert minAsphericTerm > 0 "Aspheric Terms must be A1 or higher (A$minAsphericTerm)"
            acs = zeros(T, maxAsphericTerm )
            for (i, k) in aspherics
                acs[i] = k
            end
        end
        Q = length(acs)
        new{T,3,Q}(semidiameter, SVector{Q,T}(acs), 1 / radius, conic, normradius, Cylinder(semidiameter, interface = opaqueinterface(T))) # TODO!! incorrect interface on cylinder
    end

end
export AsphericSurface

uvrange(::Type{AsphericSurface{T,N,Q}}) where {T<:Real,N,Q} = ((zero(T), one(T)), (-T(π), T(π))) # ρ and ϕ

semidiameter(z::AsphericSurface{T}) where {T<:Real} = z.semidiameter
halfsizeu(z::AsphericSurface{T}) where {T<:Real} = semidiameter(z)
halfsizev(z::AsphericSurface{T}) where {T<:Real} = semidiameter(z)

boundingobj(z::AsphericSurface{T}) where {T<:Real} = z.boundingcylinder

function point(z::AsphericSurface{T,3,Q}, ρ::T, ϕ::T)::SVector{3,T} where {T<:Real,Q}
    rad = z.semidiameter
    r = ρ * rad
    r2 = r^2
    t = one(T) - z.curvature^2 * (z.conic + one(T)) * r^2
    if t < zero(T)
        return SVector{3,T}(NaN, NaN, NaN)
    end
    h = z.curvature * r2 / (one(T) + sqrt(t))
    # sum aspheric
    prod = one(T)
    for asp in z.aspherics
        prod *= r
        h += asp * prod
    end
    return SVector{3,T}(r * cos(ϕ), r * sin(ϕ), h)
end

function partials(z::AsphericSurface{T,3,Q}, ρ::T, ϕ::T)::Tuple{SVector{3,T},SVector{3,T}} where {T<:Real,Q}
    rad = z.semidiameter
    r = ρ * rad
    t = one(T) - z.curvature^2 * (z.conic + one(T)) * r^2
    if t < zero(T)
        return SVector{3,T}(NaN, NaN, NaN), SVector{3,T}(NaN, NaN, NaN)
    end
    dhdρ = rad * z.curvature * r * sqrt(t) / t
    # sum aspherics partial
    if length(z.aspherics) != 0
        ((m, asp), rest) = Iterators.peel(enumerate(z.aspherics))
        dhdρ += rad * asp  #first term m=1 and prod = one(T)
        prod = one(T)
        for (m,asp) in rest
            prod *= r 
            dhdρ += rad * m * asp * prod 
        end
    end
    dhdϕ = zero(T)
    cosϕ = cos(ϕ)
    sinϕ = sin(ϕ)
    pu = SVector{3,T}(rad * cosϕ, rad * sinϕ, dhdρ)
    pv = SVector{3,T}(r * -sinϕ, r * cosϕ, dhdϕ)
    return pu, pv
end

function normal(z::AsphericSurface{T,3,Q}, ρ::T, ϕ::T)::SVector{3,T} where {T<:Real,Q}
    du, dv = partials(z, ρ, ϕ)
    if ρ == zero(T) && norm(dv) == zero(T)
        # in cases where there is no δϕ at ρ = 0 (i.e. anything which is rotationally symetric)
        # then we get some big problems, hardcoding this case solves the problems
        return SVector{3,T}(0, 0, 1)
    end
    return normalize(cross(du, dv))
end

function uv(z::AsphericSurface{T,3,Q}, p::SVector{3,T}) where {T<:Real,Q}
    # avoid divide by zero for ForwardDiff
    ϕ = NaNsafeatan(p[2], p[1])
    if p[1] == zero(T) && p[2] == zero(T)
        ρ = zero(T)
    else
        ρ = sqrt(p[1]^2 + p[2]^2) / semidiameter(z)
    end
    return SVector{2,T}(ρ, ϕ)
end

function onsurface(surf::AsphericSurface{T,3,Q}, p::SVector{3,T}) where {T<:Real,Q}
    ρ, ϕ = uv(surf, p)
    if ρ > one(T)
        return false
    else
        surfpoint = point(surf, ρ, ϕ)
        return samepoint(p[3], surfpoint[3])
    end
end

function inside(surf::AsphericSurface{T,3,Q}, p::SVector{3,T}) where {T<:Real,Q}
    ρ, ϕ = uv(surf, p)
    if ρ > one(T)
        return false
    else
        surfpoint = point(surf, ρ, ϕ)
        return p[3] < surfpoint[3]
    end
end

#########################################################################################################

# Assumes the ray has been transformed into the canonical coordinate frame which has the vertical axis passing through (0,0,0) and aligned with the z axis.
function surfaceintersection(surf::AcceleratedParametricSurface{T,3,AsphericSurface{T,3,Q}}, r::AbstractRay{T,3}) where {T<:Real,Q}
    cylint = surfaceintersection(surf.surface.boundingcylinder, r)
    if cylint isa EmptyInterval{T}
        return EmptyInterval(T)
    else
        if doesintersect(surf.triangles_bbox, r) || inside(surf.triangles_bbox, origin(r))
            surfint = triangulatedintersection(surf, r)
            if !(surfint isa EmptyInterval{T})
                return intervalintersection(cylint, surfint)
            end
        end
        # hasn't hit the surface
        if lower(cylint) isa RayOrigin{T} && upper(cylint) isa Infinity{T}
            if inside(surf.surface, origin(r))
                return Interval(RayOrigin(T), Infinity(T))
            else
                return EmptyInterval(T)
            end
            # otherwise check that the intersection is underneath the surface
        else
            p = point(closestintersection(cylint, false))
            ρ, ϕ = uv(surf, p)
            surfpoint = point(surf.surface, ρ, ϕ)
            if p[3] < surfpoint[3]
                return cylint # TODO!! UV (and interface) issues?
            else
                return EmptyInterval(T)
            end
        end
    end
end

function AcceleratedParametricSurface(surf::T, numsamples::Int = 17; interface::NullOrFresnel{S} = NullInterface(S)) where {S<:Real,N,T<:AsphericSurface{S,N}}
    # Zernike uses ρ, ϕ uv space so need to modify extension of triangulation
    a = AcceleratedParametricSurface(surf, triangulate(surf, numsamples, true, false, true, false), interface = interface)
    emptytrianglepool!(S)
    return a
end

function BoundingBox(surf::AsphericSurface{T,3,Q}) where {T<:Real,Q}
    xmin = -semidiameter(surf)
    xmax = semidiameter(surf)
    ymin = -semidiameter(surf)
    ymax = semidiameter(surf)
    # aspherics can be more comlpicated, so just take sum of all negative and all positives
    amin = Q > 0 ? sum(k < zero(T) ? k * surf.semidiameter^(m) : zero(T) for (m, k) in surf.aspherics) : zero(T)
    amax = Q > 0 ? sum(k > zero(T) ? k * surf.semidiameter^(m) : zero(T) for (m, k) in surf.aspherics) : zero(T)
    zmin = amin
    zmax = amax
    # curvature only goes one way
    q = one(T) - (one(T) + surf.conic) * surf.curvature^2 * surf.semidiameter^2
    if q < zero(T)
        throw(ErrorException("The surface is invalid, no bounding box can be constructed"))
    end
    hmax = surf.curvature * surf.semidiameter^2 / (one(T) + sqrt(q))
    if hmax > zero(T)
        zmax += hmax
    else
        zmin += hmax
    end
    return BoundingBox(xmin, xmax, ymin, ymax, zmin, zmax)
end
