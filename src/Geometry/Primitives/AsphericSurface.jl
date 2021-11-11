# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.



#########################################################################################################

#"pseudo-types" of aspheres
"""
AphericSurfaces polynomial evaluation is optimized for no terms `CONIC`, odd terms `ODD`, even terms `EVEN` or any `ODDEVEN`
"""
@enum AsphericSurfaceType CONIC ODD EVEN ODDEVEN


"""
    AsphericSurface{T,N,Q,M} <: ParametricSurface{T,N}

Surface incorporating an aspheric polynomial - radius, conic and aspherics are defined relative to absolute semi-diameter,.
`T` is the datatype, `N` is the dimensionality, `Q` is the number of aspheric terms, and `M` is the type of aspheric polynomial. 


```julia
AsphericSurface(semidiameter; radius, conic, aspherics=nothing, normradius = semidiameter)
```

The surface is centered at the origin and treated as being the cap of an infinite cylinder, thus creating a true half-space.
Outside of `0 <= ρ <= 1` the height of the surface is not necessarily well defined, so NaN may be returned.

`aspherics`  aspherics should be vectors containing tuples of the form (i, v) where i is the polynomial power of the aspheric term

The sag is defined by the equation

```math
z(r,\\phi) = \\frac{cr^2}{1 + \\sqrt{1 - (1+k)c^2r^2}} + \\sum_{i}^{Q}\\alpha_ir^{i}
```

where ``\\rho = \\frac{r}{\\texttt{normradius}}``, ``c = \\frac{1}{\\texttt{radius}}``, and ``k = \\texttt{conic}`` .

The function checks if the aspheric terms are missing, even, odd or both and uses an efficient polynomial evaluation strategy.

"""
struct AsphericSurface{T,N,Q,M} <: ParametricSurface{T,N} 
    semidiameter::T
    curvature::T
    conic::T 
    aspherics::SVector{Q,T}   
    normradius::T
    boundingcylinder::Cylinder{T,N}

    function AsphericSurface(M::AsphericSurfaceType, semidiameter::T, curvature::T, conic::T, aspherics::SVector{Q,T}, normradius::T, boundingcylinder) where {T<:Real,Q}
        new{T,3,Q,M}(semidiameter, curvature, conic, aspherics, normradius, boundingcylinder)
    end

    function AsphericSurface(semidiameter::T; radius::T = typemax(T), conic::T = zero(T), aspherics::Union{Nothing,Vector{Tuple{Int,T}}} = nothing, normradius::T = semidiameter) where {T<:Real}
        @assert semidiameter > 0
        @assert !isnan(semidiameter) && !isnan(radius) && !isnan(conic)
        @assert one(T) - (1 / radius)^2 * (conic + one(T)) * semidiameter^2 > 0 "Invalid surface (conic/radius combination: $radius, $conic)"

        acs = []
        if aspherics===nothing
            M = CONIC
         else
            asphericTerms = [i for (i, ) in aspherics]
            minAsphericTerm = minimum(asphericTerms)
            maxAsphericTerm = maximum(asphericTerms)
            @assert minAsphericTerm > 0 "Aspheric Terms must be Order 1 or higher ($minAsphericTerm)"
            acs = zeros(T, maxAsphericTerm )
            for (i, k) in aspherics
                acs[i] = k
            end
            odd = any([isodd(i) && a != zero(T) for (i, a) in enumerate(acs)])
            even = any([iseven(i) && a != zero(T) for (i, a) in enumerate(acs)])
            if odd && even
                M = ODDEVEN
            elseif even
                M = EVEN
                acs = [acs[2i] for i in 1:(maxAsphericTerm ÷ 2)] 
            elseif odd
                M = ODD
                acs = [acs[2i-1] for i in 1:((maxAsphericTerm+1) ÷ 2)]
            else #there are no nonzero aspherics terms in the list
                M = CONIC
                acs = []
            end
        end
        Q = length(acs)
        surf = new{T,3, Q, M}(semidiameter, 1/radius, conic, acs, normradius, Cylinder(semidiameter, interface = opaqueinterface(T)))
        return surf
    end

end

"""
    EvenAsphericSurface(semidiameter, curvature::T, conic::T, aspherics::Vector{T}; normradius::T=semidiameter)

Surface incorporating an aspheric polynomial - radius, conic and aspherics are defined relative to absolute semi-diameter.

`aspherics` should be an array of the even coefficients of the aspheric polynomial starting with A2


"""
function EvenAsphericSurface(semidiameter::T, curvature::T, conic::T, aspherics::Vector{T}; normradius::T=semidiameter) where T<:Real
    Q=length(aspherics)
    AsphericSurface(EVEN, semidiameter, curvature, conic, SVector{Q,T}(aspherics), normradius, Cylinder(semidiameter, interface = opaqueinterface(T)))
end

"""
    OddAsphericSurface(semidiameter, curvature::T, conic::T, aspherics::Vector{T}; normradius::T=semidiameter)

Surface incorporating an aspheric polynomial - radius, conic and aspherics are defined relative to absolute semi-diameter.

`aspherics`  should be an array of the odd coefficients of the aspheric polynomial starting with A1

"""
function OddAsphericSurface(semidiameter::T, curvature::T, conic::T, aspherics::Vector{T}; normradius::T=semidiameter) where T<:Real
    Q=length(aspherics)
    AsphericSurface(ODD, semidiameter, curvature, conic, SVector{Q,T}(aspherics), normradius, Cylinder(semidiameter, interface = opaqueinterface(T)))
end

"""
    OddEvenAsphericSurface(semidiameter, curvature::T, conic::T, aspherics::Vector{T}; normradius::T=semidiameter)

Surface incorporating an aspheric polynomial - radius, conic and aspherics are defined relative to absolute semi-diameter.

`aspherics` should be an array of the both odd and even coefficients of the aspheric polynomial starting with A1

"""
function OddEvenAsphericSurface(semidiameter::T, curvature::T, conic::T, aspherics::Vector{T}; normradius::T=semidiameter) where T<:Real
    Q=length(aspherics)
    AsphericSurface(ODDEVEN, semidiameter, curvature, conic, SVector{Q,T}(aspherics), normradius, Cylinder(semidiameter, interface = opaqueinterface(T)))
end

#don't have a function for CONIC as it would better to directly solve for the surface intersection instead of the rootfinding of a ParametricSurface

"""
    asphericType(surf::AsphericSurface)

Query the polynomial type of `asp.  Returns CONIC, ODD, EVEN, or ODDEVEN. CONIC corresponds to no aspheric terms, ODD 
means it only has odd aspheric terms, EVEN means only even aspheric terms and ODDEVEN means both even and odd terms. 

This function is to enable proper interpretation of `surf.aspherics` by any optimization routines that directly query the aspheric coefficients.
"""
asphericType(z::AsphericSurface{T,3,Q,M}) where {T<:Real,Q,M} = M 

export AsphericSurface, EvenAsphericSurface, OddAsphericSurface, OddEvenAsphericSurface, asphericType, EVEN, CONIC, ODD, ODDEVEN


uvrange(::AsphericSurface{T,N}) where {T<:Real,N} = ((zero(T), one(T)), (-T(π), T(π))) # ρ and ϕ

semidiameter(z::AsphericSurface{T}) where {T<:Real} = z.semidiameter
halfsizeu(z::AsphericSurface{T}) where {T<:Real} = semidiameter(z)
halfsizev(z::AsphericSurface{T}) where {T<:Real} = semidiameter(z)

boundingobj(z::AsphericSurface{T}) where {T<:Real} = z.boundingcylinder

prod_step(z::AsphericSurface{T,N,Q,ODDEVEN}, r, r2) where {T<:Real,N,Q} = r, r
prod_step(z::AsphericSurface{T,N,Q,ODD}, r, r2) where {T<:Real,N,Q} = r, r2
prod_step(z::AsphericSurface{T,N,Q,EVEN}, r, r2) where {T<:Real,N,Q} = r2, r2

function point(z::AsphericSurface{T,3,Q,M}, ρ::T, ϕ::T)::SVector{3,T} where {T<:Real,Q,M}
    rad = z.semidiameter
    r = ρ * rad
    r2 = r^2
    t = one(T) - z.curvature^2 * (z.conic + one(T)) * r^2
    if t < zero(T)
        return SVector{3,T}(NaN, NaN, NaN)
    end
    h = z.curvature * r2 / (one(T) + sqrt(t))
   # sum aspheric
    if M != CONIC
        prod, step = prod_step(z, r, r2)  #multiple dispatch on M
        asp,rest = Iterators.peel(z.aspherics)
        h += asp * prod
        for asp in rest
            prod *= step
            h += asp * prod
        end
    end
    return SVector{3,T}(r * cos(ϕ), r * sin(ϕ), h)
end

partial_prod_step(z::AsphericSurface{T,3,Q,EVEN}, r::T, r2::T) where {T<:Real,Q} = r, r2, 2:2:2Q
partial_prod_step(z::AsphericSurface{T,3,Q,ODD}, r::T, r2::T) where {T<:Real,Q} = one(T), r2, 1:2:(2Q-1)
partial_prod_step(z::AsphericSurface{T,3,Q,ODDEVEN}, r::T, r2::T) where {T<:Real,Q} = one(T), r, 1:1:Q

function partials(z::AsphericSurface{T,3,Q,M}, ρ::T, ϕ::T)::Tuple{SVector{3,T},SVector{3,T}} where {T<:Real,Q,M}
    rad = z.semidiameter
    r = ρ * rad
    r2 = r*r
    t = one(T) - z.curvature^2 * (z.conic + one(T)) * r^2
    if t < zero(T)
        return SVector{3,T}(NaN, NaN, NaN), SVector{3,T}(NaN, NaN, NaN)
    end
    dhdρ = rad * z.curvature * r * sqrt(t) / t
    # sum aspherics partial
    if M != CONIC
        prod, step, mIter = partial_prod_step(z, r, r2)
    
        ((m, asp), rest) = Iterators.peel(zip(mIter,z.aspherics))
        dhdρ += rad * asp * 2 * prod #first term m=1*2 and prod = r
        for (m,asp) in rest
            prod *= step 
            dhdρ += rad * m * asp * prod 
        end
    end
    dhdϕ = zero(T)
    cosϕ = cos(ϕ)
    sinϕ = sin(ϕ)
    return SVector{3,T}(rad * cosϕ, rad * sinϕ, dhdρ), SVector{3,T}(r * -sinϕ, r * cosϕ, dhdϕ)
end

function normal(z::AsphericSurface{T,3}, ρ::T, ϕ::T)::SVector{3,T} where {T<:Real}
    du, dv = partials(z, ρ, ϕ)
    if ρ == zero(T) && norm(dv) == zero(T)
        # in cases where there is no δϕ at ρ = 0 (i.e. anything which is rotationally symetric)
        # then we get some big problems, hardcoding this case solves the problems
        return SVector{3,T}(0, 0, 1)
    end
    return normalize(cross(du, dv))
end

function uv(z::AsphericSurface{T,3}, p::SVector{3,T}) where {T<:Real}
    # avoid divide by zero for ForwardDiff
    ϕ = NaNsafeatan(p[2], p[1])
    if p[1] == zero(T) && p[2] == zero(T)
        ρ = zero(T)
    else
        ρ = sqrt(p[1]^2 + p[2]^2) / semidiameter(z)
    end
    return SVector{2,T}(ρ, ϕ)
end

function onsurface(surf::AsphericSurface{T,3}, p::SVector{3,T}) where {T<:Real}
    ρ, ϕ = uv(surf, p)
    if ρ > one(T)
        return false
    else
        surfpoint = point(surf, ρ, ϕ)
        return samepoint(p[3], surfpoint[3])
    end
end

function inside(surf::AsphericSurface{T,3}, p::SVector{3,T}) where {T<:Real}
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
function surfaceintersection(surf::AcceleratedParametricSurface{T,3,AsphericSurface{T,3}}, r::AbstractRay{T,3}) where {T<:Real}
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

function asphericSag(surf::AsphericSurface{T,3,Q, EVEN}) where {T<:Real,Q}
    amin = sum(k < zero(T) ? k * surf.semidiameter^(2m) : zero(T) for (m, k) in enumerate(surf.aspherics)) 
    amax = sum(k > zero(T) ? k * surf.semidiameter^(2m) : zero(T) for (m, k) in enumerate(surf.aspherics)) 
    return amin, amax
end

function asphericSag(surf::AsphericSurface{T,3,Q, ODD}) where {T<:Real,Q}
    amin = sum(k < zero(T) ? k * surf.semidiameter^(2m-1) : zero(T) for (m, k) in enumerate(surf.aspherics)) 
    amax = sum(k > zero(T) ? k * surf.semidiameter^(2m-1) : zero(T) for (m, k) in enumerate(surf.aspherics)) 
    return amin, amax
end

function asphericSag(surf::AsphericSurface{T,3,Q, ODDEVEN}) where {T<:Real,Q}
    amin = sum(k < zero(T) ? k * surf.semidiameter^(m) : zero(T) for (m, k) in enumerate(surf.aspherics)) 
    amax = sum(k > zero(T) ? k * surf.semidiameter^(m) : zero(T) for (m, k) in enumerate(surf.aspherics)) 
    return amin, amax    
end

function asphericSag(surf::AsphericSurface{T,3,Q, CONIC}) where {T<:Real,Q}
    amin = zero(T)
    amax = zero(T)
    return amin, amax    
end

function BoundingBox(surf::AsphericSurface{T,3}) where {T<:Real}
    xmin = -semidiameter(surf)
    xmax = semidiameter(surf)
    ymin = -semidiameter(surf)
    ymax = semidiameter(surf)
    # aspherics can be more comlpicated, so just take sum of all negative and all positives
    amin, amax = asphericSag(surf)
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
