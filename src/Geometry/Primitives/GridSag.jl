# MIT License

# Copyright (c) Microsoft Corporation.

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE

"""
Either `GridSagLinear` or `GridSagBicubic` - determines the interpolation between sample points in the grid for a [`GridSagSurface`](@ref).
"""
@enum GridSagInterpolation GridSagLinear GridSagBicubic
export GridSagInterpolation, GridSagLinear, GridSagBicubic

"""
    GridSagSurface{T,N,S<:Union{ZernikeSurface{T,N},ChebyshevSurface{T,N}},Nu,Nv} <: ParametricSurface{T,N}

Either a Zernike (circular) or Chebyshev (rectangular) surface with grid sag height added to the base sag.
The surface shape is determined by either a linear or a bicubic spline interpolation of the `Nu×Nv` grid of sag values, set by the `interpolation` argument taking either `GridSagLinear` or `GridSagBicubic`.

Each entry in the grid is a vector of the form ``[z, \\frac{\\partial z}{\\partial x}, \\frac{\\partial z}{\\partial y}, \\frac{\\partial^2 z}{\\partial x \\partial y}]``.
The first data item corresponds to the lower left corner of the surface, that is, the corner defined by the -u and -v limit.
Each point that follows is read across the face of the surface from left to right moving upwards.
If zero is given for the partials (and using bicubic interpolation) then the partials will be approximated using finite differences.

The sag grid can be decentered from the surface in uv space, if so the surface may become wild outside of the area over which the grid is defined.
It is advised to clip the surface to the valid area using CSG operations in this case.

A surface can also be generated from a `.GRD` file by passing in the filename as the first and only positional argument. In this case the surface will be rectangular with optional radius and conic.

See docs for [`ZernikeSurface`](@ref) and [`ChebyshevSurface`](@ref) for details of the base surface.

```julia
GridSagSurface(basesurface::Union{ZernikeSurface{T,N},ChebyshevSurface{T,N}}, sag_grid::AbstractArray{T,3}; interpolation = GridSagBicubic, decenteruv = (0, 0))
GridSagSurface{T}(filename::String; radius = Inf, conic = 0, interpolation = GridSagBicubic)
```
"""
struct GridSagSurface{T,N,S<:Union{ZernikeSurface{T,N},ChebyshevSurface{T,N}},Nu,Nv} <: ParametricSurface{T,N}
    basesurface::S
    gridsag::SVector{Nu,SVector{Nv,SVector{4,T}}}
    interpolation::GridSagInterpolation
    decenteruv::Tuple{T,T}

    function GridSagSurface(basesurface::S, sag_grid::AbstractArray{T,2}; interpolation::GridSagInterpolation = GridSagBicubic, decenteruv::Tuple{T,T} = (zero(T), zero(T))) where {T<:Real,N,S<:Union{ZernikeSurface{T,N},ChebyshevSurface{T,N}}}
        grid = []
        Nv, Nu, d = size(sag_grid)
        @assert d == 4
        for r in 1:Nu
            row = []
            for c in 1:Nv
                push!(row, SVector{4,T}(sag_grid[c, r], 0, 0, 0))
            end
            push!(grid, row)
        end
        grid = SVector{Nu,SVector{Nv,SVector{4,T}}}(grid)
        return GridSagSurface(basesurface, grid, interpolation = interpolation, decenteruv = decenteruv)
    end

    function GridSagSurface(basesurface::S, sag_grid::AbstractArray{T,3}; interpolation::GridSagInterpolation = GridSagBicubic, decenteruv::Tuple{T,T} = (zero(T), zero(T))) where {T<:Real,N,S<:Union{ZernikeSurface{T,N},ChebyshevSurface{T,N}}}
        grid = []
        Nv, Nu, d = size(sag_grid)
        @assert d == 4
        for r in 1:Nu
            row = []
            for c in 1:Nv
                push!(row, SVector{4,T}(sag_grid[c, r, :]))
            end
            push!(grid, row)
        end
        grid = SVector{Nu,SVector{Nv,SVector{4,T}}}(grid)
        return GridSagSurface(basesurface, grid, interpolation = interpolation, decenteruv = decenteruv)
    end

    function GridSagSurface(basesurface::S, sag_grid::AbstractArray{<:AbstractVector{T},2}; interpolation::GridSagInterpolation = GridSagBicubic, decenteruv::Tuple{T,T} = (zero(T), zero(T))) where {T<:Real,N,S<:Union{ZernikeSurface{T,N},ChebyshevSurface{T,N}}}
        grid = []
        Nv, Nu = size(sag_grid)
        @assert length(sag_grid[1, 1]) == 4
        for r in 1:Nu
            push!(grid, sag_grid[:, r])
        end
        grid = SVector{Nu,SVector{Nv,SVector{4,T}}}(grid)
        return GridSagSurface(basesurface, grid, interpolation = interpolation, decenteruv = decenteruv)
    end

    function GridSagSurface(basesurface::S, sag_grid::SVector{Nu,SVector{Nv,SVector{4,T}}}; interpolation::GridSagInterpolation = GridSagBicubic, decenteruv::Tuple{T,T} = (zero(T), zero(T))) where {T<:Real,N,S<:Union{ZernikeSurface{T,N},ChebyshevSurface{T,N}},Nu,Nv}
        # TODO can probably simplify this?
        if interpolation === GridSagLinear
            @warn "Linear interpolation doesn't ensure C1 continuity, may cause problems..."
        end
        Δx = 2 * halfsizeu(basesurface) / (Nu - 1)
        Δy = 2 * halfsizev(basesurface) / (Nv - 1)
        kx = 1 / 2Δx
        ky = 1 / 2Δy
        kxy = kx * ky
        dxdκ = 2 * halfsizeu(basesurface) / (Nu - 1)
        dydγ = 2 * halfsizev(basesurface) / (Nv - 1)
        grid = MVector{Nu,MVector{Nv,SVector{4,T}}}(sag_grid)
        # if no derivs are provided (and we need them, i.e. bicubic) then approximate them here
        calcderivs = interpolation === GridSagBicubic && all(all.(iszero.(getindex.(v, [2:4])) for v in sag_grid))
        for u in 1:Nu
            for v in 1:Nv
                if calcderivs
                    # use finite differences if missing derivatives
                    z = sag_grid[u][v][1]
                    # at the edges just assume the derivative continues unchanged from that between the adjacent and current points
                    # in this case we need to half Δx and/or Δy hence umod and vmod below
                    u₊₁ = u + 1
                    u₋₁ = u - 1
                    umod = 1.0
                    if u === Nu
                        u₊₁ = u
                        umod = 2
                    elseif u === 1
                        u₋₁ = 1
                        umod = 2
                    end
                    v₊₁ = v + 1
                    v₋₁ = v - 1
                    vmod = 1.0
                    if v === Nv
                        v₊₁ = v
                        vmod = 2
                    elseif v === 1
                        v₋₁ = 1
                        vmod = 2
                    end
                    dzdx = kx * (sag_grid[u₊₁][v][1] - sag_grid[u₋₁][v][1]) * umod
                    dzdy = ky * (sag_grid[u][v₊₁][1] - sag_grid[u][v₋₁][1]) * vmod
                    dzdxdy = kxy * ((sag_grid[u₊₁][v₊₁][1] - sag_grid[u₊₁][v₋₁][1]) - (sag_grid[u₋₁][v₊₁][1] - sag_grid[u₋₁][v₋₁][1])) * umod * vmod
                else
                    z, dzdx, dzdy, dzdxdy = sag_grid[u][v]
                end
                # transform the derivatives from xy to κγ (i.e [0,1] in each patch)
                dzdκ = dzdx * dxdκ
                dzdγ = dzdy * dydγ
                dzdκdγ = dzdxdy * dxdκ * dydγ
                grid[u][v] = SVector{4,T}(z, dzdκ, dzdγ, dzdκdγ)
            end
        end
        return new{T,N,S,Nu,Nv}(basesurface, SVector{Nu,SVector{Nv,SVector{4,T}}}(grid), interpolation, decenteruv)
    end

    function GridSagSurface{T}(filename::String; radius::T = typemax(T), conic::T = zero(T), interpolation::GridSagInterpolation = GridSagBicubic) where {T<:Real}
        @assert !isnan(radius) && !isnan(conic)
        @assert isfile(filename)
        if !isvalid(readuntil(filename, " "))
            f = open(filename, enc"UTF-16LE", "r")
        else
            f = open(filename, "r")
        end
        spec = readline(f)
        spec = replace(spec, "\ufeff" => "")
        specs = split(spec)
        Nv = parse(Int, specs[1])
        Nu = parse(Int, specs[2])
        Δx = parse(Float64, specs[3])
        Δy = parse(Float64, specs[4])
        units = parse(Int, specs[5])
        if units == 0
            unitmod = 1
        elseif units == 1
            unitmod = 10
        elseif units == 2
            unitmod = 25.4
        elseif units == 3
            unitmod = 1000
        else
            throw(ErrorException("Invalid units"))
        end
        halfsizeu = (Nu - 1) * Δx / 2 / unitmod
        halfsizev = (Nv - 1) * Δy / 2 / unitmod
        xdec = parse(Float64, specs[6]) / unitmod
        ydec = parse(Float64, specs[7]) / unitmod
        decenteruv = (xdec / halfsizeu, ydec / halfsizev)
        grid = Array{Float64,3}(undef, Nv, Nu, 4)
        thisrowcount = 1
        thisrownum = 1
        for l in readlines(f)
            vs = [parse(Float64, x) for x in split(l)]
            if vs[5] == 1
                z, dzdx, dzdy, d2zdxdy = zeros(4) # nodata
            else
                z, dzdx, dzdy, d2zdxdy = vs[1:4]
            end
            # flip the y direction as GRD goes from -x, +y right and down while we go from -x, -y right and up
            grid[Nv - thisrownum + 1, thisrowcount, :] = [z / unitmod, dzdx / unitmod, dzdy / unitmod, d2zdxdy / unitmod^2]
            if thisrowcount == Nu
                thisrownum += 1
                thisrowcount = 1
            else
                thisrowcount += 1
            end
        end
        close(f)
        return GridSagSurface(ChebyshevSurface(halfsizeu, halfsizev, nothing, radius = radius, conic = conic), grid, interpolation = interpolation, decenteruv = decenteruv)
    end
end
export GridSagSurface

uvrange(::Type{GridSagSurface{T,N,S,Nu,Nv}}) where {T<:Real,N,S<:Union{ZernikeSurface{T,N},ChebyshevSurface{T,N}},Nu,Nv} = uvrange(S)

@inline function gridsag(s::GridSagSurface{T,3,S,Nu,Nv}, ui::Int, vi::Int, i::Int = 0) where {T<:Real,S,Nu,Nv}
    @assert 0 <= ui < Nu
    @assert 0 <= vi < Nv
    @assert i <= 4
    if i <= 0
        return @inbounds s.gridsag[ui + 1][vi + 1]
    else
        return @inbounds s.gridsag[ui + 1][vi + 1][i]
    end
end

@inline function bicubicα(surf::GridSagSurface{T,3,S,Nu,Nv}, uli::Int, uui::Int, vli::Int, vui::Int) where {T<:Real,S,Nu,Nv}
    f00, fu00, fv00, fuv00 = gridsag(surf, uli, vli)
    f01, fu01, fv01, fuv01 = gridsag(surf, uli, vui)
    f11, fu11, fv11, fuv11 = gridsag(surf, uui, vui)
    f10, fu10, fv10, fuv10 = gridsag(surf, uui, vli)
    return SMatrix{4,4,T,16}(1, 0, -3, 2, 0, 0, 3, -2, 0, 1, -2, 1, 0, 0, -1, 1) * SMatrix{4,4,T,16}(f00, f10, fu00, fu10, f01, f11, fu01, fu11, fv00, fv10, fuv00, fuv10, fv01, fv11, fuv01, fuv11) * SMatrix{4,4,T,16}(1, 0, 0, 0, 0, 0, 1, 0, -3, 3, -2, -1, 2, -2, 1, 1)
end

@inline function bicubicζ(surf::GridSagSurface{T,3,S,Nu,Nv}, uli::Int, uui::Int, vli::Int, vui::Int, κ::T, γ::T) where {T<:Real,S,Nu,Nv}
    return (SMatrix{1,4,T,4}(1, κ, κ^2, κ^3) * bicubicα(surf, uli, uui, vli, vui) * SMatrix{4,1,T,4}(1, γ, γ^2, γ^3))[1]
end

@inline function bicubicdζ(surf::GridSagSurface{T,3,S,Nu,Nv}, uli::Int, uui::Int, vli::Int, vui::Int, κ::T, γ::T) where {T<:Real,S,Nu,Nv}
    α = bicubicα(surf, uli, uui, vli, vui)
    dκ = (SMatrix{1,4,T,4}(0, 1, 2κ, 3κ^2) * α * SMatrix{4,1,T,4}(1, γ, γ^2, γ^3))[1]
    dγ = (SMatrix{1,4,T,4}(1, κ, κ^2, κ^3) * α * SMatrix{4,1,T,4}(0, 1, 2γ, 3γ^2))[1]
    return dκ, dγ
end

@inline function normalizedcoords(Nu::Int, Nv::Int, u::T, v::T, decenter::Tuple{T,T}) where {T<:Real}
    u = u - decenter[1]
    v = v - decenter[2]
    û = (u + 1) / 2
    v̂ = (v + 1) / 2
    uli = max(min(floor(Int, û * (Nu - 1)), Nu - 2), 0)
    uui = min(uli + 1, Nu - 1)
    vli = max(min(floor(Int, v̂ * (Nv - 1)), Nv - 2), 0)
    vui = min(vli + 1, Nv - 1)
    κ = û * (Nu - 1) - uli
    γ = v̂ * (Nv - 1) - vli
    return uli, uui, vli, vui, κ, γ
end

@inline function ζ(surf::GridSagSurface{T,3,S,Nu,Nv}, u::T, v::T) where {T<:Real,S,Nu,Nv}
    uli, uui, vli, vui, κ, γ = normalizedcoords(Nu, Nv, u, v, surf.decenteruv)
    if surf.interpolation === GridSagLinear
        sagll = gridsag(surf, uli, vli, 1)
        saglu = gridsag(surf, uli, vui, 1)
        sagul = gridsag(surf, uui, vli, 1)
        saguu = gridsag(surf, uui, vui, 1)
        z = κ * (γ * saguu + (1 - γ) * sagul) + (1 - κ) * (γ * saglu + (1 - γ) * sagll)
    elseif surf.interpolation === GridSagBicubic
        z = bicubicζ(surf, uli, uui, vli, vui, κ, γ)
    end
    return z
end

@inline function dζ(surf::GridSagSurface{T,3,S,Nu,Nv}, u::T, v::T) where {T<:Real,S,Nu,Nv}
    uli, uui, vli, vui, κ, γ = normalizedcoords(Nu, Nv, u, v, surf.decenteruv)
    # calculate sag derivative at point
    if surf.interpolation === GridSagLinear
        sagll = gridsag(surf, uli, vli, 1)
        saglu = gridsag(surf, uli, vui, 1)
        sagul = gridsag(surf, uui, vli, 1)
        saguu = gridsag(surf, uui, vui, 1)
        dζdκ = (γ * saguu + (1 - γ) * sagul) - (γ * saglu + (1 - γ) * sagll)
        dζdγ = (κ * saguu + (1 - κ) * saglu) - (κ * sagul + (1 - κ) * sagll)
        # NOT C1 CONTINUOUS!!!
    elseif surf.interpolation === GridSagBicubic
        dζdκ, dζdγ = bicubicdζ(surf, uli, uui, vli, vui, κ, γ)
        # C1 continuous, but not C2 continuous (at grid boundaries)
    end
    # change variables
    dζdu = dζdκ * (Nu - 1) / 2
    dζdv = dζdγ * (Nv - 1) / 2
    return dζdu, dζdv
end

function point(surf::GridSagSurface{T,3,S}, ρ::T, ϕ::T)::SVector{3,T} where {T<:Real,S<:ZernikeSurface{T,3}}
    x, y, z = point(surf.basesurface, ρ, ϕ)
    u = ρ * cos(ϕ)
    v = ρ * sin(ϕ)
    return SVector{3,T}(x, y, z + ζ(surf, u, v))
end

function point(surf::GridSagSurface{T,3,S}, u::T, v::T)::SVector{3,T} where {T<:Real,S<:ChebyshevSurface{T,3}}
    x, y, z = point(surf.basesurface, u, v)
    return SVector{3,T}(x, y, z + ζ(surf, u, v))
end

function partials(surf::GridSagSurface{T,3,S}, ρ::T, ϕ::T)::Tuple{SVector{3,T},SVector{3,T}} where {T<:Real,S<:ZernikeSurface{T,3}}
    pρ, pϕ = partials(surf.basesurface, ρ, ϕ)
    u = ρ * cos(ϕ)
    v = ρ * sin(ϕ)
    dζdu, dζdv = dζ(surf, u, v)
    dζdρ = dζdu * cos(ϕ) + dζdv * sin(ϕ)
    dζdϕ = dζdv * u - dζdu * v
    return SVector{3,T}(pρ[1], pρ[2], pρ[3] + dζdρ), SVector{3,T}(pϕ[1], pϕ[2], pϕ[3] + dζdϕ)
end

function partials(surf::GridSagSurface{T,3,S}, u::T, v::T)::Tuple{SVector{3,T},SVector{3,T}} where {T<:Real,S<:ChebyshevSurface{T,3}}
    pu, pv = partials(surf.basesurface, u, v)
    du, dv = dζ(surf, u, v)
    return SVector{3,T}(pu[1], pu[2], pu[3] + du), SVector{3,T}(pv[1], pv[2], pv[3] + dv)
end

function normal(z::GridSagSurface{T,N}, u::T, v::T)::SVector{N,T} where {T<:Real,N}
    du, dv = partials(z, u, v)
    return normalize(cross(du, dv))
end

uv(s::GridSagSurface{T,3}, p::SVector{3,T}) where {T<:Real} = uv(s.basesurface, p)

function onsurface(surf::GridSagSurface{T,3,S}, p::SVector{3,T}) where {T<:Real,S}
    u, v = uv(surf, p)
    if abs(u) > one(T) || (abs(v) > one(T) && S <: ChebyshevSurface{T,3})
        return false
    else
        # not sure it's really necessary/correct to include the cylinder here
        surfpoint = point(surf, u, v)
        return @inbounds samepoint(p[3], surfpoint[3]) # || (ρ == 1 && p[3] < z)
    end
end

function inside(surf::GridSagSurface{T,3,S}, p::SVector{3,T}) where {T<:Real,S}
    u, v = uv(surf, p)
    if abs(u) > one(T) || (abs(v) > one(T) && S <: ChebyshevSurface{T,3})
        return false
    else
        surfpoint = point(surf, u, v)
        return @inbounds p[3] < surfpoint[3]
    end
end

#########################################################################################################

# Assumes the ray has been transformed into the canonical coordinate frame which has the vertical axis passing through (0,0,0) and aligned with the z axis.
function surfaceintersection(surf::AcceleratedParametricSurface{T,3,GridSagSurface{T,3,S,Nu,Nv}}, r::AbstractRay{T,3}) where {T<:Real,S,Nu,Nv}
    boundingintvl = surfaceintersection(boundingobj(surf.surface.basesurface), r)
    if boundingintvl isa EmptyInterval{T}
        return EmptyInterval(T)
    else
        if doesintersect(surf.triangles_bbox, r) || inside(surf.triangles_bbox, origin(r))
            surfint = triangulatedintersection(surf, r)
            if !(surfint isa EmptyInterval{T})
                return intervalintersection(boundingintvl, surfint)
            end
        end
        # hasn't hit the surface
        if lower(boundingintvl) isa RayOrigin{T} && upper(boundingintvl) isa Infinity{T}
            if inside(surf.surface, origin(r))
                return Interval(RayOrigin(T), Infinity(T))
            else
                return EmptyInterval(T)
            end
            # otherwise check that the intersection is underneath the surface
        else
            p = point(closestintersection(boundingintvl, false))
            ρ, ϕ = uv(surf, p)
            surfpoint = point(surf.surface, ρ, ϕ)
            if @inbounds p[3] < surfpoint[3]
                return boundingintvl # TODO!! UV (and interface) issues?
            else
                return EmptyInterval(T)
            end
        end
    end
end

function AcceleratedParametricSurface(surf::GridSagSurface{T,N,Z}, numsamples::Int = 17; interface::NullOrFresnel{T} = NullInterface(T)) where {T<:Real,N,Z<:ZernikeSurface{T,N}}
    # Zernike uses ρ, ϕ uv space so need to modify extension of triangulation
    a = AcceleratedParametricSurface(surf, triangulate(surf, numsamples, true, false, true, false), interface = interface)
    emptytrianglepool!(T)
    return a
end

function BoundingBox(surf::GridSagSurface{T,3}) where {T<:Real}
    bbox = BoundingBox(surf.basesurface)
    maxsag = maximum(maximum.(getindex.(surf.gridsag[i], 1) for i in 1:length(surf.gridsag)))
    minsag = minimum(minimum.(getindex.(surf.gridsag[i], 1) for i in 1:length(surf.gridsag)))
    return BoundingBox(bbox.xmin, bbox.xmax, bbox.ymin, bbox.ymax, bbox.zmin + minsag, bbox.zmax + maxsag)
end
