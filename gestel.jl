using OpticSim
using OpticSim.GlassCat
using OpticSim.Geometry

using StaticArrays
using UUIDs
using Base.Iterators

R = 125. # radius of curvature of inner spherical surface
r = 3. # radial thickness of spherical surface
θ1 = deg2rad(-49) # left limit of lens
θ2 = deg2rad(-15) # right limit of lens
ϕ1 = deg2rad(-10) # bottom limit of lens
ϕ2 = deg2rad(10) # top limit of lens
material = SCHOTT.N_BK7

δθ = deg2rad(1)
δϕ = deg2rad(1.5)

# create the spherical base plate
function sphericalsurface(R::T, r::T, θ1::T, θ2::T, ϕ1::T, ϕ2::T, material::GlassCat.AbstractGlass) where {T<:Real}
    i1 = FresnelInterface{T}(material, Air; interfacemode = Reflect) # air outside
    i2 = FresnelInterface{T}(Air, material; interfacemode = Transmit) # air inside

    s_inner = Sphere(R; interface=i2)
    s_outer = Sphere(R + r; interface=i1)

    origin = zero(SVector{3,T})
    p_left = Plane(SVector{3,T}(sin(θ1 - π/2), 0, cos(θ1 - π/2)), origin; interface=i1)
    p_right = Plane(SVector{3,T}(sin(θ2 + π/2), 0, cos(θ2 + π/2)), origin; interface=i1)
    p_bottom = Plane(SVector{3,T}(0, sin(ϕ1 - π/2), cos(ϕ1 - π/2)), origin; interface=i1)
    p_top = Plane(SVector{3,T}(0, sin(ϕ2 + π/2), cos(ϕ2 + π/2)), origin; interface=i1)

    return Object(((s_outer - s_inner) ∩ p_left ∩ p_right ∩ p_bottom ∩ p_top)())
end

# this is essentially an AxisymmetricOpticalSystem
function lensstack(radius::T, height::T) where {T<:Real}
    # aspherics = [
    #     (4, -4.61356730446999984E-04),
    #     (6, -4.39168278590500033E-02),
    #     (8, 5.87764170412300030E-02),
    #     (10, -3.60351635835699999E-02),
    #     (12, 8.41976189388700044E-03),
    # ]
    interface = FresnelInterface{T}(SCHOTT.N_BK7, Air)

    topsurface = AcceleratedParametricSurface(ZernikeSurface(radius; radius=T(-4.11), conic=T(-11.846)); interface)
    barrel = Cylinder(radius, T(2*height); interface)
    botsurface = Plane(
        SVector{3,T}(0.0, 0.0, -1.0),
        SVector{3,T}(0.0, 0.0, -height);
        vishalfsizeu=radius,
        vishalfsizev=radius,
        interface
    )

    return topsurface ∩ barrel ∩ botsurface
end

# create a bunch of lens stacks, transform them into a spherical hex grid, and group them by lens type (HEX7)
function tiledlensstacks(
    nx, ny, hexradius::T, δθ::T, δϕ::T, R::T, r::T, θ1::T, θ2::T, ϕ1::T, ϕ2::T, material::GlassCat.AbstractGlass
) where {T<:Real}
    ρ = R + 1.5r
    origin = [θ1 + θ2; ϕ1 + ϕ2] / 2
    coords = Repeat.hexcellsinbox(nx, ny)

    lensstacks::Vector{Vector{Object{<:CSGTree{T}}}} = fill([], 7)

    for (i, j) in eachcol(coords)
        lenstype = mod(i - 2j, 7) + 1

        θ, ϕ = origin + Repeat.HexBasis1()[i, j] .* [δθ, δϕ]
        transform = rotation(zero(T), θ, ϕ) * translation(zero(T), zero(T), ρ)
        transformedlensstack = lensstack(hexradius, r)(transform)

        push!(lensstacks[lenstype], Object(transformedlensstack))
    end

    return lensstacks
end

lensstacks = tiledlensstacks(5, 4, 2., δθ, δϕ, R, r, θ1, θ2, ϕ1, ϕ2, material)
baseplate = sphericalsurface(R, r, θ1, θ2, ϕ1, ϕ2, material)

colors = [color(c) for c in split("white red orange green blue cyan purple")]
properties = Properties(
    [object.id => Dict("color" => colors[i]) for i in 1:7 for object in lensstacks[i]]...,
    baseplate.id => Dict("color" => colorant"black")
)

assembly = LensAssembly(reduce(vcat, lensstacks)..., baseplate)
Vis.draw(assembly; properties)

# Vis.draw(sphericalsurface(R, r, θ1, θ2, ϕ1, ϕ2, material))
# Vis.draw(Object(lensstack(5., 2.))(Transform{Float64}(0.0, Float64(π/4), 0.0, 0.0, 0.0, -5.0)))
