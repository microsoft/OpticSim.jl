abstract type EmissionShape end
abstract type EllipseEmission <: EmissionShape end
abstract type RectEmission <: EmissionShape end
abstract type PointEmission <: EmissionShape end

"""
    RayOriginGenerator{T<:Real}

Generates 3D points in world space which serve as origins for rays.
"""
abstract type RayOriginGenerator{T<:Real} end

abstract type AbstractRayGenerator{T<:Real} end

"""
    GeometricRayGenerator{T,O<:RayOriginGenerator{T}} <: AbstractRayGenerator{T}

Generates geometric [`Ray`](@ref)s according to the specific implementation of the subclass.
"""
abstract type GeometricRayGenerator{T,O<:RayOriginGenerator{T}} <: AbstractRayGenerator{T} end

"""
    OpticalRayGenerator{T} <: AbstractRayGenerator{T}

Generates [`OpticalRay`](@ref)s according to the specific implementation of the subclass.
"""
abstract type OpticalRayGenerator{T} <: AbstractRayGenerator{T} end
export AbstractRayGenerator, GeometricRayGenerator, OpticalRayGenerator

Base.iterate(a::AbstractRayGenerator, state = 1) = state > length(a) ? nothing : (generateray(a, state - 1), state + 1)
Base.getindex(a::AbstractRayGenerator, index) = generateray(a, index)
Base.firstindex(a::AbstractRayGenerator) = 0
Base.lastindex(a::AbstractRayGenerator) = length(a) - 1
Base.copy(a::AbstractRayGenerator) = a # most don't have any heap allocated stuff so don't really need copying

"""
    randinsolidangle(direction::SVector{3,T}, uvec::SVector{3,T}, vvec::SVector{3,T}, θmax::T)

Generates a unit vector pointing somewhere within the cone with half angle `θmax` around `direction`. `uvec` and `vvec` should be orthogonal to each other and `direction`.
"""
function randinsolidangle(direction::SVector{3,T}, uvec::SVector{3,T}, vvec::SVector{3,T}, θmax::T)::SVector{3,T} where {T<:Real}
    ϕ = rand(T) * 2π
    θ = NaNsafeacos(one(T) + rand(T) * (cos(θmax) - 1))
    return normalize(sin(θ) * (cos(ϕ) * uvec + sin(ϕ) * vvec) + cos(θ) * direction)
end

################################################################

function getuvvecs(direction::SVector{3,T}, rotationvec::SVector{3,T}) where {T<:Real}
    direction = normalize(direction)
    if abs(dot(rotationvec, direction)) == one(T)
        rotationvec = SVector{3,T}(1.0, 0.0, 0.0)
    end
    uvec = normalize(cross(normalize(rotationvec), direction))
    vvec = normalize(cross(direction, uvec))
    return direction, uvec, vvec
end

"""
    RandomRectOriginPoints{T} <: RayOriginGenerator{T}

Generates ray origins randomly within a rectangle.

```julia
RandomRectOriginPoints(numrays, halfsizeu, halfsizev; position = (0.0, 0.0, 0.0), direction = (0.0, 0.0, -1.0), rotationvec = (0.0, 1.0, 0.0))
```
"""
struct RandomRectOriginPoints{T} <: RayOriginGenerator{T}
    position::SVector{3,T}
    direction::SVector{3,T}
    uvec::SVector{3,T}
    vvec::SVector{3,T}
    halfsizeu::T
    halfsizev::T
    nrays::Int
    function RandomRectOriginPoints(numrays::Int, halfsizeu::T, halfsizev::T; position::SVector{3,T} = SVector{3,T}(0.0, 0.0, 0.0), direction::SVector{3,T} = SVector{3,T}(0.0, 0.0, -1.0), rotationvec::SVector{3,T} = SVector{3,T}(0.0, 1.0, 0.0)) where {T<:Real}
        new{T}(position, getuvvecs(direction, rotationvec)..., halfsizeu, halfsizev, numrays)
    end
end

"""
    GridRectOriginPoints{T} <: RayOriginGenerator{T}

Generates ray origins on a rectangular grid.

```julia
GridRectOriginPoints(numraysu, numraysv, halfsizeu, halfsizev; position = (0.0, 0.0, 0.0), direction = (0.0, 0.0, -1.0), rotationvec = (0.0, 1.0, 0.0))
```
"""
struct GridRectOriginPoints{T} <: RayOriginGenerator{T}
    position::SVector{3,T}
    direction::SVector{3,T}
    uvec::SVector{3,T}
    vvec::SVector{3,T}
    halfsizeu::T
    halfsizev::T
    nraysu::Int
    nraysv::Int
    function GridRectOriginPoints(numraysu::Int, numraysv::Int, halfsizeu::T, halfsizev::T; position::SVector{3,T} = SVector{3,T}(0.0, 0.0, 0.0), direction::SVector{3,T} = SVector{3,T}(0.0, 0.0, -1.0), rotationvec::SVector{3,T} = SVector{3,T}(0.0, 1.0, 0.0)) where {T<:Real}
        new{T}(position, getuvvecs(direction, rotationvec)..., halfsizeu, halfsizev, numraysu, numraysv)
    end
end

"""
    RandomEllipseOriginPoints{T} <: RayOriginGenerator{T}

Generates ray origins randomly across an ellipse.

```julia
RandomEllipseOriginPoints(numrays, halfsizeu, halfsizev; position = (0.0, 0.0, 0.0), direction = (0.0, 0.0, -1.0), rotationvec = (0.0, 1.0, 0.0))
```
"""
struct RandomEllipseOriginPoints{T} <: RayOriginGenerator{T}
    position::SVector{3,T}
    direction::SVector{3,T}
    uvec::SVector{3,T}
    vvec::SVector{3,T}
    halfsizeu::T
    halfsizev::T
    nrays::Int
    function RandomEllipseOriginPoints(numrays::Int, halfsizeu::T, halfsizev::T; position::SVector{3,T} = SVector{3,T}(0.0, 0.0, 0.0), direction::SVector{3,T} = SVector{3,T}(0.0, 0.0, -1.0), rotationvec::SVector{3,T} = SVector{3,T}(0.0, 1.0, 0.0)) where {T<:Real}
        new{T}(position, getuvvecs(direction, rotationvec)..., halfsizeu, halfsizev, numrays)
    end
end

"""
    HexapolarOriginPoints{T} <: RayOriginGenerator{T}

Generates ray origins in a hexapolar pattern.

```julia
HexapolarOriginPoints(nrings::Int, halfsizeu::T, halfsizev::T; position = (0.0, 0.0, 0.0), direction = (0.0, 0.0, -1.0), rotationvec = (0.0, 1.0, 0.0))
```
"""
struct HexapolarOriginPoints{T} <: RayOriginGenerator{T}
    position::SVector{3,T}
    direction::SVector{3,T}
    uvec::SVector{3,T}
    vvec::SVector{3,T}
    halfsizeu::T
    halfsizev::T
    nrings::Int
    function HexapolarOriginPoints(nrings::Int, halfsizeu::T, halfsizev::T; position::SVector{3,T} = SVector{3,T}(0.0, 0.0, 0.0), direction::SVector{3,T} = SVector{3,T}(0.0, 0.0, -1.0), rotationvec::SVector{3,T} = SVector{3,T}(0.0, 1.0, 0.0)) where {T<:Real}
        new{T}(position, getuvvecs(direction, rotationvec)..., halfsizeu, halfsizev, nrings)
    end
end

"""
    OriginPoint{T} <: RayOriginGenerator{T}

Single point origin for a source.

```julia
OriginPoint{T}(numrays; position = (0.0, 0.0, 0.0), direction = (0.0, 0.0, -1.0), rotationvec = (0.0, 1.0, 0.0))
```
"""
struct OriginPoint{T} <: RayOriginGenerator{T}
    position::SVector{3,T}
    direction::SVector{3,T}
    uvec::SVector{3,T}
    vvec::SVector{3,T}
    nrays::Int
    function OriginPoint{T}(numrays::Int; position::SVector{3,T} = SVector{3,T}(0.0, 0.0, 0.0), direction::SVector{3,T} = SVector{3,T}(0.0, 0.0, -1.0), rotationvec::SVector{3,T} = SVector{3,T}(0.0, 1.0, 0.0)) where {T<:Real}
        new{T}(position, getuvvecs(direction, rotationvec)..., numrays)
    end
end

export RandomRectOriginPoints, GridRectOriginPoints, HexapolarOriginPoints, RandomEllipseOriginPoints, OriginPoint

Base.length(a::RayOriginGenerator) = a.nrays
Base.length(a::GridRectOriginPoints) = a.nraysu * a.nraysv
Base.length(a::HexapolarOriginPoints) = 1 + round(Int, (a.nrings * (a.nrings + 1) / 2) * 6)

direction(a::RayOriginGenerator) = a.direction
position(a::RayOriginGenerator) = a.position
uvec(a::RayOriginGenerator) = a.uvec
vvec(a::RayOriginGenerator) = a.vvec

"""
    genorigin(o::RayOriginGenerator{T}, n::Int) -> SVector{3,T}

Generate origin positions for rays based on the type of the generator, e.g., randomly within a rectangle or ellipse. `n` is the index of the point being generated, starting from 0.
This has little meaning for random generators, but is important for `HexapolarOriginPoints` and `GridRectOriginPoints`.
"""
genorigin(a::OriginPoint, ::Int) = a.position
function genorigin(a::RandomRectOriginPoints{T}, ::Int) where {T<:Real}
    u = 2 * rand(T) * a.halfsizeu - a.halfsizeu
    v = 2 * rand(T) * a.halfsizev - a.halfsizev
    return a.position + u * a.uvec + v * a.vvec
end
function genorigin(a::RandomEllipseOriginPoints{T}, ::Int) where {T<:Real}
    ϕ = rand(T) * 2π
    ρ = sqrt(rand(T))
    u = cos(ϕ) * a.halfsizeu
    v = sin(ϕ) * a.halfsizev
    return a.position + ρ * (u * a.uvec + v * a.vvec)
end
function genorigin(a::HexapolarOriginPoints{T}, n::Int) where {T<:Real}
    n = mod(n, length(a))
    if n == 0
        return a.position
    else
        t = 1
        ringi = 1
        for i in 1:(a.nrings)
            t += 6 * i
            if n < t
                ringi = i
                break
            end
        end
        ρ = ringi / a.nrings
        pind = n - (t - 6 * ringi)
        ϕ = (pind / (6 * ringi)) * 2π
        u = cos(ϕ) * a.halfsizeu
        v = sin(ϕ) * a.halfsizev
        return a.position + ρ * (u * a.uvec + v * a.vvec)
    end
end
function genorigin(a::GridRectOriginPoints{T}, n::Int) where {T<:Real}
    n = mod(n, length(a))
    v = a.nraysv == 1 ? zero(T) : 2 * Int(floor(n / a.nraysu)) / (a.nraysv - 1) - 1.0
    u = a.nraysu == 1 ? zero(T) : 2 * mod(n, a.nraysu) / (a.nraysu - 1) - 1.0
    return a.position + (a.halfsizeu * u * a.uvec) + (a.halfsizev * v * a.vvec)
end

################################################################

"""
    CollimatedSource{T,O} <: GeometricRayGenerator{T,O}

Source which generates collimated rays.

```julia
CollimatedSource(generator::RayOriginGenerator)
```
"""
struct CollimatedSource{T,O<:RayOriginGenerator{T}} <: GeometricRayGenerator{T,O}
    generator::O
    function CollimatedSource(generator::O) where {T<:Real,O<:RayOriginGenerator{T}}
        new{T,O}(generator)
    end
end

"""
    GridSource{T,O} <: GeometricRayGenerator{T,O}

Source which generates rays in directions which fall on an even grid on a rectangle subtended by angles halfangleu, halfanglev.

```julia
GridSource(generator::RayOriginGenerator, numraysu, numraysv, halfangleu, halfanglev)
```
"""
struct GridSource{T,O<:RayOriginGenerator{T}} <: GeometricRayGenerator{T,O}
    generator::O
    nraysu::Int
    nraysv::Int
    halfangleu::T
    halfanglev::T
    function GridSource(generator::O, numraysu::Int, numraysv::Int, halfangleu::T, halfanglev::T) where {T<:Real,O<:RayOriginGenerator{T}}
        new{T,O}(generator, numraysu, numraysv, halfangleu, halfanglev)
    end
end

"""
    RandomSource{T,O} <: GeometricRayGenerator{T,O}

Source which generates rays in directions sampled randomly from the solid angle ±θmax centred on the direction of the source.

```julia
RandomSource(generator::RayOriginGenerator, numrays = 1, θmax = T(π / 2))
```
"""
struct RandomSource{T,O<:RayOriginGenerator{T}} <: GeometricRayGenerator{T,O}
    generator::O
    nrays::Int
    θmax::T
    function RandomSource(generator::O, numrays::Int = 1, θmax::T = T(π / 2)) where {T<:Real,O<:RayOriginGenerator{T}}
        @assert θmax <= π
        new{T,O}(generator, numrays, θmax)
    end
end

export CollimatedSource, GridSource, RandomSource

Base.length(a::RandomSource) = a.nrays * length(a.generator)
Base.length(a::GridSource) = a.nraysu * a.nraysv * length(a.generator)
Base.length(a::GeometricRayGenerator) = length(a.generator)

direction(a::GeometricRayGenerator) = direction(a.generator)
position(a::GeometricRayGenerator) = position(a.generator)
uvec(a::GeometricRayGenerator) = uvec(a.generator)
vvec(a::GeometricRayGenerator) = vvec(a.generator)


"""
    gendirection(o::GeometricRayGenerator{T}, n::Int) -> SVector{3,T}

Generate directions for rays based on the type of the generator, e.g., randomly within a cone or collimated. `n` is the index of the point being generated, starting from 0.
This has little meaning for random generators, but is important for `GridSource`.
"""
gendirection(a::RandomSource, ::Int) = randinsolidangle(a.generator.direction, a.generator.uvec, a.generator.vvec, a.θmax)
gendirection(a::CollimatedSource, ::Int) = a.generator.direction
function gendirection(a::GridSource{T}, n::Int) where {T<:Real}
    # distributing evenly across the area of the rectangle which subtends the given angle (*not* evenly across the angle range)
    dindex = mod(n, a.nraysu * a.nraysv)
    v = a.nraysv == 1 ? zero(T) : 2 * Int(floor(dindex / a.nraysu)) / (a.nraysv - 1) - 1.0
    u = a.nraysu == 1 ? zero(T) : 2 * mod(dindex, a.nraysu) / (a.nraysu - 1) - 1.0
    θu = atan(u * tan(a.halfangleu) / 2) * a.halfangleu
    θv = atan(v * tan(a.halfanglev) / 2) * a.halfanglev
    dir = cos(θv) * (cos(θu) * direction(a) + sin(θu) * uvec(a)) + sin(θv) * vvec(a)
    return dir
end

"""
    generateray(o::GeometricRayGenerator{T}, n::Int) -> Ray{T,3}

Generate geometric rays distributed according to the type of the generator. `n` is the index of the point being generated, starting from 0.
This has little meaning for random generators, but is important for `GridSource`, for example.
"""
function generateray(a::GeometricRayGenerator{T,O}, n::Int) where {T<:Real,O}
    origin = genorigin(a.generator, n)::SVector{3,T}
    direction = gendirection(a, n)::SVector{3,T}
    return Ray(origin, direction)
end

origingen(a::GeometricRayGenerator) = a.generator

################################################################

"""
    UniformOpticalSource{T,O,P} <: OpticalRayGenerator{T}

Source of OpticalRays with uniform power.

```julia
UniformOpticalSource(generator::GeometricRayGenerator, centralwavelength, power = 1.0; sourcenum = 0)
```
"""
struct UniformOpticalSource{T,O,P<:GeometricRayGenerator{T,O}} <: OpticalRayGenerator{T}
    generator::P
    power::T
    centralwavelength::T
    sourcenum::Int

    function UniformOpticalSource(generator::P, centralwavelength::Unitful.Length, power::T = one(T); sourcenum::Int = 0) where {T<:Real,O<:RayOriginGenerator{T},P<:GeometricRayGenerator{T,O}}
        λ = Unitful.ustrip(Unitful.u"μm", centralwavelength, sourcenum)
        new{T,O,P}(generator, power, λ)
    end

    function UniformOpticalSource(generator::P, centralwavelength::T, power::T = one(T); sourcenum::Int = 0) where {T<:Real,O<:RayOriginGenerator{T},P<:GeometricRayGenerator{T,O}}
        new{T,O,P}(generator, power, centralwavelength, sourcenum)
    end
end

"""
    CosineOpticalSource{T,O,P} <: OpticalRayGenerator{T}

Source of OpticalRays with power defined by:
``I(\\theta) \\approx I_0(\\cos\\theta)^C``
Where ``\\theta`` is the angle of the ray to the central direction of the source, and ``C`` is the `cosineexp` parameter.

```julia
CosineOpticalSource(generator::GeometricRayGenerator, cosineexp, centralwavelength, power = 1.0; sourcenum = 0)
```
"""
struct CosineOpticalSource{T,O,P<:GeometricRayGenerator{T,O}} <: OpticalRayGenerator{T}
    generator::P
    power::T
    centralwavelength::T
    cosineexp::T
    sourcenum::Int

    function CosineOpticalSource(generator::P, cosineexp::T, centralwavelength::Unitful.Length, power::T = one(T); sourcenum::Int = 0) where {T<:Real,O<:RayOriginGenerator{T},P<:GeometricRayGenerator{T,O}}
        λ = Unitful.ustrip(Unitful.u"μm", centralwavelength)
        new{T,O,P}(generator, power, λ, cosineexp, sourcenum)
    end

    function CosineOpticalSource(generator::P, cosineexp::T, centralwavelength::T, power::T = one(T); sourcenum::Int = 0) where {T<:Real,O<:RayOriginGenerator{T},P<:GeometricRayGenerator{T,O}}
        new{T,O,P}(generator, power, centralwavelength, cosineexp, sourcenum)
    end
end

"""
    GaussianOpticalSource{T,O,P} <: OpticalRayGenerator{T}

Source of OpticalRays with power defined by:
``I(\\theta) \\approx I_0e^{-(G_ul^2 + G_vm^2)}``
Where ``l`` and ``m`` are the direction cosines in the u and v directions to the central direction of the source, and ``G_u`` and ``G_v`` and the `gaussianu` and `gaussianv` parameters.

```julia
GaussianOpticalSource(generator::GeometricRayGenerator, gaussianu, gaussianv, centralwavelength, power = 1.0; sourcenum = 0)
```
"""
struct GaussianOpticalSource{T,P<:GeometricRayGenerator{T}} <: OpticalRayGenerator{T}
    generator::P
    power::T
    centralwavelength::T
    gaussianu::T
    gaussianv::T
    sourcenum::Int

    function GaussianOpticalSource(generator::P, gaussianu::T, gaussianv::T, centralwavelength::Unitful.Length, power::T = one(T); sourcenum::Int = 0) where {T<:Real,P<:GeometricRayGenerator{T}}
        λ = Unitful.ustrip(Unitful.u"μm", centralwavelength)
        new{T,P}(generator, power, λ, gaussianu, gaussianv, sourcenum)
    end

    function GaussianOpticalSource(generator::P, gaussianu::T, gaussianv::T, centralwavelength::T, power::T = one(T); sourcenum::Int = 0) where {T<:Real,P<:GeometricRayGenerator{T}}
        new{T,P}(generator, power, centralwavelength, gaussianu, gaussianv, sourcenum)
    end
end

export UniformOpticalSource, CosineOpticalSource, GaussianOpticalSource

function generateray(a::GaussianOpticalSource{T,P}, n::Int) where {T<:Real,P}
    r = generateray(a.generator, n)::Ray{T,3}
    l = dot(direction(r), uvec(a.generator))
    m = dot(direction(r), vvec(a.generator))
    power = a.power * exp(-(a.gaussianu * l^2 + a.gaussianv * m^2))
    return OpticalRay(r, power, a.centralwavelength, sourcenum = a.sourcenum)
end

function generateray(a::CosineOpticalSource{T,P}, n::Int) where {T<:Real,P}
    r = generateray(a.generator, n)::Ray{T,3}
    cosanglebetween = dot(direction(r), direction(a.generator))
    power = a.power * cosanglebetween^a.cosineexp
    return OpticalRay(r, power, a.centralwavelength, sourcenum = a.sourcenum)
end

"""
    generateray(o::OpticalRayGenerator{T}, n::Int) -> OpticalRay{T,3}

Generate optical rays distributed according to the type of the generator. `n` is the index of the point being generated, starting from 0.
This has little meaning for random generators, but is important for generators using `GridSource` or `GridRectOriginPoints`, for example.
"""
function generateray(a::UniformOpticalSource{T,P}, n::Int) where {T<:Real,P}
    r = generateray(a.generator, n)::Ray{T,3}
    return OpticalRay(r, one(T), a.centralwavelength, sourcenum = a.sourcenum)
end

Base.length(a::OpticalRayGenerator) = length(a.generator)
origingen(a::OpticalRayGenerator) = origingen(a.generator)
direction(a::OpticalRayGenerator) = direction(a.generator)
position(a::OpticalRayGenerator) = position(a.generator)
uvec(a::OpticalRayGenerator) = uvec(a.generator)
vvec(a::OpticalRayGenerator) = vvec(a.generator)

################################################################################################################################

"""
    PixelSource{T,C} <: RayGenerator{T}

Ray generator which encapsulates a number of subpixels, rays are generated in each subpixel in an interleaved manner.
The subpixels must be positioned correctly relative to each other when input to the constructor.

All subpixels must be coplanar, the orientation of the pixel and any display made using the pixel is taken from the first subpixel.
All subpixels must have the same number of rays (if we make subpixels be handled sequentially we wouldn't need this).
All subpixels must be of the same type.

```julia
PixelSource(subpixels::Vector{OpticalRayGenerator}; colormap = nothing, position = (0, 0, 0))
PixelSource(subpixels::Vararg{OpticalRayGenerator}; colormap = nothing, position = (0, 0, 0))
```
"""
struct PixelSource{T,P<:OpticalRayGenerator{T},C} <: OpticalRayGenerator{T}
    subpixels::Vector{P}
    colormap::SVector{C,Int}
    position::SVector{3,T}

    function PixelSource(subpixels::Vararg{P}; colormap = nothing, position::SVector{3,T} = SVector{3,T}(0, 0, 0)) where {T<:Real,P<:OpticalRayGenerator{T}}
        PixelSource(collect(subpixels), colormap = colormap, position = position)
    end

    function PixelSource(subpixels::Vector{P}; colormap = nothing, position::SVector{3,T} = SVector{3,T}(0, 0, 0)) where {T<:Real,P<:OpticalRayGenerator{T}}
        C = length(subpixels)
        # we need to map from subpixel index to color in image (e.g. for RGB pentile where there are multiple subpixels for one color)
        if (colormap === nothing)
            colormap = SVector{C,Int}(1:C)
        end
        for s in subpixels[1:end]
            @assert uvec(s) == uvec(subpixels[1]) && vvec(s) == vvec(subpixels[1])
            @assert length(s) == length(subpixels[1])
        end
        new{T,P,C}(subpixels, colormap, position)
    end
end
export PixelSource

"""
    generateray(a::PixelSource{T}, n::Int) -> OpticalRay{T,3}

Generates optical rays from all subpixels in the pixel. One ray is generated from each subpixel sequentially before looping back to the start.
"""
function generateray(a::PixelSource{T,P,C}, n::Int) where {T<:Real,P,C}
    subpixelindex = mod(n, C) + 1
    subpixeliter = Int(floor(n / C))
    ray = generateray(a.subpixels[subpixelindex], subpixeliter)::OpticalRay{T,3}
    return OpticalRay(origin(ray) + a.position, direction(ray), power(ray), wavelength(ray), sourcenum = sourcenum(ray))
end

uvec(a::PixelSource) = uvec(a.subpixels[1])
vvec(a::PixelSource) = vvec(a.subpixels[1])
colorindex(a::PixelSource{T,P,C}, i::Int) where {T<:Real,P,C} = a.colormap[mod(i, C) + 1]
Base.length(a::PixelSource) = sum(length.(a.subpixels))


"""
    OpticalSourceArray{T} <: RayGenerator{T}

Generates rays from an array of the given source at the specified locations.

```julia
OpticalSourceArray(generator::OpticalRayGenerator, positions::Vector{SVector{3,T}})
```
"""
struct OpticalSourceArray{T,S<:OpticalRayGenerator{T}} <: OpticalRayGenerator{T}
    generator::S
    positions::Vector{SVector{3,T}}

    function OpticalSourceArray(generator::S, positions::Vector{SVector{3,T}}) where {T<:Real,S<:OpticalRayGenerator{T}}
        new{T,S}(generator, positions)
    end
end
export OpticalSourceArray

Base.copy(a::OpticalSourceArray) = OpticalSourceArray(copy(a.generator), copy(a.positions))
Base.length(a::OpticalSourceArray) = length(a.generator) * length(a.positions)
Base.show(io::IO, a::OpticalSourceArray) = print(io, "OpticalSourceArray($(a.generator), $(length(a.positions)) positions)")

"""
    generateray(a::OpticalSourceArray{T}, n::Int) -> OpticalRay{T,3}

Generates optical rays from all generators in the array. One ray is generated from each element sequentially before looping back to the start of the array.
"""
function generateray(a::OpticalSourceArray{T,S}, n::Int) where {T<:Real,S}
    position = mod(n, length(a.positions)) + 1
    positioniter = Int(floor(n / length(a.positions)))
    ray = generateray(a.generator::S, positioniter)::OpticalRay{T,3}
    return OpticalRay(origin(ray) + a.positions[position], direction(ray), power(ray), wavelength(ray), sourcenum = sourcenum(ray))
end


"""
    BasicDisplayPanel{T,C} <: RayGenerator{T}

Ray generator representing a simple panel display. The panel is flat and pixles are on a regular rectangular grid. Each pixel corresponds to one image pixel.

```julia
BasicDisplayPanel(pixel::PixelSource, pitchx, pitchy, image, position = (0, 0, 0))
BasicDisplayPanel(generator::OpticalSourceArray{T,PixelSource}, pixelvals::Vector{Union{T,Vector{T}}})
```
"""
struct BasicDisplayPanel{T,P<:PixelSource{T}} <: OpticalRayGenerator{T}
    generator::OpticalSourceArray{T,P}
    pixelvals::Vector{Union{T,Vector{T}}}

    function BasicDisplayPanel(pixel::P, pitchx::T, pitchy::T, image::Array, position::SVector{3,T} = SVector{3,T}(0, 0, 0), gamma::T = 2.2) where {T<:Real,P<:PixelSource{T}}
        @assert pitchx > 0 && pitchy > 0
        h, w = size(image)
        positions = Vector{SVector{3,T}}(undef, 0)
        pixelvals = []
        for colindex in 1:w
            for rowindex in 1:h
                if size(channelview(image))[1] == h
                    k = T(channelview(image)[rowindex, colindex])^gamma
                else
                    k = T.(channelview(image)[:, rowindex, colindex]) .^ gamma
                end
                if sum(k) > 0.001
                    o = position + (colindex - (w + 1) / 2) * pitchx * uvec(pixel) + (rowindex - (h + 1) / 2) * pitchy * vvec(pixel)
                    push!(positions, o)
                    push!(pixelvals, k)
                end
            end
        end
        new{T,P}(OpticalSourceArray(pixel, positions), pixelvals)
    end

    function BasicDisplayPanel(generator::OpticalSourceArray{T,P}, pixelvals::Vector{Union{T,Vector{T}}}) where {T<:Real,P<:PixelSource{T}}
        new{T,P}(generator, pixelvals)
    end
end

export BasicDisplayPanel

Base.copy(a::BasicDisplayPanel) = BasicDisplayPanel(copy(a.generator), copy(a.pixelvals))
Base.length(a::BasicDisplayPanel) = length(a.generator)
Base.show(io::IO, a::BasicDisplayPanel) = print(io, "BasicDisplayPanel($(a.generator))")

"""
    generateray(a::BasicDisplayPanel{T}, n::Int) -> OpticalRay{T,3}

Generates optical rays from all pixels in the display. One ray is generated from each pixel sequentially before looping back to the start of the display.
"""
function generateray(a::BasicDisplayPanel{T,P}, n::Int) where {T<:Real,P}
    r = generateray(a.generator, n)::OpticalRay{T,3}
    position = mod(n, length(a.pixelvals)) + 1
    positioniter = Int(floor(n / length(a.pixelvals)))
    k = a.pixelvals[position]
    if length(k) == 1
        pow = power(r) * k
    else
        pow = power(r) * k[colorindex(a.generator.generator::P, positioniter)]
    end
    return OpticalRay(ray(r), pow, wavelength(r), sourcenum = sourcenum(r))
end

uvec(a::BasicDisplayPanel) = uvec(a.generator)
vvec(a::BasicDisplayPanel) = vvec(a.generator)

"""
    OpticalSourceGroup{T} <: RayGenerator{T}

Wrapper to group a number of separate optical sources into a single iterable object.
If the sources of the same type then performance will be improved.

```julia
OpticalSourceGroup(generators::Vector{<:OpticalRayGenerator})
OpticalSourceGroup(generators::Vararg{OpticalRayGenerator})
```
"""
abstract type OpticalSourceGroup{T} <: OpticalRayGenerator{T} end

struct OpticalSourceGroupG{T} <: OpticalSourceGroup{T}
    generators::Vector{OpticalRayGenerator{T}}
    rayspergen::Vector{Int}
end

struct OpticalSourceGroupS{T,S<:OpticalRayGenerator{T}} <: OpticalSourceGroup{T}
    generators::Vector{S}
    rayspergen::Vector{Int}
end

Base.show(io::IO, ::Type{<:OpticalSourceGroup}) = print(io, "OpticalSourceGroup")
Base.show(io::IO, a::OpticalSourceGroupS{T}) where {T<:Real} = print(io, "OpticalSourceGroupS{$T}($(length(a.generators)), $(a.rayspergen))")
Base.show(io::IO, a::OpticalSourceGroupG{T}) where {T<:Real} = print(io, "OpticalSourceGroupG{$T}($(length(a.generators)), $(a.rayspergen))")

OpticalSourceGroup(generators::Vararg{OpticalRayGenerator{T}}) where {T<:Real} = OpticalSourceGroup(collect(generators))
function OpticalSourceGroup(generators::Vector{<:OpticalRayGenerator{T}}) where {T<:Real}
    if length(generators) == 1
        return generators[1]
    end
    S = typeof(generators[1])
    for a in generators[2:end]
        if typeof(a) != S
            return OpticalSourceGroupG{T}(generators, length.(generators))
        end
    end
    OpticalSourceGroupS{T,S}(generators, length.(generators))
end
export OpticalSourceGroup

Base.copy(a::OpticalSourceGroup) = OpticalSourceGroup(a.generators)
Base.length(a::OpticalSourceGroup) = sum(a.rayspergen)

"""
    generateray(a::OpticalSourceGroup{T}, n::Int) -> OpticalRay{T,3}

Generate optical rays for each source in the group. All rays are generated for the first source, then all for the second source and so on as `n` increases.
"""
function generateray(a::OpticalSourceGroup{T}, n::Int) where {T<:Real}
    t = 0
    for (i, k) in enumerate(a.rayspergen)
        if n < t + k
            return generateray(a.generators[i], n - t)::OpticalRay{T,3}
            break
        end
        t += k
    end
    return nothing
end

"""
    RayListSource{T} <: OpticalRayGenerator{T}

Ray generator constructed manually from a list of rays which just outputs those rays in order.

```julia
RayListSource(rays::Vararg{OpticalRay})
RayListSource(rays::Vector{OpticalRay})
```
"""
struct RayListSource{T} <: OpticalRayGenerator{T}
    rays::Vector{OpticalRay{T,3}}
    RayListSource(rays::Vararg{OpticalRay{T,3}}) where {T<:Real} = new{T}(collect(rays))
    RayListSource(rays::Vector{OpticalRay{T,3}}) where {T<:Real} = new{T}(rays)
end
export RayListSource

Base.length(a::RayListSource) = length(a.rays)
origingen(::RayListSource) = nothing

function generateray(a::RayListSource, n::Int)
    return a.rays[n + 1]
end
