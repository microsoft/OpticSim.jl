# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

module Origins
export Point, RectUniform, RectGrid, Hexapolar, RectJitterGrid

using ....OpticSim
using ...Emitters
using ...Geometry
using LinearAlgebra
using Distributions
using Random

abstract type AbstractOriginDistribution{T<:Real} end

Base.iterate(a::AbstractOriginDistribution, state = 1) = state > length(a) ? nothing : (generate(a, state - 1), state + 1)
Base.getindex(a::AbstractOriginDistribution, index) = generate(a, index)
Base.firstindex(a::AbstractOriginDistribution) = 0
Base.lastindex(a::AbstractOriginDistribution) = length(a) - 1
Base.copy(a::AbstractOriginDistribution) = a # most don't have any heap allocated stuff so don't really need copying

"""
    Point{T} <: AbstractOriginDistribution{T}

Encapsulates a single point origin.

```julia
Point(position::Vec3{T}) where {T<:Real}
Point(x::T, y::T, z::T) where {T<:Real}
Point(::Type{T} = Float64) where {T<:Real}
```
"""
struct Point{T} <: AbstractOriginDistribution{T}
    origin::Vec3{T}

    function Point(position::Vec3{T}) where {T<:Real}
        return new{T}(position)
    end

    function Point(x::T, y::T, z::T) where {T<:Real}
        return new{T}(Vec3{T}(x, y, z))
    end

    function Point(::Type{T} = Float64) where {T<:Real}
        return new{T}(zero(Vec3))
    end
end

Base.length(::Point) = 1
Emitters.visual_size(::Point) = 1
Emitters.generate(o::Point, ::Int64) = o.origin

"""
    RectUniform{T} <: AbstractOriginDistribution{T}

Encapsulates a uniformly sampled rectangle with user defined number of samples.

```julia
RectUniform(width::T, height::T, samples_count::Int64) where {T<:Real}
```
"""
struct RectUniform{T} <: AbstractOriginDistribution{T}
    width::T
    height::T
    samples_count::Int64
    rng::Random.AbstractRNG

    function RectUniform(width::T, height::T, samples_count::Int64; rng=Random.GLOBAL_RNG) where {T<:Real}
        return new{T}(width, height, samples_count, rng)
    end
end

Base.length(o::RectUniform) = o.samples_count
Emitters.visual_size(o::RectUniform) = max(o.width, o.height)

# generate origin on the agrid
function Emitters.generate(o::RectUniform{T}, n::Int64) where {T<:Real}
    n = mod(n, length(o))
    u = rand(o.rng, Distributions.Uniform(-one(T), one(T)))
    v = rand(o.rng, Distributions.Uniform(-one(T), one(T)))
    return zero(Vec3{T}) + ((o.width / 2) * u * unitX3(T)) + ((o.height/2) * v * unitY3(T))
end

"""
    RectGrid{T} <: AbstractOriginDistribution{T}

Encapsulates a rectangle sampled in a grid fashion.

```julia
RectGrid(width::T, height::T, usamples::Int64, vsamples::Int64) where {T<:Real} 
```
"""
struct RectGrid{T} <: AbstractOriginDistribution{T}
    width::T
    height::T
    usamples::Int64
    vsamples::Int64
    ustep::T
    vstep::T

    function RectGrid(width::T, height::T, usamples::Int64, vsamples::Int64) where {T<:Real}
        return new{T}(width, height, usamples, vsamples, width / (usamples - 1), height / (vsamples - 1))
    end
end

Base.length(o::RectGrid) = o.usamples * o.vsamples
Emitters.visual_size(o::RectGrid) = max(o.width, o.height)

# generate origin on the grid
function Emitters.generate(o::RectGrid{T}, n::Int64) where {T<:Real}
    n = mod(n, length(o))
    v = o.vsamples == 1 ? zero(T) : 2 * Int64(floor(n / o.usamples)) / (o.vsamples - 1) - 1.0
    u = o.usamples == 1 ? zero(T) : 2 * mod(n, o.usamples) / (o.usamples - 1) - 1.0
    return zeros(Vec3{T}) + ((o.width / 2) * u * unitX3(T)) + ((o.height/2) * v * unitY3(T))
end

"""
    RectJitterGrid{T} <: AbstractOriginDistribution{T}

Encapsulates a rectangle sampled in a grid fashion with jitter.

```julia
RectGrid(width::T, height::T, ures::Int64, vres::Int64, samplesPerRegion::Int64) where {T<:Real} 
```
"""
struct RectJitterGrid{T} <: AbstractOriginDistribution{T}
    width::T
    height::T
    uResolution::Int64
    vResolution::Int64
    samplesPerRegion::Int64
    ustep::T
    vstep::T
    rng::Random.AbstractRNG

    function RectJitterGrid(width::T, height::T, ures::Int64, vres::Int64, samplesPerRegion::Int64; rng=Random.GLOBAL_RNG) where {T<:Real}
        return new{T}(width, height, ures, vres, samplesPerRegion, width / ures, height / vres, rng)
    end
end

Base.length(o::RectJitterGrid) = o.uResolution * o.vResolution * o.samplesPerRegion
Emitters.visual_size(o::RectJitterGrid) = max(o.width, o.height)

# generate origin on the grid
function Emitters.generate(o::RectJitterGrid{T}, n::Int64) where {T<:Real}
    n = mod(n, length(o))

    uu = rand(o.rng, Distributions.Uniform(zero(T), o.ustep))
    vv = rand(o.rng, Distributions.Uniform(zero(T), o.vstep))

    nn = Int64(floor(n / o.samplesPerRegion))
    v = (o.vResolution == 1) ? zero(T) : Int64(floor(nn / o.uResolution)) 
    u = (o.uResolution == 1) ? zero(T) : Int64(floor(mod(nn, o.uResolution)))

    v = v * o.vstep + vv - (o.height / 2.0)
    u = u * o.ustep + uu - (o.width / 2.0)

    return zeros(Vec3{T}) + u * unitX3(T) + v * unitY3(T)
end


"""
    Hexapolar{T} <: AbstractOriginDistribution{T}

Encapsulates an ellipse (or a circle where halfsizeu=halfsizev) sampled in an hexapolar fashion (rings).

```julia
Hexapolar(nrings::Int64, halfsizeu::T, halfsizev::T) where {T<:Real} 
```
"""
struct Hexapolar{T} <: AbstractOriginDistribution{T}
    halfsizeu::T
    halfsizev::T
    nrings::Int64

    function Hexapolar(nrings::Int64, halfsizeu::T, halfsizev::T) where {T<:Real} 
        return new{T}(halfsizeu, halfsizev, nrings)
    end
end

Base.length(o::Hexapolar) = 1 + round(Int64, (o.nrings * (o.nrings + 1) / 2) * 6)
Emitters.visual_size(o::Hexapolar) = max(o.halfsizeu*2, o.halfsizev*2)

function Emitters.generate(o::Hexapolar{T}, n::Int64) where {T<:Real}
    n = mod(n, length(o))
    if n == 0
        return zeros(Vec3{T})
    else
        t = 1
        ringi = 1
        for i in 1:(o.nrings)
            t += 6 * i
            if n < t
                ringi = i
                break
            end
        end
        ρ = ringi / o.nrings
        pind = n - (t - 6 * ringi)
        ϕ = (pind / (6 * ringi)) * 2π
        u = cos(ϕ) * o.halfsizeu
        v = sin(ϕ) * o.halfsizev
        return zeros(Vec3{T}) + ρ * (u * unitX3(T) + v * unitY3(T))
    end
end

end # module Origins
