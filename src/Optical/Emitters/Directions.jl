# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

module Directions
export Constant, RectGrid, UniformCone, HexapolarCone

using ....OpticSim
using ...Emitters
using ...Geometry
using LinearAlgebra

abstract type AbstractDirectionDistribution{T<:Real} end

Base.iterate(a::AbstractDirectionDistribution, state = 1) = state > length(a) ? nothing : (generate(a, state - 1), state + 1)
Base.getindex(a::AbstractDirectionDistribution, index) = generate(a, index)
Base.firstindex(a::AbstractDirectionDistribution) = 0
Base.lastindex(a::AbstractDirectionDistribution) = length(a) - 1
Base.copy(a::AbstractDirectionDistribution) = a # most don't have any heap allocated stuff so don't really need copying

"""
    Constant{T} <: AbstractDirectionDistribution{T}

Encapsulates a single ray direction, where the default direction is unitZ3 [0, 0, 1].

```julia
Constant(direction::Vec3{T}) where {T<:Real}
Constant(::Type{T} = Float64) where {T<:Real}
```
"""
struct Constant{T} <: AbstractDirectionDistribution{T}
    direction::Vec3{T}

    function Constant(direction::T...) where{T<:Real}
        return new{T}(direction)
    end

    function Constant(direction::Vec3{T}) where {T<:Real}
        return new{T}(direction)
    end

    function Constant(::Type{T} = Float64) where {T<:Real}
        direction = unitZ3(T)
        return new{T}(direction)
    end
end

Base.length(d::Constant) = 1
Emitters.generate(d::Constant, ::Integer) = d.direction

"""
    RectGrid{T} <: AbstractDirectionDistribution{T}

Encapsulates a single ray direction, where the default direction is unitZ3 [0, 0, 1].

```julia
Constant(direction::Vec3{T}) where {T<:Real}
Constant(::Type{T} = Float64) where {T<:Real}
```
"""
struct RectGrid{T} <: AbstractDirectionDistribution{T}
    direction::Vec3{T}
    halfangleu::T
    halfanglev::T
    nraysu::Integer
    nraysv::Integer
    uvec::Vec3{T}
    vvec::Vec3{T}

    function RectGrid(direction::Vec3{T}, halfangleu::T, halfanglev::T, numraysu::Integer, numraysv::Integer) where {T<:Real}
        (uvec, vvec) = get_uv_vectors(direction)
        return new{T}(direction, halfangleu, halfanglev, numraysu, numraysv, uvec, vvec)
    end

    function RectGrid(halfangleu::T, halfanglev::T, numraysu::Integer, numraysv::Integer) where {T<:Real}
        direction, uvec, vvec = (unitZ3(T), unitX3(T), unitY3(T))
        return new{T}(direction, halfangleu, halfanglev, numraysu, numraysv, uvec, vvec)
    end
end

Base.length(d::RectGrid) = d.nraysu * d.nraysv
function Emitters.generate(d::RectGrid{T}, n::Integer) where {T<:Real}
    direction = d.direction
    uvec = d.uvec
    vvec = d.vvec

    # distributing evenly across the area of the rectangle which subtends the given angle (*not* evenly across the angle range)
    dindex = mod(n, d.nraysu * d.nraysv)
    v = d.nraysv == 1 ? zero(T) : 2 * Integer(floor(dindex / d.nraysu)) / (d.nraysv - 1) - 1.0
    u = d.nraysu == 1 ? zero(T) : 2 * mod(dindex, d.nraysu) / (d.nraysu - 1) - 1.0
    θu = atan(u * tan(d.halfangleu) / 2) * d.halfangleu
    θv = atan(v * tan(d.halfanglev) / 2) * d.halfanglev
    dir = cos(θv) * (cos(θu) * direction + sin(θu) * uvec) + sin(θv) * vvec
    return dir
end

"""
    UniformCone{T} <: AbstractDirectionDistribution{T}

Encapsulates `numsamples` rays sampled uniformly from a cone with max angle θmax.

```julia
UniformCone(direction::Vec3{T}, θmax::T, numsamples::Integer) where {T<:Real}
UniformCone(θmax::T, numsamples::Integer) where {T<:Real}
```
"""
struct UniformCone{T} <: AbstractDirectionDistribution{T}
    direction::Vec3{T}
    θmax::T
    numsamples::Integer
    uvec::Vec3{T}
    vvec::Vec3{T}

    function UniformCone(direction::Vec3{T}, θmax::T, numsamples::Integer) where {T<:Real}
        (uvec, vvec) = get_uv_vectors(direction)
        return new{T}(direction, θmax, numsamples, uvec, vvec)
    end

    function UniformCone(θmax::T, numsamples::Integer) where {T<:Real}
        direction, uvec, vvec = (unitZ3(T), unitX3(T), unitY3(T))
        return new{T}(direction, θmax, numsamples, uvec, vvec)
    end
end

Base.length(d::UniformCone) = d.numsamples
function Emitters.generate(d::UniformCone{T}, ::Integer) where {T<:Real}
    direction = d.direction
    θmax = d.θmax
    uvec = d.uvec
    vvec = d.vvec

    ϕ = rand(T) * 2π
    θ = acos(clamp(one(T) + rand(T) * (cos(θmax) - 1), -one(T), one(T)))
    return normalize(sin(θ) * (cos(ϕ) * uvec + sin(ϕ) * vvec) + cos(θ) * direction)
end

"""
    HexapolarCone{T} <: AbstractDirectionDistribution{T}

Rays are generated by sampling a cone with θmax angle in an hexapolar fashion. The number of rays depends on the requested rings and is computed using the following formula:
`1 + round(Integer, (nrings * (nrings + 1) / 2) * 6)`

```julia
HexapolarCone(direction::Vec3{T}, θmax::T, nrings::Integer) where {T<:Real}
HexapolarCone(θmax::T, nrings::Integer = 3) where {T<:Real}
```
"""
struct HexapolarCone{T} <: AbstractDirectionDistribution{T}
    direction::Vec3{T}
    θmax::T
    nrings::Integer
    uvec::Vec3{T}
    vvec::Vec3{T}

    function HexapolarCone(direction::Vec3{T}, θmax::T, nrings::Integer) where {T<:Real}
        (uvec, vvec) = get_orthogonal_vectors(direction)
        return new{T}(direction, θmax, nrings, uvec, vvec)
    end

    # assume canonical directions
    function HexapolarCone(θmax::T, nrings::Integer = 3) where {T<:Real}
        direction, uvec, vvec = (unitZ3(T), unitX3(T), unitY3(T))
        return new{T}(direction, θmax, nrings, uvec, vvec)
    end
end

Base.length(d::HexapolarCone) = 1 + round(Integer, (d.nrings * (d.nrings + 1) / 2) * 6)
function Emitters.generate(d::HexapolarCone{T}, n::Integer) where {T<:Real}
    dir = d.direction
    θmax = d.θmax
    uvec = d.uvec
    vvec = d.vvec

    n = mod(n, length(d))
    if n == 0
        return normalize(dir)
    else
        t = 1
        ringi = 1
        for i in 1:(d.nrings)
            t += 6 * i
            if n < t
                ringi = i
                break
            end
        end
        ρ = ringi / d.nrings
        pind = n - (t - 6 * ringi)
        
        ϕ = (pind / (6 * ringi)) * 2π
        # elevation calculated as ring fraction multipled by max angle
        θ = acos(clamp(one(T) +  (cos(ρ * θmax) - 1), -one(T), one(T)))
        return normalize(sin(θ) * (cos(ϕ) * uvec + sin(ϕ) * vvec) + cos(θ) * dir)
    end
end

end # module Directions
