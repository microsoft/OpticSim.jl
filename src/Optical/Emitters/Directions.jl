module Directions

using ....OpticSim, ...Geometry
using ...Emitters
using LinearAlgebra

abstract type AbstractDirectionDistribution{T<:Real} end

#---------------------------------------
# Direction Distrubution Common Utilities
#---------------------------------------

export generate

Base.iterate(a::AbstractDirectionDistribution, state = 1) = state > length(a) ? nothing : (generate(a, state - 1), state + 1)
Base.getindex(a::AbstractDirectionDistribution, index) = generate(a, index)
Base.firstindex(a::AbstractDirectionDistribution) = 0
Base.lastindex(a::AbstractDirectionDistribution) = length(a) - 1
Base.copy(a::AbstractDirectionDistribution) = a # most don't have any heap allocated stuff so don't really need copying


#---------------------------------------
# Constant Ray Direction
#---------------------------------------
"""
Encapsulates a single ray direction, where the default direction is unitZ3 [0, 0, 1].

```julia
Constant(direction::Vec3{T}) where {T<:Real}
Constant(::Type{T} = Float64) where {T<:Real}
```
"""
struct Constant{T} <: AbstractDirectionDistribution{T}
    direction::Vec3{T}

    function Constant(direction::Vec3{T}) where {T<:Real}
        return new{T}(direction)
    end
    
    function Constant(::Type{T} = Float64) where {T<:Real}
        direction = unitZ3(T)
        return new{T}(direction)
    end
end
export Constant

Base.length(d::Constant{T}) where {T} = 1
function Emitters.generate(d::Constant{T}, n::Int) where {T<:Real}
    return d.direction
end

#---------------------------------------
# Grid  Ray Distribution
#---------------------------------------
"""
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
    nraysu::Int
    nraysv::Int
    uvec::Vec3{T}
    vvec::Vec3{T}

    function RectGrid(direction::Vec3{T}, halfangleu::T, halfanglev::T, numraysu::Int, numraysv::Int) where {T<:Real}
        (uvec, vvec) = get_uv_vectors(direction)
        return new{T}(direction, halfangleu, halfanglev, numraysu, numraysv, uvec, vvec)
    end

    function RectGrid(halfangleu::T, halfanglev::T, numraysu::Int, numraysv::Int) where {T<:Real}
        direction, uvec, vvec = (unitZ3(T), unitX3(T), unitY3(T))
        return new{T}(direction, halfangleu, halfanglev, numraysu, numraysv, uvec, vvec)
    end
end
export RectGrid

Base.length(d::RectGrid{T}) where {T} = d.nraysu * d.nraysv
function Emitters.generate(d::RectGrid{T}, n::Int) where {T<:Real}
    direction = d.direction
    uvec = d.uvec
    vvec = d.vvec

    # distributing evenly across the area of the rectangle which subtends the given angle (*not* evenly across the angle range)
    dindex = mod(n, d.nraysu * d.nraysv)
    v = d.nraysv == 1 ? zero(T) : 2 * Int(floor(dindex / d.nraysu)) / (d.nraysv - 1) - 1.0
    u = d.nraysu == 1 ? zero(T) : 2 * mod(dindex, d.nraysu) / (d.nraysu - 1) - 1.0
    θu = atan(u * tan(d.halfangleu) / 2) * d.halfangleu
    θv = atan(v * tan(d.halfanglev) / 2) * d.halfanglev
    dir = cos(θv) * (cos(θu) * direction + sin(θu) * uvec) + sin(θv) * vvec
    return dir
end


#---------------------------------------
# Cone Ray Distribution
#---------------------------------------
"""
Encapsulates a numsamples rays sampled uniformlly from a cone with max angle θmax.

```julia
UniformCone(direction::Vec3{T}, θmax::T, numsamples::Int) where {T<:Real}
UniformCone(θmax::T, numsamples::Int) where {T<:Real}
```
"""
struct UniformCone{T} <: AbstractDirectionDistribution{T}
    direction::Vec3{T}
    θmax::T
    numsamples::Int
    uvec::Vec3{T}
    vvec::Vec3{T}

    function UniformCone(direction::Vec3{T}, θmax::T, numsamples::Int) where {T<:Real}
        (uvec, vvec) = get_uv_vectors(direction)
        return new{T}(direction, θmax, numsamples, uvec, vvec)
    end

    function UniformCone(θmax::T, numsamples::Int) where {T<:Real}
        direction, uvec, vvec = (unitZ3(T), unitX3(T), unitY3(T))
        return new{T}(direction, θmax, numsamples, uvec, vvec)
    end
end
export UniformCone

Base.length(d::UniformCone{T}) where {T} = d.numsamples
function Emitters.generate(d::UniformCone{T}, n::Int) where {T<:Real}
    direction = d.direction
    θmax = d.θmax
    uvec = d.uvec
    vvec = d.vvec

    ϕ = rand(T) * 2π
    θ = acos(clamp(one(T) + rand(T) * (cos(θmax) - 1), -one(T), one(T)))
    return normalize(sin(θ) * (cos(ϕ) * uvec + sin(ϕ) * vvec) + cos(θ) * direction)
end

#----------------------------------------------------------------------------
# Cone Ray Hexapolar Distribution
#-----------------------------------------------------------------------------
"""
Rays are generated by sampling a cone with θmax angle in an hexapolar fasion. The number of rays depends on the requested rings and is computed using the following formula:
`1 + round(Int, (nrings * (nrings + 1) / 2) * 6)`

```julia
HexapolarCone(direction::Vec3{T}, θmax::T, nrings::Int) where {T<:Real}
HexapolarCone(θmax::T, nrings::Int = 3) where {T<:Real}
```
"""
struct HexapolarCone{T} <: AbstractDirectionDistribution{T}
    direction::Vec3{T}
    θmax::T
    nrings::Int
    uvec::Vec3{T}
    vvec::Vec3{T}

    function HexapolarCone(direction::Vec3{T}, θmax::T, nrings::Int) where {T<:Real}
        (uvec, vvec) = get_orthogonal_vectors(direction)
        return new{T}(direction, θmax, nrings, uvec, vvec)
    end

    # assume canonical directions
    function HexapolarCone(θmax::T, nrings::Int = 3) where {T<:Real}
        direction, uvec, vvec = (unitZ3(T), unitX3(T), unitY3(T))
        return new{T}(direction, θmax, nrings, uvec, vvec)
    end
end
export HexapolarCone


Base.length(d::HexapolarCone{T}) where {T} = 1 + round(Int, (d.nrings * (d.nrings + 1) / 2) * 6)
function Emitters.generate(d::HexapolarCone{T}, n::Int) where {T<:Real}
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
