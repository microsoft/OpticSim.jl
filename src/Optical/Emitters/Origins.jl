module Origins

using ....OpticSim, ...Geometry
using ...Emitters
using LinearAlgebra
using Distributions


abstract type AbstractOriginDistribution{T<:Real} end

#---------------------------------------
# Origin Distrubution Common Utilities
#---------------------------------------

Base.iterate(a::AbstractOriginDistribution, state = 1) = state > length(a) ? nothing : (generate(a, state - 1), state + 1)
Base.getindex(a::AbstractOriginDistribution, index) = generateray(a, index)
Base.firstindex(a::AbstractOriginDistribution) = 0
Base.lastindex(a::AbstractOriginDistribution) = length(a) - 1
Base.copy(a::AbstractOriginDistribution) = a # most don't have any heap allocated stuff so don't really need copying

#--------------------------------------
# Point Origin
#--------------------------------------
"""
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
export Point

Base.length(o::Point{T}) where {T} = 1
Emitters.visual_size(o::Point{T}) where {T} = 1

function Emitters.generate(o::Point{T}, n::Int) where {T<:Real}
    return o.origin
end

#-------------------------------------
# Random Rectangle Origin
#-------------------------------------
"""
Encapsulates a uniformly sampled rectangle with user defined number of samples.

```julia
RectUniform(width::T, height::T, count::Int) where {T<:Real}
```
"""
struct RectUniform{T} <: AbstractOriginDistribution{T}
    width::T
    height::T
    samples_count::Int

    function RectUniform(width::T, height::T, count::Int) where {T<:Real}
        return new{T}(width, height, count)
    end
end
export RectUniform


Base.length(o::RectUniform{T}) where {T} = o.samples_count
Emitters.visual_size(o::RectUniform{T}) where {T} = max(o.width, o.height)

# generate origin on the agrid
function Emitters.generate(o::RectUniform{T}, n::Int) where {T<:Real}
    n = mod(n, length(o))
    u = rand(Distributions.Uniform(T(-1), T(1)))
    v = rand(Distributions.Uniform(T(-1), T(1)))
    return zero(Vec3{T}) + ((o.width / 2) * u * unitX3(T)) + ((o.height/2) * v * unitY3(T))
end


#-------------------------------------
# Grid Rectangle Origin
#-------------------------------------
"""
Encapsulates a rectangle sampled in a grid fashion.

```julia
RectGrid(width::T, height::T, usamples::Int, vsamples::Int) where {T<:Real} 
```
"""
struct RectGrid{T} <: AbstractOriginDistribution{T}
    width::T
    height::T
    usamples::Int
    vsamples::Int
    ustep::T
    vstep::T

    function RectGrid(width::T, height::T, usamples::Int, vsamples::Int) where {T<:Real} 
        return new{T}(width, height, usamples, vsamples, width / (usamples - 1), height / (vsamples - 1))
    end
end
export RectGrid

Base.length(o::RectGrid{T}) where {T} = o.usamples * o.vsamples
Emitters.visual_size(o::RectGrid{T}) where {T} = max(o.width, o.height)

# generate origin on the agrid
function Emitters.generate(o::RectGrid{T}, n::Int) where {T<:Real}
    n = mod(n, length(o))
    v = o.vsamples == 1 ? zero(T) : 2 * Int(floor(n / o.usamples)) / (o.vsamples - 1) - 1.0
    u = o.usamples == 1 ? zero(T) : 2 * mod(n, o.usamples) / (o.usamples - 1) - 1.0
    return zeros(Vec3{T}) + ((o.width / 2) * u * unitX3(T)) + ((o.height/2) * v * unitY3(T))
end


#-------------------------------------
# Hexapolar Origin
#-------------------------------------
"""
Encapsulates an ellipse (or a circle where halfsizeu=halfsizev) sampled in an hexapolar fasion (rings)

```julia
Hexapolar(nrings::Int, halfsizeu::T, halfsizev::T) where {T<:Real} 
```
"""
struct Hexapolar{T} <: AbstractOriginDistribution{T}
    halfsizeu::T
    halfsizev::T
    nrings::Int

    function Hexapolar(nrings::Int, halfsizeu::T, halfsizev::T) where {T<:Real} 
        return new{T}(halfsizeu, halfsizev, nrings)
    end
end
export Hexapolar

Base.length(o::Hexapolar{T}) where {T} = 1 + round(Int, (o.nrings * (o.nrings + 1) / 2) * 6)
Emitters.visual_size(o::Hexapolar{T}) where {T} = max(o.halfsizeu*2, o.halfsizev*2)

function Emitters.generate(o::Hexapolar{T}, n::Int) where {T<:Real}
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
