module AngularPower

import OpticSim
using LinearAlgebra
using ..Emitters
using ..Geometry

abstract type AbstractAngularPowerDistribution{T<:Real} end

#---------------------------------------
# Lambertian
#---------------------------------------

struct Lambertian{T} <: AbstractAngularPowerDistribution{T} 
    Lambertian(::Type{T} = Float64) where {T<:Real} = new{T}()
end

function Emitters.apply(d::Lambertian{T}, Tr::Transform{T}, power::T, ray::OpticSim.Ray{T,3}) where {T}
    return power
end
export Lambertian

#---------------------------------------
# Cosine Power Distribution
#---------------------------------------

struct Cosine{T} <: AbstractAngularPowerDistribution{T} 
    cosine_exp::T

    function Cosine(cosine_exp::T = one(T)) where {T<:Real}
        new{T}(cosine_exp)
    end
end
export Cosine

# returning ray power
function Emitters.apply(d::Cosine{T}, Tr::Transform{T}, power::T, ray::OpticSim.Ray{T,3}) where {T}
    cosanglebetween = dot(OpticSim.direction(ray), forward(Tr))
    power = power * cosanglebetween ^ d.cosine_exp
    return power
end

#---------------------------------------
# Gaussian Power Distribution
#---------------------------------------

struct Gaussian{T} <: AbstractAngularPowerDistribution{T} 
    gaussianu::T
    gaussianv::T

    function Gaussian(gaussianu::T, gaussianv::T) where {T<:Real}
        new{T}(gaussianu, gaussianv)
    end
end
export Gaussian

# returning ray power
function Emitters.apply(d::Gaussian{T}, Tr::Transform{T}, power::T, ray::OpticSim.Ray{T,3}) where {T}

    l = dot(OpticSim.direction(ray), right(Tr))
    m = dot(OpticSim.direction(ray), up(Tr))
    power = power * exp(-(d.gaussianu * l^2 + d.gaussianv * m^2))

    return power
end


end # module Angular Power
