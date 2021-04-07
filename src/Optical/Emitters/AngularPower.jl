module AngularPower

using ....OpticSim, ...Geometry
using ...Emitters
using LinearAlgebra

abstract type AbstractAngularPowerDistribution{T<:Real} end

#---------------------------------------
# Lambertian
#---------------------------------------
"""
Ray power is unaffected by angle.
"""
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
"""
    Cosine(cosine_exp::T = one(T)) where {T<:Real}

Cosine power distribution. Ray power is calculated by:

`power = power * (cosine_angle ^ cosine_exp)`
"""
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
"""
    Gaussian(gaussianu::T, gaussianv::T) where {T<:Real}

GGaussian power distribution. Ray power is calculated by:

`power = power * exp(-(gaussianu * l^2 + gaussianv * m^2))`
where l and m are the cos_angles between the two axes respectivly.
"""
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
