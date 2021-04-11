module AngularPower
export Lambertian, Cosine, Gaussian

using ....OpticSim
using ...Emitters
using ...Geometry
using LinearAlgebra

abstract type AbstractAngularPowerDistribution{T<:Real} end

"""
    Lambertian{T} <: AbstractAngularPowerDistribution{T} 

Ray power is unaffected by angle.
"""
struct Lambertian{T} <: AbstractAngularPowerDistribution{T}
    Lambertian(::Type{T} = Float64) where {T<:Real} = new{T}()
end

Emitters.apply(::Lambertian, ::Transform, power::Real, ::OpticSim.Ray{<:Real,3}) = power

"""
    Cosine{T} <: AbstractAngularPowerDistribution{T} 

Cosine power distribution. Ray power is calculated by:

`power = power * (cosine_angle ^ cosine_exp)`
"""
struct Cosine{T} <: AbstractAngularPowerDistribution{T} 
    cosine_exp::T

    function Cosine(cosine_exp::T = one(T)) where {T<:Real}
        new{T}(cosine_exp)
    end
end

function Emitters.apply(d::Cosine, tr::Transform, power::Real, ray::OpticSim.Ray{<:Real,3})
    cosanglebetween = dot(OpticSim.direction(ray), forward(tr))
    return power * cosanglebetween ^ d.cosine_exp
end

"""
    Gaussian{T} <: AbstractAngularPowerDistribution{T} 

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

function Emitters.apply(d::Gaussian, tr::Transform, power::Real, ray::OpticSim.Ray{<:Real,3})
    l = dot(OpticSim.direction(ray), right(tr))
    m = dot(OpticSim.direction(ray), up(tr))
    return power * exp(-(d.gaussianu * l^2 + d.gaussianv * m^2))
end

end # module Angular Power
