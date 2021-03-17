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

#=
Emitters are defined by Pixels and SpatialLayouts. A emitter has a spectrum, and an optical power distribution over the hemisphere. These are
intrinsic physical properties of the emitter.
=#
#Ray distribution types
using Distributions

"""
Each AbstractSpectrum type defines a spectrumSample function which returns a uniformly sampled point from the spectrum of the light source and the power at that wavelength.
"""
abstract type AbstractSpectrum{T<:Real} end
""" Every ray direction distribution must implement the functions direction(a::Source)"""
abstract type AbstractRayDirectionDistribution{T<:Real} end
abstract type AbstractAngularPowerDistribution{T<:Real} end
abstract type AbstractRayOriginDistribution{T<:Real} end

"""
Create an emitter by choosing a type for the spectral distribution, ray origin distribution, ray direction distribution, and angular light power distribution.
For example, to create an emitter with uniform spectrum, rays originating on a regular grid, ray directions distributed around a vector v up to a maximum angle θmax, with Lambertian power distribution call the Source construction with these types and arguments:

a = Source{UniformSpectrum,GridOrigin, ConeDistribution, Lambertian}(UniformSpectrum(),GridOrigin(),ConeDistribution(θmax),Lambertian())
"""
struct Source{T<:Real,S<:AbstractSpectrum{T},O<:AbstractRayOriginDistribution{T},D<:AbstractRayDirectionDistribution{T},P<:AbstractAngularPowerDistribution{T}}
    spectrum::S
    rayorigin::O
    raydirection::D
    power::P

    Source(spectrum::S, origins::O, directions::D, power::P, ::Type{T} = Float64) where {S<:AbstractSpectrum,O<:AbstractRayOriginDistribution,D<:AbstractRayDirectionDistribution,P<:AbstractAngularPowerDistribution,T<:Real} = new{T,S,O,D,P}(spectrum, origins, directions, power)
end
export Source

const UNIFORMSHORT = 0.450 #um
const UNIFORMLONG = 0.680 #um

"""
flat spectrum from 450nm to 680nm
"""
struct UniformSpectrum{T} <: AbstractSpectrum{T} end
export UniformSpectrum
"""returns a tuple (power,wavelength) """
spectrumsample(::Source{T,UniformSpectrum{T},O,D,P}) where {O,D,P,T} = (T(1), rand(Distributions.Uniform(UNIFORMSHORT, UNIFORMLONG)))

struct DeltaFunctionSpectrum{T<:Real} <: AbstractSpectrum{T}
    λ::T
end
export DeltaFunctionSpectrum
spectrumsample(a::Source{T,DeltaFunctionSpectrum{T},O,D,P}) where {O,D,P,T} = (T(1), a.spectrum.λ)

"""
Use measured spectrum to compute emitter power. Create spectrum by reading CSV files.

Evaluate spectrum at arbitrary wavelength with [`spectrumpower`](@ref)
"""
struct MeasuredSpectrum{T<:Real} <: AbstractSpectrum{T}
    lowwavelength::Int
    highwavelength::Int
    wavelengthstep::Int
    powersamples::Vector{T}

    function MeasuredSpectrum(samples::DataFrame)
        colnames = names(samples)
        @assert "Wavelength" in colnames
        @assert "Power" in colnames
        wavelengths = samples[:Wavelength]
        @assert eltype(wavelengths) <: Integer

        power::Vector{T} where {T<:Real} = samples[:Power] # no missing values allowed and must be real numbers
        maxpower = maximum(power)
        T = eltype(power)
        p = sortperm(wavelengths)
        wavelengths = wavelengths[p]
        power = power[p] ./ maxpower #make sure power values have the same sorted permutation as wavelengths normalized with maximum power equal to 1
        λmin = wavelengths[p[1]]
        λmax = wavelengths[p[end]]
        step = wavelengths[2] - wavelengths[1]

        for i in 3:length(wavelengths)
            @assert wavelengths[i] - wavelengths[i - 1] == step  # no missing values allowed and step size must be consistent
        end

        return new{T}(λmin, λmax, step, power)
    end
end
export MeasuredSpectrum

"""expects wavelength in nm not um"""
function spectrumpower(spectrum::MeasuredSpectrum{T}, λ::T)::Union{Nothing,T} where {T<:Real}
    λ = λ #convert from um to nm
    if λ < spectrum.lowwavelength || λ > spectrum.highwavelength
        return nothing
    end

    lowindex = floor(Int, (λ - spectrum.lowwavelength)) ÷ spectrum.wavelengthstep + 1

    if lowindex == length(spectrum.powersamples)
        return T(spectrum.powersamples[end])
    else
        highindex = lowindex + 1
        α = mod(λ, spectrum.wavelengthstep) / spectrum.wavelengthstep

        return T((1 - α) * spectrum.powersamples[lowindex] + α * spectrum.powersamples[highindex])
    end
end

# Every spectrum type must define a function spectrumsample which randomly selects a wavelength and returns the power at that wavelength. For delta function spectra the same wavelength will be returned every time.

function spectrumsample(a::Source{T,MeasuredSpectrum{T},O,D,P}) where {O,D,P,T}
    spectrum = a.spectrum::MeasuredSpectrum{T}
    λ = rand(Distributions.Uniform{T}(T(spectrum.lowwavelength), T(spectrum.highwavelength)))
    power = spectrumpower(spectrum, λ)
    if power === nothing        #added this condition because the compiler was allocating if just returned values directly.
        return (nothing, nothing)
    else
        return (power, λ / T(1000)) #convert from nm to um
    end
end

"""stores the state of ray generation. Keeps track of ray number so you can have a PointOrigin which emits multiple rays. This can be used in the future to do phase summations on the image plane, since only point sources are perfectly spatially coherent."""
struct RayState
    originnumber::Int
    directionnumber::Int

    RayState() = new(1, 1)
    RayState(originnum::Int, directionnum::Int) = new(originnum, directionnum)
end

struct PointOrigin{T<:Real} <: AbstractRayOriginDistribution{T}
    position::SVector{3,T}
end
export PointOrigin

numsamples(a::PointOrigin{T}) where {T} = 1

origin(a::Source{T,S,PointOrigin{T},D,P}, ::Int) where {S,D,P,T} = a.rayorigin.position

struct UniformRandomOrigin{T} <: AbstractRayOriginDistribution{T}
    width::T
    height::T
    numsamples::Int
end
export UniformRandomOrigin

numsamples(a::UniformRandomOrigin{T}) where {T} = a.numsamples

function origin(a::Source{T,S,UniformRandomOrigin{T},D,P}, ::Int) where {S,D,P,T}
    origin = a.origin
    x = rand(Distributions.Uniform(-origin.width / 2.0, origin.width / 2.0))
    y = rand(Distributions.Uniform(-origin.height / 2.0, origin.height / 2.0))
    z = T(0)
end

struct GridOrigin{T} <: AbstractRayOriginDistribution{T}
    width::T
    height::T
    usamples::Int
    vsamples::Int
    ustep::T
    vstep::T

    GridOrigin(width::T, height::T, usamples::Int, vsamples::Int) where {T<:Real} = new{T}(width, height, usamples, vsamples, width / (usamples - 1), height / (vsamples - 1))
end
export GridOrigin

numsamples(a::GridOrigin{T}) where {T<:Real} = a.usamples * a.vsamples

function origin(a::Source{T,S,GridOrigin{T},D,P}, originnumber::Int) where {S,D,P,T}
    @assert originnumber <= numsamples(a.rayorigin)
    origin = a.rayorigin

    u = originnumber ÷ origin.vsamples
    v = mod(originnumber, origin.usamples)
    return SVector{3,T}(T(-origin.width / 2.0 + origin.ustep * u), T(-origin.height / 2.0 + origin.vstep * v), T(0))
end

struct HexapolarOrigin{T} <: AbstractRayOriginDistribution{T}
    position::SVector{3,T}
    direction::SVector{3,T}
    uvec::SVector{3,T}
    vvec::SVector{3,T}
    halfsizeu::T
    halfsizev::T
    nrings::Int
    function HexapolarOrigin{T}(nrings::Int, halfsizeu::T, halfsizev::T; position::SVector{3,T} = SVector{3,T}(0.0, 0.0, 0.0), direction::SVector{3,T} = SVector{3,T}(0.0, 0.0, -1.0), rotationvec::SVector{3,T} = SVector{3,T}(0.0, 1.0, 0.0)) where {T<:Real}
        new{T}(position, getuvvecs(direction, rotationvec)..., halfsizeu, halfsizev, nrings)
    end
end

struct ConstantRayDirection{T} <: AbstractRayDirectionDistribution{T}
    direction::SVector{3,T}
end
export ConstantRayDirection

numsamples(a::ConstantRayDirection{T}) where {T} = 1

direction(a::Source{T,S,O,ConstantRayDirection{T},P}, directionnumber::Int) where {S,O,P,T} = a.raydirection.direction

struct ConeDistribution{T} <: AbstractRayDirectionDistribution{T}
    θmax::T
    numsamples::Int
end
export ConeDistribution

numsamples(a::ConeDistribution{T}) where {T} = a.numsamples

"""
Generates a unit vector pointing somewhere within the cone with half angle `θmax` around `direction` which is the normal to the emitter surface.
"""
function direction(a::Source{T,S,O,ConeDistribution{T},P}, directionnumber::Int) where {S,O,P,T<:Real}
    direction = SVector{3,T}(T(0), T(0), T(1))
    θmax = a.raydirection.θmax
    uvec = SVector{3,T}(T(1), T(0), T(0))
    vvec = SVector{3,T}(T(0), T(1), T(0))

    ϕ = rand(T) * 2π
    θ = NaNsafeacos(one(T) + rand(T) * (cos(θmax) - 1))
    return normalize(sin(θ) * (cos(ϕ) * uvec + sin(ϕ) * vvec) + cos(θ) * direction)
end



Base.length(a::Source{T,S,O,D,P}) where {T,S,O,D,P} = numsamples(a.rayorigin) * numsamples(a.raydirection)

"""raysample is the function that can couple origin and direction generation if necessary. The default function couples them in a simple way but more complex coupling should be possible. For each origin numsamples(a.raydirection) direction samples are taken with identical origin. Then the origin number is incremented. This repeats till all rays have been generated. The origin and direction functions receive an integer indicating the origin or direction number so regular patterns such as rectangular and hexapolar grids can be generated properly."""
function raysample(a::Source{T,S,O,D,P}, r::RayState) where {T,S,O,D,P}
    newray = Ray(origin(a, r.originnumber), direction(a, r.directionnumber))
    dirnum = r.directionnumber + 1
    if dirnum > numsamples(a.raydirection)
        return (newray, RayState(r.originnumber + 1, 1))  #reset counter for direction number and move on to next origin
    else
        return (newray, RayState(r.originnumber, dirnum))
    end

end

struct Lambertian{T} <: AbstractAngularPowerDistribution{T} end
export Lambertian

directionpower(::Source{T,S,O,D,Lambertian{T}}, ::SVector{3,T}) where {S,O,D,T} = T(1)

function Base.iterate(a::Source)::Tuple{OpticalRay,RayState}
    state = RayState()
    generateray(a, state)
end

function Base.iterate(a::Source, b::RayState)::Union{Nothing,Tuple{OpticalRay,RayState}}
    if b.originnumber > numsamples(a.rayorigin)
        return nothing
    else
        return generateray(a, b)
    end
end

# """generates an optical ray with origin and spectrum distribution dictated by OriginDistribution and Spectrum type. Ray is generated in the emitter local coordinate frame."""
function generateray(a::Source, state::RayState, sourcenumber = 0)
    spectralpower, λ = spectrumsample(a)
    ray, newstate = raysample(a, state)
    optray = OpticalRay(ray, spectralpower * directionpower(a, direction(ray)), λ, sourcenum = sourcenumber)
    return (optray, newstate)
end
