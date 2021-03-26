module Spectrum     

import ...OpticSim
using ..Emitters

using DataFrames
using Distributions

abstract type AbstractSpectrum{T<:Real} end

const UNIFORMSHORT = 0.450 #um
const UNIFORMLONG = 0.680 #um

#---------------------------------------
# Uniform Spectrum
#---------------------------------------
struct Uniform{T} <: AbstractSpectrum{T}
    low_end::T
    high_end::T

    # user defined range of spectrum
    function Uniform(low::T, high::T) where {T<:Real}
        return new{T}(low, high)
    end

    # with no specific range we will use the constants' values
    function Uniform(::Type{T} = Float64) where {T<:Real}
        return new{T}(UNIFORMSHORT, UNIFORMLONG)
    end
end
export Uniform

Emitters.generate(s::Uniform{T}) where {T} = (T(1), rand(Distributions.Uniform(s.low_end, s.high_end)))

#---------------------------------------
# Delta Function Spectrum
#---------------------------------------
struct DeltaFunction{T<:Real} <: AbstractSpectrum{T}
    λ::T
end
export DeltaFunction

Emitters.generate(s::DeltaFunction{T}) where {T} = (T(1), s.λ)

#---------------------------------------
# Measured Spectrum
#---------------------------------------

struct Measured{T<:Real} <: AbstractSpectrum{T}
    low_wave_length::Int
    high_wave_length::Int
    wave_length_step::Int
    power_samples::Vector{T}

    function Measured(samples::DataFrame)
        colnames = names(samples)
        @assert "Wavelength" in colnames
        @assert "Power" in colnames
        wavelengths = samples[!, :Wavelength]
        @assert eltype(wavelengths) <: Integer

        power::Vector{T} where {T<:Real} = samples[!, :Power] # no missing values allowed and must be real numbers
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
export Measured

"""expects wavelength in nm not um"""
function spectrumpower(spectrum::Measured{T}, λ::T)::Union{Nothing,T} where {T<:Real}
    λ = λ #convert from um to nm
    if λ < spectrum.low_wave_length || λ > spectrum.high_wave_length
        return nothing
    end

    lowindex = floor(Int, (λ - spectrum.low_wave_length)) ÷ spectrum.wave_length_step + 1

    if lowindex == length(spectrum.power_samples)
        return T(spectrum.power_samples[end])
    else
        highindex = lowindex + 1
        α = mod(λ, spectrum.wave_length_step) / spectrum.wave_length_step

        return T((1 - α) * spectrum.power_samples[lowindex] + α * spectrum.power_samples[highindex])
    end
end


function Emitters.generate(s::Measured{T}) where {T}
    spectrum = s
    λ = rand(Distributions.Uniform{T}(T(spectrum.low_wave_length), T(spectrum.high_wave_length)))
    power = spectrumpower(spectrum, λ)
    if power === nothing        #added this condition because the compiler was allocating if just returned values directly.
        return (nothing, nothing)
    else
        return (power, λ / T(1000)) #convert from nm to um
    end
end


end # module Spectrum 
