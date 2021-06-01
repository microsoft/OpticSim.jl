# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

module Spectrum     
export Uniform, DeltaFunction, Measured

using ....OpticSim
using ...Emitters
using DataFrames
using Distributions
import Unitful: Length, ustrip
using Unitful.DefaultSymbols

const UNIFORMSHORT = 0.450 #um
const UNIFORMLONG = 0.680 #um

abstract type AbstractSpectrum{T<:Real} end

"""
    Uniform{T} <: AbstractSpectrum{T}

Encapsulates a flat spectrum range which is sampled uniformly. Unless stated diferrently, the range used will be 450nm to 680nm.

```julia
Uniform(low_end::T, high_end::T) where {T<:Real}
Uniform(::Type{T} = Float64) where {T<:Real}
```
"""
struct Uniform{T} <: AbstractSpectrum{T}
    low_end::T
    high_end::T 

    # user defined range of spectrum
    function Uniform(low_end::T, high_end::T) where {T<:Real}
        return new{T}(low_end, high_end)
    end

    # with no specific range we will use the constants' values
    function Uniform(::Type{T} = Float64) where {T<:Real}
        return new{T}(UNIFORMSHORT, UNIFORMLONG)
    end
end

Emitters.generate(s::Uniform{T}) where {T<:Real} = (one(T), rand(Distributions.Uniform(s.low_end, s.high_end)))

"""
    DeltaFunction{T} <: AbstractSpectrum{T}

Encapsulates a constant spectrum.

```julia
DeltaFunction{T<:Real}
```
"""
struct DeltaFunction{T} <: AbstractSpectrum{T}
    λ::T
end

DeltaFunction(λ::Length) = DeltaFunction{Float64}(ustrip(μm, λ))

Emitters.generate(s::DeltaFunction{T}) where {T<:Real} = (one(T), s.λ)

"""
    Measured{T} <: AbstractSpectrum{T}

Encapsulates a measured spectrum to compute emitter power. Create spectrum by reading CSV files.
Evaluate spectrum at arbitrary wavelength with [`spectrumpower`](@ref) (**more technical details coming soon**)

```julia
Measured(samples::DataFrame)
```
"""
struct Measured{T} <: AbstractSpectrum{T}
    low_wave_length::Integer
    high_wave_length::Integer
    wave_length_step::Integer
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

function Emitters.generate(spectrum::Measured{T}) where {T<:Real}
    λ = rand(Distributions.Uniform(convert(T, spectrum.low_wave_length), convert(T, spectrum.high_wave_length)))
    power = spectrumpower(spectrum, λ)
    if power === nothing        #added this condition because the compiler was allocating if just returned values directly.
        return (nothing, nothing)
    else
        return (power, λ / convert(T, 1000)) #convert from nm to um
    end
end

"""expects wavelength in nm not um"""
function spectrumpower(spectrum::Measured{T}, λ::T) where {T<:Real}
    if λ < spectrum.low_wave_length || λ > spectrum.high_wave_length
        return nothing
    end

    lowindex = floor(Integer, (λ - spectrum.low_wave_length)) ÷ spectrum.wave_length_step + 1

    if lowindex == length(spectrum.power_samples)
        return convert(T, spectrum.power_samples[end])
    else
        highindex = lowindex + 1
        α = mod(λ, spectrum.wave_length_step) / spectrum.wave_length_step

        return convert(T, (1 - α) * spectrum.power_samples[lowindex] + α * spectrum.power_samples[highindex])
    end
end

end # module Spectrum 
