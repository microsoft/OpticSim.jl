# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

module Sources
export Source, CompositeSource

using ....OpticSim
using ...Emitters
using ...Geometry
using ..Spectrum
using ..Directions
using ..Origins
using ..AngularPower
using LinearAlgebra
using Distributions

abstract type AbstractSource{T} <: OpticSim.OpticalRayGenerator{T} end

Base.getindex(s::AbstractSource, index) = generate(s, index)
Base.firstindex(s::AbstractSource) = 0
Base.lastindex(s::AbstractSource) = length(s) - 1
Base.copy(a::AbstractSource) = a # most don't have any heap allocated stuff so don't really need copying

""" Simple ray generator that takes a vector of OpticalRay objects."""
struct RayListSource{T,N} <: AbstractSource{T}
    rays::Vector{OpticalRay{T,N}}
end

Base.length(raylist::RayListSource) = length(raylist.rays)
Base.iterate(raylist::RayListSource) = length(raylist.rays) == 0 ? nothing : (raylist.rays[1],1)
Base.iterate(raylist::RayListSource, state) = state > length(raylist.rays) ? nothing : (raylist.rays[state+1],state+1)

"""
    Source{T<:Real, Tr<:Transform{T}, S<:Spectrum.AbstractSpectrum{T}, O<:Origins.AbstractOriginDistribution{T}, D<:Directions.AbstractDirectionDistribution{T}, P<:AngularPower.AbstractAngularPowerDistribution{T}} <: AbstractSource{T}

This data-type represents the basic emitter (Source), which is a combination of a Spectrum, Angular Power Distribution, Origins and Directions distibution and a 3D Transform.

```julia
Source(::Type{T} = Float64;
       transform::Tr = Transform(),
       spectrum::S = Spectrum.Uniform(),
       origins::O = Origins.Point(),
       directions::D = Directions.Constant(),
       power::P = AngularPower.Lambertian(),
       sourcenum::Int64 = 0) where {
            Tr<:Transform,
            S<:Spectrum.AbstractSpectrum,
            O<:Origins.AbstractOriginDistribution,
            D<:Directions.AbstractDirectionDistribution,
            P<:AngularPower.AbstractAngularPowerDistribution,
            T<:Real}

Source(transform::Tr, spectrum::S, origins::O, directions::D, power::P, ::Type{T} = Float64; sourcenum::Int64 = 0) where {   
            Tr<:Transform,
            S<:Spectrum.AbstractSpectrum,
            O<:Origins.AbstractOriginDistribution,
            D<:Directions.AbstractDirectionDistribution,
            P<:AngularPower.AbstractAngularPowerDistribution,
            T<:Real}
```
"""
struct Source{
    T<:Real, 
    Tr<:Transform{T},
    S<:Spectrum.AbstractSpectrum{T},
    O<:Origins.AbstractOriginDistribution{T},
    D<:Directions.AbstractDirectionDistribution{T},
    P<:AngularPower.AbstractAngularPowerDistribution{T}
} <: AbstractSource{T}
    transform::Tr
    spectrum::S
    origins::O
    directions::D
    power_distribution::P
    sourcenum::Int64

    function Source(
        ::Type{T} = Float64;
        transform::Tr = Transform(T),
        spectrum::S = Spectrum.Uniform(T),
        origins::O = Origins.Point(T),
        directions::D = Directions.Constant(T),
        power::P = AngularPower.Lambertian(T),
        sourcenum::Int64 = 0
    ) where {
        Tr<:Transform,
        S<:Spectrum.AbstractSpectrum,
        O<:Origins.AbstractOriginDistribution,
        D<:Directions.AbstractDirectionDistribution,
        P<:AngularPower.AbstractAngularPowerDistribution,
        T<:Real
    }
        new{T, Tr, S, O, D, P}(transform, spectrum, origins, directions, power, sourcenum)
    end

    function Source(
        transform::Tr,
        spectrum::S,
        origins::O,
        directions::D,
        power::P,
        ::Type{T} = Float64;
        sourcenum::Int64 = 0
    ) where {
        Tr<:Transform,
        S<:Spectrum.AbstractSpectrum,
        O<:Origins.AbstractOriginDistribution,
        D<:Directions.AbstractDirectionDistribution,
        P<:AngularPower.AbstractAngularPowerDistribution,
        T<:Real
    }
         new{T, Tr, S, O, D, P}(transform, spectrum, origins, directions, power, sourcenum)
    end
end

# used to not generate new origin points if we can use last points - mainly to keep consistency of origin generation when randomness is involved
struct SourceGenerationState{T<:Real}
    n::Int64
    last_index_used::Int64
    last_point_generated::Vec3{T}
end

Base.length(s::Source) = length(s.origins) * length(s.directions)
Base.iterate(s::Source{T}) where {T<:Real} = iterate(s, SourceGenerationState(1, -1, zeros(Vec3{T})))
Base.iterate(s::Source, state::SourceGenerationState) = state.n > length(s) ? nothing : generate(s, state)

function Emitters.generate(s::Source{T}, n::Int64) where {T<:Real}
    return generate(s, SourceGenerationState(n+1, -1, zeros(Vec3{T})))[1]
end

function Emitters.generate(s::Source{T}, state::SourceGenerationState{T} = SourceGenerationState(-1, -1, zeros(Vec3{T}))) where {T<:Real}
    # @info "New State Generation"
    n::Int64 = state.n - 1
    origin_index = floor(Int64, (n / length(s.directions))) 
    direction_index = mod(n, length(s.directions)) 

    if (origin_index == state.last_index_used)
        o = state.last_point_generated
    else
        o = generate(s.origins, origin_index)
        o = s.transform * o
    end

    # o = generate(s.origins, origin_index)
    d = generate(s.directions, direction_index)
    
    # rotate the direction according to the orientation of the ray-originator
    d = rotation(s.transform) * d

    ray = OpticSim.Ray(o, d)
    power, wavelength = generate(s.spectrum)
    
    power = apply(s.power_distribution, s.transform, power, ray)

    # return (OpticSim.Ray(o, d), SourceGenerationState(state.n+1, origin_index, o))
    return (OpticSim.OpticalRay(ray, power, wavelength; sourcenum=s.sourcenum), SourceGenerationState(state.n+1, origin_index, o))
end

"""
    CompositeSource{T} <: AbstractSource{T}

This data-type represents the composite emitter (Source) which is constructed with a list of basic or composite emitters and a 3D Transform.

```julia
CompositeSource(transform::Transform{T}, sources::Vector{<:AbstractSource}) where {T<:Real} 
```
"""
struct CompositeSource{T} <: AbstractSource{T}
    transform::Transform{T}
    sources::Vector{<:AbstractSource}

    uniform_length::Int
    total_length::Int
    start_indexes::Vector{Int}

    function CompositeSource(transform::Transform{T}, sources::Vector{<:AbstractSource}) where {T<:Real}
        lens = [length(src) for src in sources]
        size1 = findmin(lens)[1]
        size2 = findmax(lens)[1]
        if (size1 == size2)
            uniform_length = size1
            start_indexes = Vector{Int}()
            total_length = length(sources) * uniform_length
            # @info "Uniform Length $size1 total=$total_length"
        else
            uniform_length = -1

            start_indexes = Vector{Int}()
            start = 0
            for s in sources
                push!(start_indexes, start)
                start = start + length(s)
            end
            push!(start_indexes, start)
            # @info start_indexes
            total_length = start 
        end

        new{T}(transform, sources, uniform_length, total_length, start_indexes)
    end
end

Base.length(s::CompositeSource) = s.total_length
Base.iterate(s::CompositeSource{T}) where {T<:Real} = iterate(s, SourceGenerationState(1, -1, zeros(Vec3{T})))
Base.iterate(s::CompositeSource, state::SourceGenerationState) = state.n > length(s) ? nothing : generate(s, state)

function Emitters.generate(s::CompositeSource{T}, n::Int64) where {T<:Real}
    return generate(s, SourceGenerationState(n+1, -1, zeros(Vec3{T})))[1]
end

function Emitters.generate(s::CompositeSource{T}, state::SourceGenerationState{T} = SourceGenerationState(-1, -1, zeros(Vec3{T}))) where {T<:Real}
    # @info "New Composite State Generation"
    n::Int64 = state.n - 1

    if (s.uniform_length != -1)
        source_index = floor(Int64, (n / s.uniform_length)) 
        ray_index = mod(n, s.uniform_length) 
        # @info "New Composite State Generation: source=$source_index  ray=$ray_index  ($(length(s.sources)))"
    else
        n = mod(n, s.total_length)

        ## perform a binary search to find out which source can generate the spesific ray
        low = 1
        high = length(s.start_indexes) - 1    # last cell in this vector is a dummy cell containing the index after the last ray
        source_index = -1
        while (low <= high) 
            mid = floor(Int64, (low + high) / 2)

            if (n < s.start_indexes[mid])
                high = mid - 1
            elseif (n >= s.start_indexes[mid+1])
                low = mid + 1
            else
                source_index = mid
                break;
            end
        end
        @assert source_index != -1
        # @info "Bin Search n=$n source_index=$source_index"
        source_index = source_index - 1 
        ray_index = n - s.start_indexes[source_index+1]
    end

    ray, new_state = generate(s.sources[source_index+1], SourceGenerationState(ray_index+1, state.last_index_used, state.last_point_generated))
    ray = s.transform * ray                 
    return (ray, SourceGenerationState(state.n+1, new_state.last_index_used, new_state.last_point_generated))
end

end # module Sources
