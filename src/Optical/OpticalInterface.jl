# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

"""
    OpticalInterface{T<:Real}

Any subclass of OpticalInterface **must** implement the following:

```julia
processintersection(opticalinterface::OpticalInterface{T}, point::SVector{N,T}, normal::SVector{N,T}, incidentray::OpticalRay{T,N}, temperature::T, pressure::T, ::Bool, firstray::Bool = false) -> Tuple{SVector{N,T}, T, T}
```

See documentation for [`processintersection`](@ref) for details.

These methods are also commonly implemented, but not essential:
```julia
insidematerialid(i::OpticalInterface{T}) -> OpticSim.GlassCat.AbstractGlass
outsidematerialid(i::OpticalInterface{T}) -> OpticSim.GlassCat.AbstractGlass
reflectance(i::OpticalInterface{T}) -> T
transmission(i::OpticalInterface{T}) -> T
```
"""
abstract type OpticalInterface{T<:Real} end
export OpticalInterface

"""
Valid modes for deterministic raytracing
"""
@enum InterfaceMode Reflect Transmit ReflectOrTransmit
export InterfaceMode, Reflect, Transmit, ReflectOrTransmit

######################################################################

"""
    NullInterface{T} <: OpticalInterface{T}

Interface which will be ignored totally by any rays, used only in construction of CSG objects.

```julia
NullInterface(T = Float64)
NullInterface{T}()
```
"""
struct NullInterface{T} <: OpticalInterface{T}
    NullInterface(::Type{T} = Float64) where {T<:Real} = new{T}()
    NullInterface{T}() where {T<:Real} = new{T}()
end

insidematerialid(::NullInterface{T}) where {T<:Real} = glassid(OpticSim.GlassCat.Air)
outsidematerialid(::NullInterface{T}) where {T<:Real} = glassid(OpticSim.GlassCat.Air)
reflectance(::NullInterface{T}) where {T<:Real} = zero(T)
transmission(::NullInterface{T}) where {T<:Real} = one(T)

export NullInterface

######################################################################

"""
    ParaxialInterface{T} <: OpticalInterface{T}

Interface describing an idealized planar lens, i.e. one that is thin and with no aberrations.

**In general this interface should not be constructed directly, the `ParaxialLensEllipse` and `ParaxialLensRect` functions should be used to create a [`ParaxialLens`](@ref) object directly.**

```julia
ParaxialInterface(focallength::T, centroid::SVector{3,T}, outsidematerial::Y)
```
"""
struct ParaxialInterface{T} <: OpticalInterface{T}
    focallength::T
    outsidematerial::GlassID
    centroid::SVector{3,T}
    function ParaxialInterface(focallength::T, centroid::SVector{3,T}, outsidematerial::Y) where {Y<:OpticSim.GlassCat.AbstractGlass,T<:Real}
        return new{T}(focallength, glassid(outsidematerial), centroid)
    end
end
export ParaxialInterface

function Base.show(io::IO, a::ParaxialInterface{R}) where {R<:Real}
    print(io, "ParaxialInterface($(a.focallength), $(glassname(a.outsidematerial)))")
end

insidematerialid(a::ParaxialInterface{T}) where {T<:Real} = a.outsidematerial
outsidematerialid(a::ParaxialInterface{T}) where {T<:Real} = a.outsidematerial
reflectance(::ParaxialInterface{T}) where {T<:Real} = zero(T)
transmission(::ParaxialInterface{T}) where {T<:Real} = one(T)

######################################################################

"""
    FresnelInterface{T} <: OpticalInterface{T}

Interface between two materials with behavior defined according to the [Fresnel equations](https://en.wikipedia.org/wiki/Fresnel_equations), with a specified reflectance and transmission.
Assumes unpolarized light.

```julia
FresnelInterface{T}(insidematerial, outsidematerial; reflectance = 0, transmission = 1, interfacemode = ReflectOrTransmit)
```

The interfacemode can be used to trace rays deterministically. Valid values are defined in the InterfaceMode enum.
Reflect means that all values are reflected, Transmit means that all values are transmitted. ReflectOrTransmit will randomly
reflect and transmit rays with the distribution given by the reflection and transmission arguments. This is also the default.
In all cases the power recorded with the ray is correctly updated. This can be used to fake sequential raytracing. For
example a beamsplitter surface may be set to either Reflect or Transmit to switch between the two outgoing ray paths.

"""
struct FresnelInterface{T} <: OpticalInterface{T}
    # storing glasses as IDs (integer) rather than the whole thing seems to improve performance significantly, even when the Glass type is a fixed size (i.e. the interface is pointer-free)
    insidematerial::GlassID
    outsidematerial::GlassID
    reflectance::T
    transmission::T
    interfacemode::InterfaceMode

    function FresnelInterface{T}(insidematerial::Z, outsidematerial::Y; reflectance::T = zero(T), transmission::T = one(T), interfacemode = ReflectOrTransmit) where {Z<:OpticSim.GlassCat.AbstractGlass,Y<:OpticSim.GlassCat.AbstractGlass,T<:Real}
        return FresnelInterface{T}(glassid(insidematerial), glassid(outsidematerial), reflectance = reflectance, transmission = transmission, interfacemode = interfacemode)
    end

    function FresnelInterface{T}(insidematerialid::GlassID, outsidematerialid::GlassID; reflectance::T = zero(T), transmission::T = one(T), interfacemode = ReflectOrTransmit) where {T<:Real}
        @assert zero(T) <= reflectance <= one(T)
        @assert zero(T) <= transmission <= one(T)
        @assert reflectance + transmission <= one(T)
        return new{T}(insidematerialid, outsidematerialid, reflectance, transmission, interfacemode)
    end
end
export FresnelInterface

function Base.show(io::IO, a::FresnelInterface{R}) where {R<:Real}
    print(io, "FresnelInterface($(glassname(a.insidematerial)), $(glassname(a.outsidematerial)), $(a.reflectance), $(a.transmission), $(a.interfacemode))")
end

insidematerialid(a::FresnelInterface{T}) where {T<:Real} = a.insidematerial
outsidematerialid(a::FresnelInterface{T}) where {T<:Real} = a.outsidematerial
reflectance(a::FresnelInterface{T}) where {T<:Real} = a.reflectance
transmission(a::FresnelInterface{T}) where {T<:Real} = a.transmission
interfacemode(a::FresnelInterface{T}) where {T<:Real} = a.interfacemode

transmissiveinterface(::Type{T}, insidematerial::X, outsidematerial::Y) where {T<:Real,X<:OpticSim.GlassCat.AbstractGlass,Y<:OpticSim.GlassCat.AbstractGlass} = FresnelInterface{T}(insidematerial, outsidematerial, reflectance = zero(T), transmission = one(T))
reflectiveinterface(::Type{T}, insidematerial::X, outsidematerial::Y) where {T<:Real,X<:OpticSim.GlassCat.AbstractGlass,Y<:OpticSim.GlassCat.AbstractGlass} = FresnelInterface{T}(insidematerial, outsidematerial, reflectance = one(T), transmission = zero(T))
opaqueinterface(::Type{T} = Float64) where {T<:Real} = FresnelInterface{T}(OpticSim.GlassCat.Air, OpticSim.GlassCat.Air, reflectance = zero(T), transmission = zero(T))
export opaqueinterface

######################################################################

"""
    ThinGratingInterface{T} <: OpticalInterface{T}

Interface representing an idealized thin grating. `period` is in microns, `vector` should lie in the plane of the surface.
Transmission and reflectance can be specified for an arbitrary number of orders up to 10, selected using the `maxorder` and `minorder` parameters.
If `nothing` then `reflectance` is assumed to be **0** and `transmission` is assumed to be **1**.

```julia
ThinGratingInterface(vector, period, insidematerial, outsidematerial; maxorder = 1, minorder = -1, reflectance = nothing, transmission = nothing)
```
"""
struct ThinGratingInterface{T} <: OpticalInterface{T}
    insidematerial::GlassID
    outsidematerial::GlassID
    vector::SVector{3,T}
    period::T
    maxorder::Int
    minorder::Int
    transmission::SVector{GRATING_MAX_ORDERS,T}
    reflectance::SVector{GRATING_MAX_ORDERS,T}

    function ThinGratingInterface(vector::SVector{3,T}, period::T, insidematerial::Z, outsidematerial::Y; maxorder::Int = 1, minorder::Int = -1, reflectance::Union{Nothing,AbstractVector{T}} = nothing, transmission::Union{Nothing,AbstractVector{T}} = nothing) where {T<:Real,Z<:OpticSim.GlassCat.AbstractGlass,Y<:OpticSim.GlassCat.AbstractGlass}
        @assert maxorder >= minorder
        norders = maxorder - minorder + 1
        @assert norders <= GRATING_MAX_ORDERS "Thin grating is limited to $GRATING_MAX_ORDERS orders"
        @assert zero(T) <= sum(reflectance) <= one(T)
        @assert zero(T) <= sum(transmission) <= one(T)
        @assert zero(T) <= sum(reflectance) + sum(transmission) <= one(T)
        @assert ((reflectance === nothing) || length(reflectance) == norders) && ((transmission === nothing) || length(transmission) == norders)
        if reflectance !== nothing
            sreflectance = vcat(SVector{length(reflectance),T}(reflectance), ones(SVector{GRATING_MAX_ORDERS - length(reflectance),T}))
        else
            sreflectance = zeros(SVector{GRATING_MAX_ORDERS,T})
        end
        if transmission !== nothing
            stransmission = vcat(SVector{length(transmission),T}(transmission), ones(SVector{GRATING_MAX_ORDERS - length(transmission),T}))
        else
            stransmission = ones(SVector{GRATING_MAX_ORDERS,T})
        end
        new{T}(glassid(insidematerial), glassid(outsidematerial), normalize(vector), period, maxorder, minorder, sreflectance, stransmission)
    end
end
export ThinGratingInterface

insidematerialid(a::ThinGratingInterface{T}) where {T<:Real} = a.insidematerial
outsidematerialid(a::ThinGratingInterface{T}) where {T<:Real} = a.outsidematerial
reflectance(a::ThinGratingInterface{T}, order::Int) where {T<:Real} = a.reflectance[order - a.minorder + 1]
transmission(a::ThinGratingInterface{T}, order::Int) where {T<:Real} = a.transmission[order - a.minorder + 1]

function Base.show(io::IO, a::ThinGratingInterface{R}) where {R<:Real}
    print(io, "ThinGratingInterface($(glassname(a.insidematerial)), $(glassname(a.outsidematerial)), $(a.vector), $(a.period), $(a.transmission), $(a.reflectance))")
end

######################################################################

"""
`ConvergingBeam`, `DivergingBeam` or `CollimatedBeam`, defines the behavior of a beam in a [`HologramInterface`](@ref).
"""
@enum BeamState ConvergingBeam DivergingBeam CollimatedBeam
export BeamState, ConvergingBeam, DivergingBeam, CollimatedBeam

"""
    HologramInterface{T} <: OpticalInterface{T}

Interface representing a _thick_ hologram (though geometrically thin).
The efficiency, `η`, is calculated using Kogelnik's coupled wave theory so is only valid for the first order.
If the zero order is included then it has efficiency `1 - η`.
Also assumes that the HOE was recorded under similar conditions to the playback conditions, `thickness` is in microns.

`BeatState` arguments can be one of `ConvergingBeam`, `DivergingBeam` and `CollimatedBeam`. In the first two cases `signalpointordir` and `referencepointordir` are 3D point in global coordinate space. For `CollimatedBeam` they are normalized direction vectors.

For reference, see:
- _Coupled Wave Theory for Thick Hologram Gratings_ - H Kogelnik, 1995
- _Sequential and non-sequential simulation of volume holographic gratings_ - M Kick et al, 2018

```julia
HologramInterface(signalpointordir::SVector{3,T}, signalbeamstate::BeamState, referencepointordir::SVector{3,T}, referencebeamstate::BeamState, recordingλ::T, thickness::T, beforematerial, substratematerial, aftermaterial, signalrecordingmaterial, referencerecordingmaterial, RImodulation::T, include0order  = false)
```
"""
struct HologramInterface{T} <: OpticalInterface{T}
    beforematerial::GlassID
    substratematerial::GlassID
    aftermaterial::GlassID
    signalpointordir::SVector{3,T}
    signalbeamstate::BeamState
    referencepointordir::SVector{3,T}
    referencebeamstate::BeamState
    recordingλ::T
    signalrecordingmaterial::GlassID
    referencerecordingmaterial::GlassID
    thickness::T
    RImodulation::T
    include0order::Bool

    function HologramInterface(signalpointordir::SVector{3,T}, signalbeamstate::BeamState, referencepointordir::SVector{3,T}, referencebeamstate::BeamState, recordingλ::T, thickness::T, beforematerial::Z, substratematerial::X, aftermaterial::Y, signalrecordingmaterial::P, referencerecordingmaterial::Q, RImodulation::T, include0order::Bool = false) where {T<:Real,X<:OpticSim.GlassCat.AbstractGlass,Z<:OpticSim.GlassCat.AbstractGlass,Y<:OpticSim.GlassCat.AbstractGlass,P<:OpticSim.GlassCat.AbstractGlass,Q<:OpticSim.GlassCat.AbstractGlass}
        @assert substratematerial !== OpticSim.GlassCat.Air "Substrate can't be air"
        if signalbeamstate === CollimatedBeam
            signalpointordir = normalize(signalpointordir)
        end
        if referencebeamstate === CollimatedBeam
            referencepointordir = normalize(referencepointordir)
        end
        new{T}(glassid(beforematerial), glassid(substratematerial), glassid(aftermaterial), signalpointordir, signalbeamstate, referencepointordir, referencebeamstate, recordingλ, glassid(signalrecordingmaterial), glassid(referencerecordingmaterial), thickness, RImodulation, include0order)
    end

    # 'zero' constructor to use for blank elements in SVector - not ideal...
    HologramInterface{T}() where {T<:Real} = new{T}(NullInterface(T), 0, NullInterface(T), zeros(SVector{3,T}), CollimatedBeam, zeros(SVector{3,T}), CollimatedBeam, zero(T), 0, 0, zero(T), zero(T), false)
end
export HologramInterface

insidematerialid(a::HologramInterface{T}) where {T<:Real} = a.beforematerial
outsidematerialid(a::HologramInterface{T}) where {T<:Real} = a.aftermaterial
substratematerialid(a::HologramInterface{T}) where {T<:Real} = a.substratematerial

function Base.show(io::IO, a::HologramInterface{R}) where {R<:Real}
    print(io, "HologramInterface($(glassname(a.beforematerial)), $(glassname(a.substratematerial)), $(glassname(a.aftermaterial)), $(a.signalpointordir), $(a.signalbeamstate), $(a.referencepointordir), $(a.referencebeamstate), $(a.recordingλ), $(a.thickness), $(a.RImodulation))")
end

"""
    MultiHologramInterface{T} <: OpticalInterface{T}

Interface to represent multiple overlapped [`HologramInterface`](@ref)s on a single surface. Each ray randomly selects an interface to use.

```julia
MultiHologramInterface(interfaces::Vararg{HologramInterface{T}})
MultiHologramInterface(interfaces::Vector{HologramInterface{T}})
```
"""
struct MultiHologramInterface{T} <: OpticalInterface{T}
    interfaces::Ptr{HologramInterface{T}} # pointer to the array of hologram interfaces TODO!! fix this to not be hacky...
    numinterfaces::Int

    MultiHologramInterface(interfaces::Vararg{HologramInterface{T},N}) where {T<:Real,N} = MultiHologramInterface(collect(interfaces))

    function MultiHologramInterface(interfaces::Vector{HologramInterface{T}}) where {T<:Real}
        N = length(interfaces)
        @assert N > 1 "Don't need to use MultiHologramInterface if only 1 hologram"
        p = pointer(interfaces)
        push!(MultiHologramInterfaceRefCache, interfaces)
        new{T}(p, N)
    end
end
export MultiHologramInterface
interface(m::MultiHologramInterface{T}, i::Int) where {T<:Real} = (@assert i <= m.numinterfaces; return unsafe_load(m.interfaces, i))

# need this to keep julia references to the arrays used in MultiHologramInterface so they don't get garbage collected
const MultiHologramInterfaceRefCache = []

######################################################################

# a few things need a concrete type to prevent allocations resulting from the ambiguities introduce by an abstract type (i.e. OpticalInterface{T})
# if adding new subtypes of OpticalInterface they must be added to this definition as well (and only paramaterised by T so they are fully specified)
AllOpticalInterfaces{T} = Union{NullInterface{T},ParaxialInterface{T},FresnelInterface{T},HologramInterface{T},MultiHologramInterface{T},ThinGratingInterface{T}}
NullOrFresnel{T} = Union{NullInterface{T},FresnelInterface{T}}

# can't have more than 4 types here or we get allocations because inference on the union stops: https://github.com/JuliaLang/julia/blob/1e6e65691254a7fe81f5da8706bb30aa6cb3f8d2/base/compiler/types.jl#L113
# use the macro to get around this: @unionsplit OpticalInterface T x f(x)
import InteractiveUtils
# TODO REMOVE
macro unionsplit(type, paramtype, var, call)
    st = InteractiveUtils.subtypes(@eval $type)
    orig_expr = Expr(:if, :($var isa $(st[1]){$paramtype}), call)
    expr = orig_expr
    for t in st[2:end]
        push!(expr.args, Expr(:elseif, :($var isa $t{$paramtype}), call))
        expr = expr.args[end]
    end
    push!(expr.args, :(error("unhandled type")))
    return :($(esc(orig_expr)))
end
