# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.


"""
    OpticalRay{T,N} <: AbstractRay{T,N}

Ray with power, wavelength and optical path length.

**NOTE**: we use monte carlo integration to get accurate results on the detector, this means that all rays essentially hit the detector with power = 1 and some rays are thrown away at any interface to correctly match the reflection/transmission at that interface. For inspection purposes we also track the 'instantaneous' power of the ray in the `power` field of the `OpticalRay`.

```julia
OpticalRay(ray::Ray{T,N}, power::T, wavelength::T, opl=zero(T))
OpticalRay(origin::SVector{N,T}, direction::SVector{N,T}, power::T, wavelength::T, opl=zero(T))
```

Has the following accessor methods:
```julia
ray(r::OpticalRay{T,N}) -> Ray{T,N}
direction(r::OpticalRay{T,N}) -> SVector{N,T}
origin(r::OpticalRay{T,N}) -> SVector{N,T}
power(r::OpticalRay{T,N}) -> T
wavelength(r::OpticalRay{T,N}) -> T
pathlength(r::OpticalRay{T,N}) -> T
sourcepower(r::OpticalRay{T,N}) -> T
nhits(r::OpticalRay{T,N}) -> Int
sourcenum(r::OpticalRay{T,N}) -> Int
```
"""
struct OpticalRay{T,N,P<:Polarization.AbstractPolarization{T}} <: AbstractRay{T,N}
    ray::Ray{T,N}
    power::T
    wavelength::T
    opl::T
    nhits::Int
    sourcepower::T
    sourcenum::Int
    polarization::P #this will only work if N = 3.

    function OpticalRay(ray::Ray{T,N}, power::T, wavelength::T; opl::T = zero(T), nhits::Int = 0, sourcenum::Int = 0, sourcepower::T = power, polarization::P = Polarization.NoPolarization{T}()) where {T<:Real,N,P}
        return new{T,N,P}(ray, power, wavelength, opl, nhits, sourcepower, sourcenum, polarization)
    end

    
    function OpticalRay(origin::SVector{N,T}, direction::SVector{N,T}, power::T, wavelength::T; opl::T = zero(T), nhits::Int = 0, sourcenum::Int = 0, sourcepower::T = power, polarization::P = Polarization.NoPolarization{T}()) where {T<:Real,N,P}
        return new{T,N,P}(Ray(origin, normalize(direction)), power, wavelength, opl, nhits, sourcepower, sourcenum, polarization)
    end

    function OpticalRay(origin::AbstractArray{T,1}, direction::AbstractArray{T,1}, power::T, wavelength::T; opl::T = zero(T), nhits::Int = 0, sourcenum::Int = 0, sourcepower::T = power, polarization::P = Polarization.NoPolarization{T}()) where {T<:Real,P}
        @assert length(origin) == length(direction) "origin (dimension $(length(origin))) and direction (dimension $(length(direction))) vectors do not have the same dimension"
        N = length(origin)
        return new{T,N,P}(Ray(SVector{N,T}(origin), normalize(SVector{N,T}(direction))), power, wavelength, opl, nhits, sourcepower, sourcenum, polarization)
    end

    # Convenience constructor. Not as much typing
    OpticalRay(ox::T, oy::T, oz::T, dx::T, dy::T, dz::T; wavelength = 0.55) where {T<:Real} = OpticalRay(SVector{3,T}(ox, oy, oz), SVector{3,T}(dx, dy, dz), one(T), T(wavelength)) #doesn't have to be inside struct definition but if it is then VSCode displays hover information. If it's outside the struct definition it doesn't.
end
export OpticalRay

ray(r::OpticalRay)  = r.ray
direction(r::OpticalRay)  = direction(ray(r))
origin(r::OpticalRay)  = origin(ray(r))
power(r::OpticalRay)  = r.power
wavelength(r::OpticalRay) = r.wavelength
pathlength(r::OpticalRay)  = r.opl
nhits(r::OpticalRay) = r.nhits
sourcepower(r::OpticalRay)  = r.sourcepower
sourcenum(r::OpticalRay)  = r.sourcenum
polarization(r::OpticalRay) = r.polarization
export ray, power, wavelength, pathlength, nhits, sourcepower, sourcenum,polarization

function Base.print(io::IO, a::OpticalRay{T,N}) where {T,N}
    println(io, "$(rpad("Origin:", 25)) $(origin(a))")
    println(io, "$(rpad("Direction:", 25)) $(direction(a))")
    println(io, "$(rpad("Power:", 25)) $(power(a))")
    println(io, "$(rpad("Source Power:", 25)) $(sourcepower(a))")
    println(io, "$(rpad("Wavelength (in air):", 25)) $(wavelength(a))")
    println(io, "$(rpad("Optical Path Length:", 25)) $(pathlength(a))")
    println(io, "$(rpad("Hits:", 25)) $(nhits(a))")
    if sourcenum(a) != 0
        println(io, "$(rpad("Source Number:", 25)) $(sourcenum(a))")
    end
end

function Base.:*(a::Transform{T}, r::OpticalRay{T,N}) where {T,N}
    return OpticalRay(a * ray(r), power(r), wavelength(r), opl = pathlength(r), nhits = nhits(r), sourcenum = sourcenum(r), sourcepower = sourcepower(r))
end
