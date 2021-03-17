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
struct OpticalRay{T,N} <: AbstractRay{T,N}
    ray::Ray{T,N}
    power::T
    wavelength::T
    opl::T
    nhits::Int
    sourcepower::T
    sourcenum::Int

    function OpticalRay(ray::Ray{T,N}, power::T, wavelength::T; opl::T = zero(T), nhits::Int = 0, sourcenum::Int = 0, sourcepower::T = power) where {T<:Real,N}
        return new{T,N}(ray, power, wavelength, opl, nhits, sourcepower, sourcenum)
    end

    function OpticalRay(origin::SVector{N,T}, direction::SVector{N,T}, power::T, wavelength::T; opl::T = zero(T), nhits::Int = 0, sourcenum::Int = 0, sourcepower::T = power) where {T<:Real,N}
        return new{T,N}(Ray(origin, normalize(direction)), power, wavelength, opl, nhits, sourcepower, sourcenum)
    end

    function OpticalRay(origin::AbstractArray{T,1}, direction::AbstractArray{T,1}, power::T, wavelength::T; opl::T = zero(T), nhits::Int = 0, sourcenum::Int = 0, sourcepower::T = power) where {T<:Real}
        @assert length(origin) == length(direction)
        N = length(origin)
        return new{T,N}(Ray(SVector{N,T}(origin), normalize(SVector{N,T}(direction))), power, wavelength, opl, nhits, sourcepower, sourcenum)
    end

    # Convenience constructor. Not as much typing
    OpticalRay(ox::T, oy::T, oz::T, dx::T, dy::T, dz::T; wavelength = 0.55) where {T<:Real} = OpticalRay(SVector{3,T}(ox, oy, oz), SVector{3,T}(dx, dy, dz), one(T), T(wavelength), one(T)) #doesn't have to be inside struct definition but if it is then VSCode displays hover information. If it's outside the struct definition it doesn't.
end
export OpticalRay

ray(r::OpticalRay{T,N}) where {T<:Real,N} = r.ray
direction(r::OpticalRay{T,N}) where {T<:Real,N} = direction(ray(r))
origin(r::OpticalRay{T,N}) where {T<:Real,N} = origin(ray(r))
power(r::OpticalRay{T,N}) where {T<:Real,N} = r.power
wavelength(r::OpticalRay{T,N}) where {T<:Real,N} = r.wavelength
pathlength(r::OpticalRay{T,N}) where {T<:Real,N} = r.opl
nhits(r::OpticalRay{T,N}) where {T<:Real,N} = r.nhits
sourcepower(r::OpticalRay{T,N}) where {T<:Real,N} = r.sourcepower
sourcenum(r::OpticalRay{T,N}) where {T<:Real,N} = r.sourcenum
export ray, power, wavelength, pathlength, nhits, sourcepower, sourcenum

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

function Base.:*(a::RigidBodyTransform{T}, r::OpticalRay{T,N}) where {T,N}
    return OpticalRay(a * ray(r), power(r), wavelength(r), opl = pathlength(r), nhits = nhits(r), sourcenum = sourcenum(r), sourcepower = sourcepower(r))
end
