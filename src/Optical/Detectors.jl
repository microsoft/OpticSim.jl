# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

import Base

"""
    AbstractDetector{T<:Number}

AbstractDetector is the base class for detectors. Detectors are attached to surfaces
and define how which information from impinging rays are stored and how they are stored.

The type parameter T defines the number format of the stored detector data.

Subtypes of AbstractDetector must define the update! and name methods.
"""
abstract type AbstractDetector{T<:Number} end

"""
    update!(d::AbstractDetector, s::Surface, t::LensTrace)
    
Update detector d attached to surface s with information from the raytrace t.
"""
function update! end

"""
    name(d::AbstractDetector)

Return name of detector.
"""
function name end

#
# HierarchicalImage detector
#

"""
    ImageDetector{T<:Number} 
    
ImageDetector is a detector that stores a HierarchicalImage.
"""
struct ImageDetector{T} <: AbstractDetector{T}
    name::String
    image::HierarchicalImage{T}
    function ImageDetector(name, detectorpixelsx::Int = 1000, detectorpixelsy::Int = 1000, D = Float32)
        new{D}(name, HierarchicalImage{D}(detectorpixelsy, detectorpixelsx))
    end
end

function update!(d::ImageDetector{T}, s, t) where T
    pixu, pixv = uvtopix(s, uv(intersection(t)), size(d))
    d.image[pixv, pixu] += T(sourcepower(ray(t)))
end

function name(d::ImageDetector)
    d.name
end

function reset!(d::ImageDetector{T}) where T
    reset!(d.image)
end

function Base.size(d::ImageDetector{T}) where T
    size(d.image)
end

function Base.getindex(d::ImageDetector{T}, a, b) where T
    getindex(d.image, a, b)
end

function Base.iterate(d::ImageDetector{T}, args...) where T
    Base.iterate(d.image, args...)
end
