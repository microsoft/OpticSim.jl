# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

abstract type AbstractRayGenerator{T<:Real} end

"""
    GeometricRayGenerator{T,O<:RayOriginGenerator{T}} <: AbstractRayGenerator{T}

Generates geometric [`Ray`](@ref)s according to the specific implementation of the subclass.
"""
abstract type GeometricRayGenerator{T} <: AbstractRayGenerator{T} end

"""
    OpticalRayGenerator{T} <: AbstractRayGenerator{T}

Generates [`OpticalRay`](@ref)s according to the specific implementation of the subclass.
"""
abstract type OpticalRayGenerator{T} <: AbstractRayGenerator{T} end
export AbstractRayGenerator, GeometricRayGenerator, OpticalRayGenerator
