# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.


struct RectangularBasis{N,T} <: AbstractBasis{N,T}
    RectangularBasis(::Type{T} = Float64) where{T<:Real} = new{2,T}()
end
export RectangularBasis

basismatrix(::RectangularBasis{2,T}) where{T} = SMatrix{2,2,T}(T(1),T(0),T(0),T(1))

# SVector{2,SVector{2,T}}(hexe₁(T),hexe₂(T))

"""Returns the vertices of the unit tile polygon for the basis"""
function vertices(::RectangularBasis{2,T}) where{T}
	return SMatrix{2,4,T}(
			-.5, -.5,
			-.5, .5,
			.5, .5,
			.5, -.5)
end

rectangularlattice(ipitch::T = 1.0,jpitch::T = 1.0) where{T<:Real} = LatticeBasis(SMatrix{2,2,T}(
    ipitch, 0,
    0, jpitch))
export rectangularlattice