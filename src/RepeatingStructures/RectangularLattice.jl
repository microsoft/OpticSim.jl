# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

rectangularlattice(ipitch::T = 1.0,jpitch::T = 1.0) where{T<:Real} = LatticeBasis(ipitch*SVector{2,T}(1,0),jpitch*SVector{2,T}(0,1))
export rectangularlattice