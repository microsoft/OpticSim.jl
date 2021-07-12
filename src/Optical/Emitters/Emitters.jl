# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

#=
Emitters are defined by Pixels and SpatialLayouts. An emitter has a spectrum, and an optical power distribution over the hemisphere. 
These are intrinsic physical properties of the emitter.
=#

module Emitters
export generate, visual_size, apply
export right, up, forward
export Spectrum, Directions, Origins, AngularPower, Sources
export pointemitter, collimatedemitter

using ...OpticSim, ..Geometry
using LinearAlgebra
import Unitful.Length
using Unitful.DefaultSymbols

# defining name placeholders to override in nested modules

"""
    generate(???)

[TODO]
"""
function generate end

"""
    visual_size(???)

[TODO]
"""
function visual_size end

"""
    apply(???)

[TODO] Returns ray power
"""
function apply end

include("Spectrum.jl")
include("Directions.jl")
include("Origins.jl")
include("AngularPower.jl")
include("Sources.jl")
include("StandardEmitters.jl")

end # module Emitters

export Emitters
