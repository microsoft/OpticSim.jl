# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

module GlassCat

using Polynomials
using Plots
using StringEncodings
using Unitful
using StaticArrays
using Base: @.
import Unitful: Length, Temperature, Quantity, Units

include("constants.jl")

include("GlassTypes.jl")
export GlassID, info, glassid, glassname, glassforid
include("Air.jl")
export Air, isair

# include built glass cat source files
include("data/jl/AGFGlassCat.jl")
include("data/jl/CARGILLE.jl")

# include functionality for managing runtime (dynamic) glass cats: MIL_GLASSES and MODEL_GLASSES
include("runtime.jl")
export glassfromMIL, modelglass

# include functions for searching the glass cats
include("search.jl")
export glasscatalogs, glassnames, findglass

include("utilities.jl")
export plot_indices, index, polyfit_indices, absairindex, absorption

# include utility functions for maintaining the AGF source list
include("sources.jl")
export add_agf

end # module
export GlassCat
