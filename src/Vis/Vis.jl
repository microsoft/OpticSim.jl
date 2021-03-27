export Vis
module Vis

import Makie
import Makie.AbstractPlotting
import Makie.AbstractPlotting.MakieLayout

# If using precompiled system image (which we always are) you have to run AbstractPlotting.__init__() after loading Makie
# during precompilation, the display stack gets shuffled around such that the Makie display does not take priority.
# See https://discourse.julialang.org/t/makie-doesnt-display-plot-when-using-a-custom-julia-sysimage/38515.
__init__() = AbstractPlotting.__init__()

include("Visualization.jl")
include("Emitters.jl")
end # module Vis
