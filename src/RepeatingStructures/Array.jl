# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

struct UnitCell
    coordinateframe #this should be computed once when the UnitCell is created. It incorporates the parent coordinateframe transformation and the transformation from the parent to the unit cell.
    emitter
    # lens::LensAssembly
    detectors #have multiple eye positions that get computed all at once. Maybe. Lots of space so might be inefficient.
end

""" coordinateframe is a coordinate frame defined at a point in the repeating structure which all other positions will be measured with respect to. Used to define a globally consistent coordinate frame that pixels in each unit call map to."""
struct RepeatingStructure{N}
    coordinateframe #need to think about how to define this. Maybe using a Transform?
    cells::SVector{N,UnitCell} #may want to access cells not just by vector index, i.e., x,y, or something more complicated for hexagonal systems.
end

""" Maps a pixel position from an emitter in a unit cell into a globally consistent angular coordinate frame with coordinates (θ,ϕ). This allows the contributions from multiple cells to be summed correctly."""
function globalcoordinates(px,py,cell::UnitCell,repeat::RepeatingStructure)::NamedTuple{(:θ,:ϕ),Tuple{Real,Real}}
    #result will be like this (θ = val1, ϕ = val2)
end

function rectangulararray(nrows,ncols)
end

function hexagonalarray()
end


