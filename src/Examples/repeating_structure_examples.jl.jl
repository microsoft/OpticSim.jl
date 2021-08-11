# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.



#############################################################################

"""draw the 2 ring neighbors of the hex cell at coordinates (0,0)"""
drawhexneighbors() =  Luxor.@svg Vis.drawcells(Repeat.HexBasis1(),50,Repeat.neighbors(Repeat.HexBasis1,(0,0),2))
export drawneighbors

"""draw hex cell at coordinates (0,0) and the 1 and 2 ring neighbors"""
drawhexregion() = Luxor.@svg Vis.drawcells(Repeat.HexBasis1(),50,Repeat.region(Repeat.HexBasis1,(0,0),2))
export drawhexregion

"""draw hex cells that fit within a rectangular box centered at coordinates (0,0). Use fill color yellow."""
function drawhexrect() 
    cells = Repeat.hexcellsinbox(2,2)
    Luxor.@svg Vis.drawcells(Repeat.HexBasis1(),50,cells,color = repeat(["yellow"],length(cells)))
end
export drawhexrect

"""draw hex cells that fit within a rectangular box centered at coordinates (0,0). Use random fill colors selected for maximum distinguishability."""
function drawhexrectcolors()
    cells = Repeat.hexcellsinbox(4,4)
    Luxor.@svg Vis.drawcells(Repeat.HexBasis1(),30,cells)
end
export drawhexrectcolors