# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.



#############################################################################

"""draw the 2 neighbors of the hex cell at coordinates (0,0)"""
drawhexneighbors() =  Luxor.@svg Vis.drawcells(Repeat.HexBasis1(),50,Repeat.neighbors(Repeat.HexBasis1,(0,0),2))
export drawneighbors

"""draw hex cell at coordinates (0,0) and the 1 and 2 neighbors"""
drawhexregion() = Luxor.@svg Vis.drawcells(Repeat.HexBasis1(),50,Repeat.region(Repeat.HexBasis1,(0,0),2))
export drawhexregion

"""draw hex cells that fit within a rectangular box centered at coordinates (0,0). Use fill color yellow."""
drawhexrect() = Luxor.@svg Vis.drawcells(Repeat.HexBasis1(),50,Repeat.hexcellsinbox(2,2),"yellow")
export drawhexrect

"""draw hex cells that fit within a rectangular box centered at coordinates (0,0). Use random fill colors selected for maximum distinguishability."""
drawhexrectcolors() =  Luxor.@svg Vis.drawcells(Repeat.HexBasis1(),50,Repeat.hexcellsinbox(4,4))
export drawhexrectcolors