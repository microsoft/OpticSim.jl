# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.



#############################################################################

drawneighbors() = Vis.@wrapluxor Vis.drawhexcells(50,Repeat.neighbors(Repeat.HexBasis1,(0,0),2))
export drawneighbors

drawhexregion() = Vis.@wrapluxor Vis.drawhexcells(50,Repeat.hexregion((0,0),2))
export drawhexregion

drawhexrect(color) = Vis.@wrapluxor :(Vis.drawhexcells(50,Repeat.hexcellsinbox(4,4),$color)
export drawhexrect

drawhexrectcolors() = Vis.@wrapluxor Vis.drawhexcells(50,Repeat.hexcellsinbox(4,4))
export drawhexrectcolors