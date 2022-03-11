# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

export drawrectlattice, drawhexneighbors, drawhexregion, drawhexrect, drawhexrectcolors
export hex3cluster, hex3RGB, hexRGBW, hex12RGB
export drawhex3RGB, drawhex12RGB

drawrectlattice() = Vis.drawcells(Repeat.RectangularBasis(),50.0,SMatrix{2,4,Int64}(0,0,0,1,1,0,1,1))

"""draw the 2 ring neighbors of the hex cell at coordinates (0,0)"""
drawhexneighbors() =  Vis.drawcells(Repeat.HexBasis1(),50,Repeat.neighbors(Repeat.HexBasis1,(0,0),2))

"""draw hex cell at coordinates (0,0) and the 1 and 2 ring neighbors"""
drawhexregion() = Vis.drawcells(Repeat.HexBasis1(),50,Repeat.region(Repeat.HexBasis1,(0,0),2))

"""draw hex cells that fit within a rectangular box centered at coordinates (0,0). Use fill color yellow."""
function drawhexrect() 
    cells = Repeat.hexcellsinbox(2,2)
    Vis.drawcells(Repeat.HexBasis1(),50,cells,color = repeat(["yellow"],length(cells)))
end

"""draw hex cells that fit within a rectangular box centered at coordinates (0,0). Use random fill colors selected for maximum distinguishability."""
function drawhexrectcolors()
    cells = Repeat.hexcellsinbox(4,4)
    Vis.drawcells(Repeat.HexBasis1(),30,cells)
end

""" Create a LatticeCluser with three elements at (0,0),(-1,0),(-1,1) coordinates in the HexBasis1 lattice"""
function hex3cluster()
    clusterelts = SVector((0,0),(-1,0),(-1,1))
    eltlattice = HexBasis1()
    clusterbasis = LatticeBasis(( -1,2),(2,-1))
    return LatticeCluster(clusterbasis,eltlattice,clusterelts)
end


""" Create a ClusterWithProperties with four types of elements, R,G,B,W """
function hexRGBW()
    clusterelements = SVector((0,0),(-1,0),(-1,1),(0,-1))
    colors = [colorant"red",colorant"green",colorant"blue",colorant"white"]
    names = ["R","G","B","W"]
    eltlattice = HexBasis1()
    clusterbasis = LatticeBasis((0,2),(2,-2))
    lattice = LatticeCluster(clusterbasis,eltlattice,clusterelements)
    properties =  DataFrame(Color = colors, Name = names)
    return ClusterWithProperties(lattice,properties)
end


""" draw 3 repeats of hex3RGB cluster """
drawhex3RGB() = Vis.draw(hex3RGB(),[0 1 0; 0 0 1])

""" draw 3 repeats of hex12RGB cluster """
drawhex12RGB() = Vis.draw(Repeat.Multilens.hex12RGB(),[0 1 0 1; 0 0 1 1])
