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

""" Create a LatticeCluser with three elements at (0,0),(-1,0),(-1,1) coordinates in the HexBasis1 lattice"""
function hex3cluster()
    clusterelts = SVector((0,0),(-1,0),(-1,1))
    eltlattice = HexBasis1()
    clusterbasis = LatticeBasis(( -1,2),(2,-1))
    return LatticeCluster(clusterbasis,eltlattice,clusterelts)
end
export hex3cluster

""" Create a ClusterWithProperties with three types of elements, R,G,B """
function hex3RGB()
    clusterelements = SVector((0,0),(-1,0),(-1,1))
    colors = [color("red"),color("green"),color("blue")]
    names = ["R","G","B"]
    eltlattice = HexBasis1()
    clusterbasis = LatticeBasis(( -1,2),(2,-1))
    lattice = LatticeCluster(clusterbasis,eltlattice,clusterelements)
    properties =  DataFrame(Color = colors, Name = names)
    return ClusterWithProperties(lattice,properties)
end
export hex3RGB

""" Create a ClusterWithProperties with four types of elements, R,G,B,W """
function hexRGBW()
    clusterelements = SVector((0,0),(-1,0),(-1,1),(0,-1))
    colors = [color("red"),color("green"),color("blue"),color("white")]
    names = ["R","G","B","W"]
    eltlattice = HexBasis1()
    clusterbasis = LatticeBasis((0,2),(2,-2))
    lattice = LatticeCluster(clusterbasis,eltlattice,clusterelements)
    properties =  DataFrame(Color = colors, Name = names)
    return ClusterWithProperties(lattice,properties)
end
export hexRGBW

function hex12RGB()
    clusterelements = SVector(
        (-1,1),(0,1),
        (-1,0),(0,0),(1,0),
        (-1,-1),(0,-1),(1,-1),(2,-1),
        (0,-2),(1,-2),(2,-2)
    )
    red = color("red")
    grn = color("green")
    blu = color("blue")
    colors = [
        grn,blu,
        blu,red,grn,
        red,grn,blu,red,
        blu,red,grn
        ]
    names = [
        "G3","B0",
        "B2","R1","G2",
        "R0","G1","B1","R3",
        "B3","R2","G0"
    ]
    eltlattice = HexBasis3()
    clusterbasis = LatticeBasis((2,2),(-3,2))
    lattice = LatticeCluster(clusterbasis,eltlattice,clusterelements)
    properties =  DataFrame(Color = colors, Name = names)
    return ClusterWithProperties(lattice,properties)
end
export hex12RGB

""" draw 3 repeats of hex3RGB cluster """
drawhex3RGB() = Vis.draw(hex3RGB(),[0 1 0; 0 0 1])
export drawhex3RGB

""" draw 3 repeats of hex12RGB cluster """
drawhex12RGB() = Vis.draw(hex12RGB(),[0 1 0; 0 0 1])
export drawhex12RGB



