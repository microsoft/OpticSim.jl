# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

function hex3RGB()
    clusterelements = SVector((0,0),(-1,0),(-1,1))
    colors = [colorant"red",colorant"green",colorant"blue"]
    names = ["R","G","B"]
    eltlattice = Repeat.HexBasis1()
    clusterbasis = Repeat.LatticeBasis(( -1,2),(2,-1))
    lattice = Repeat.LatticeCluster(clusterbasis,eltlattice,clusterelements)
    properties =  DataFrames.DataFrame(Color = colors, Name = names)
    return Repeat.ClusterWithProperties(lattice,properties)
end
export hex3RGB

function hex4RGB()
    clusterelements = SVector((-1,0),(-1,1),(0,0),(0,1))
    r,g,b,w = colorant"red",colorant"green",colorant"blue",colorant"white"
    colors = [r,g,b,w]
    names = ["R","G","B","W"]
    eltlattice = Repeat.HexBasis1()
    clusterbasis = Repeat.LatticeBasis((1,1),(2,-1))
    lattice = Repeat.LatticeCluster(clusterbasis,eltlattice,clusterelements)
    properties =  DataFrames.DataFrame(Color = colors, Name = names)
    return Repeat.ClusterWithProperties(lattice,properties)
end
export hex4RGB

function hex7RGB()
    clusterelements = SVector((1,-1),(0,-1),(-1,0),(0,0),(-1,1),(0,1),(1,0),)
    r,g,b,w = colorant"red",colorant"green",colorant"blue", colorant"white"
    colors = [r,g,b,w,r,g,b]
    names = ["R1","G1","B1","W","R2","G2","B2"]
    eltlattice = Repeat.HexBasis1()
    clusterbasis = Repeat.LatticeBasis((3,0),(0,3))
    lattice = Repeat.LatticeCluster(clusterbasis,eltlattice,clusterelements)
    properties =  DataFrames.DataFrame(Color = colors, Name = names)
    return Repeat.ClusterWithProperties(lattice,properties)
end
export hex7RGB

function hex9RGB()
    clusterelements = SVector((-1,-1),(-1,0),(-2,1),(0,-1),(0,0),(-1,1),(1,-1),(1,0),(0,1))
    r,g,b = colorant"red",colorant"green",colorant"blue"
    colors = [r,g,b,r,g,b,r,g,b]
    names = ["R-1","G-1","B-1","R0","G0","B0","R1","G1","B1"]
    eltlattice = Repeat.HexBasis1()
    clusterbasis = Repeat.LatticeBasis((3,0),(0,3))
    lattice = Repeat.LatticeCluster(clusterbasis,eltlattice,clusterelements)
    properties =  DataFrames.DataFrame(Color = colors, Name = names)
    return Repeat.ClusterWithProperties(lattice,properties)
end
export hex9RGB

function hex12RGB()
    clusterelements = SVector(
        (-1,1),(0,1),
        (-1,0),(0,0),(1,0),
        (-1,-1),(0,-1),(1,-1),(2,-1),
        (0,-2),(1,-2),(2,-2)
    )
    red = colorant"red"
    grn = colorant"green"
    blu = colorant"blue"
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
    eltlattice = Repeat.HexBasis3()
    clusterbasis = Repeat.LatticeBasis((2,2),(-3,2))
    lattice = Repeat.LatticeCluster(clusterbasis,eltlattice,clusterelements)
    properties =  DataFrames.DataFrame(Color = colors, Name = names)
    return Repeat.ClusterWithProperties(lattice,properties)
end
export hex12RGB

function hex19RGB()
    clusterelements = SVector(
        (0,0), #ring 0
        (-1,0),(0,-1),(1,-1),(1,0),(0,1),(-1,1), #ring1
        (-2,0),(-1,-1),(0,-2),(1,-2),(2,-2),(2,-1),(2,0),(1,1),(0,2),(-1,2),(-2,2),(-2,1) #ring2
        )
    red = colorant"red"
    grn = colorant"green"
    blu = colorant"blue"
    wht = colorant"white"
    colors = [
        wht,
        blu,grn,red,blu,grn,red,
        grn,
        blu,red,grn,
        blu,red,grn,
        blu,red,grn,
        blu,red
        ]
    names = [
        "W",
        "B0","G1","R2","B3","G4","R5",
        "G0","B1","R1","G2","B2","R3","G3","B4","R4","G5","B5","R0"
    ]
    eltlattice = Repeat.HexBasis3()
    clusterbasis = Repeat.LatticeBasis((5,0),(-2,4))
    lattice = Repeat.LatticeCluster(clusterbasis,eltlattice,clusterelements)
    properties =  DataFrames.DataFrame(Color = colors, Name = names)
    # properties =  DataFrames.DataFrame(Color = colors)
    return Repeat.ClusterWithProperties(lattice,properties)
end
export hex19RGB