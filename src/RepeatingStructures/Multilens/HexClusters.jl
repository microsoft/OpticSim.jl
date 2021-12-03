# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

#Various sizes of hexagonal clusters that are useful in lenslet design. To see what they look like use the functions in HexTilings to make lattices with color properties and then draw them with Vis.draw
#example: Vis.draw(hex3RGB())

function hex3()
    clusterelements = SVector((0, 0), (-1, 0), (-1, 1))
    eltlattice = Repeat.HexBasis1()
    clusterbasis = Repeat.LatticeBasis((-1, 2), (2, -1))
    return Repeat.LatticeCluster(clusterbasis, eltlattice, clusterelements)

end
export hex3

function hex4()
    clusterelements = SVector((-1, 0), (-1, 1), (0, 0), (0, 1))
    eltlattice = Repeat.HexBasis1()
    clusterbasis = Repeat.LatticeBasis((2, 0), (0, 2))
    return Repeat.LatticeCluster(clusterbasis, eltlattice, clusterelements)
end
export hex4

function hex7()
    clusterelements = SVector((1, -1), (0, -1), (-1, 0), (0, 0), (-1, 1), (0, 1), (1, 0), )
    eltlattice = Repeat.HexBasis1()
    clusterbasis = Repeat.LatticeBasis((3, 0), (0, 3))
    clusterbasis = Repeat.LatticeBasis((2, 2), (-3, 2))
    return Repeat.LatticeCluster(clusterbasis, eltlattice, clusterelements)
end
export hex7


function hex9()
    clusterelements = SVector((-1, -1), (-1, 0), (-2, 1), (0, -1), (0, 0), (-1, 1), (1, -1), (1, 0), (0, 1))
    eltlattice = Repeat.HexBasis1()
    clusterbasis = Repeat.LatticeBasis((3, 0), (0, 3))
    return Repeat.LatticeCluster(clusterbasis, eltlattice, clusterelements)
end
export hex9

function hex12()
    clusterelements = SVector(
        (-1, 1),(0, 1),
        (-1, 0),(0, 0),(1, 0),
        (-1, -1),(0, -1),(1, -1),(2, -1),
        (0, -2),(1, -2),(2, -2)
    )

    eltlattice = Repeat.HexBasis3()
    clusterbasis = Repeat.LatticeBasis((2, 2), (-3, 2))
    return Repeat.LatticeCluster(clusterbasis, eltlattice, clusterelements)
end
export hex12

function hex19()
    clusterelements = SVector(
        (0, 0), # ring 0
        (-1, 0),(0, -1),(1, -1),(1, 0),(0, 1),(-1, 1), # ring1
        (-2, 0),(-1, -1),(0, -2),(1, -2),(2, -2),(2, -1),(2, 0),(1, 1),(0, 2),(-1, 2),(-2, 2),(-2, 1) # ring2
        )
    
    eltlattice = Repeat.HexBasis3()
    clusterbasis = Repeat.LatticeBasis((5, 0), (-2, 4))
    return Repeat.LatticeCluster(clusterbasis, eltlattice, clusterelements)
end
export hex19

# """ Returns a lattice cluster composed of the lattice elements in rings 0..n. These clusters have regular hexagonal shapes"""
# function regularhexcluster(n::Int)
#     clusterelements = Repeat.region(HexBasis1,(0,0),n)
#     eltlattice = Repeat.HexBasis1()
#     clusterbasis = Repeat.LatticeBasis((2*n-1,0))
#     return Repeat.LatticClusterA()
