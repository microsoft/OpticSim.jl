# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

#Various sizes of hexagonal clusters that are useful in lenslet design. To see what they look like use the functions in HexTilings to make lattices with color properties and then draw them with Vis.draw
#example: Vis.draw(hex3RGB())

function hex3()
    clusterelements = SVector((0, 0), (-1, 0), (-1, 1))
    eltlattice = HexBasis1()
    clusterbasis = LatticeBasis((-1, 2), (2, -1))
    return LatticeCluster(clusterbasis, eltlattice, clusterelements)

end
export hex3

function hex4()
    clusterelements = SVector((-1, 0), (-1, 1), (0, 0), (0, 1))
    eltlattice = HexBasis1()
    clusterbasis = LatticeBasis((2, 0), (0, 2))
    return LatticeCluster(clusterbasis, eltlattice, clusterelements)
end
export hex4

function hex7()
    clusterelements = SVector((1, -1), (0, -1), (-1, 0), (0, 0), (-1, 1), (0, 1), (1, 0), )
    eltlattice = HexBasis1()
    clusterbasis = LatticeBasis((2, 2), (-3, 2))
    return LatticeCluster(clusterbasis, eltlattice, clusterelements)
end
export hex7


function hex9()
    clusterelements = SVector((-1, -1), (-1, 0), (-2, 1), (0, -1), (0, 0), (-1, 1), (1, -1), (1, 0), (0, 1))
    eltlattice = HexBasis1()
    clusterbasis = LatticeBasis((3, 0), (0, 3))
    return LatticeCluster(clusterbasis, eltlattice, clusterelements)
end
export hex9

function hex12()
    clusterelements = SVector(
        (-1, 1),(0, 1),
        (-1, 0),(0, 0),(1, 0),
        (-1, -1),(0, -1),(1, -1),(2, -1),
        (0, -2),(1, -2),(2, -2)
    )

    eltlattice = HexBasis3()
    clusterbasis = LatticeBasis((2, 2), (-3, 2))
    return LatticeCluster(clusterbasis, eltlattice, clusterelements)
end
export hex12

function hex19()
    eltlattice = HexBasis1()
    clusterbasis = LatticeBasis((5, 0), (-2, 4))
    return LatticeCluster(clusterbasis, eltlattice, region(HexBasis1,(0,0),2))
end
export hex19

function hexn(regionsize::Int64)
    eltlattice = HexBasis1()
    clusterbasis = LatticeBasis((2*regionsize+1, -regionsize), (regionsize, regionsize+1))
    return LatticeCluster(clusterbasis, eltlattice, region(HexBasis1,(0,0),regionsize))

end
export hexn

