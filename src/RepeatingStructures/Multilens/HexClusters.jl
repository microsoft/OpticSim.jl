# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

function hex3RGB()
    clusterelements = SVector((0, 0), (-1, 0), (-1, 1))
    eltlattice = Repeat.HexBasis1()
    clusterbasis = Repeat.LatticeBasis((-1, 2), (2, -1))
    return Repeat.LatticeCluster(clusterbasis, eltlattice, clusterelements)

end
export hex3RGB

function hex4RGB()
    clusterelements = SVector((-1, 0), (-1, 1), (0, 0), (0, 1))
    eltlattice = Repeat.HexBasis1()
    clusterbasis = Repeat.LatticeBasis((1, 1), (2, -1))
    return Repeat.LatticeCluster(clusterbasis, eltlattice, clusterelements)
end
export hex4RGB

function hex7RGB()
    clusterelements = SVector((1, -1), (0, -1), (-1, 0), (0, 0), (-1, 1), (0, 1), (1, 0), )
    eltlattice = Repeat.HexBasis1()
    clusterbasis = Repeat.LatticeBasis((3, 0), (0, 3))
    return Repeat.LatticeCluster(clusterbasis, eltlattice, clusterelements)
end
export hex7RGB

function hex9RGB()
    clusterelements = SVector((-1, -1), (-1, 0), (-2, 1), (0, -1), (0, 0), (-1, 1), (1, -1), (1, 0), (0, 1))
    eltlattice = Repeat.HexBasis1()
    clusterbasis = Repeat.LatticeBasis((3, 0), (0, 3))
    return Repeat.LatticeCluster(clusterbasis, eltlattice, clusterelements)
end
export hex9RGB

function hex12RGB()
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
export hex12RGB

function hex19RGB()
    clusterelements = SVector(
        (0, 0), # ring 0
        (-1, 0),(0, -1),(1, -1),(1, 0),(0, 1),(-1, 1), # ring1
        (-2, 0),(-1, -1),(0, -2),(1, -2),(2, -2),(2, -1),(2, 0),(1, 1),(0, 2),(-1, 2),(-2, 2),(-2, 1) # ring2
        )
    eltlattice = Repeat.HexBasis3()
    clusterbasis = Repeat.LatticeBasis((5, 0), (-2, 4))
    return Repeat.LatticeCluster(clusterbasis, eltlattice, clusterelements)
end
export hex19RGB