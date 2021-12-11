# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

#Various sizes of hexagonal clusters that are useful in lenslet design. To see what they look like use the functions in HexTilings to make lattices with color properties and then draw them with Vis.draw
#example: Vis.draw(hex3RGB())

function hex3(scale::T =  1.0) where{T<:Real}
    clusterelements = SVector((0, 0), (-1, 0), (-1, 1))
    eltlattice = HexBasis1(scale)
    clusterbasis = LatticeBasis((-1, 2), (2, -1))
    return LatticeCluster(clusterbasis, eltlattice, clusterelements)

end
export hex3

function hex4(scale::T =  1.0) where{T<:Real}
    clusterelements = SVector((-1, 0), (-1, 1), (0, 0), (0, 1))
    eltlattice = HexBasis1(scale)
    clusterbasis = LatticeBasis((2, 0), (0, 2))
    return LatticeCluster(clusterbasis, eltlattice, clusterelements)
end
export hex4

hex7(scale::T) where{T<:Real} = hexn(1,scale)
export hex7


function hex9(scale::T =  1.0) where{T<:Real}
    clusterelements = SVector((-1, -1), (-1, 0), (-2, 1), (0, -1), (0, 0), (-1, 1), (1, -1), (1, 0), (0, 1))
    eltlattice = HexBasis1(scale)
    clusterbasis = LatticeBasis((3, 0), (0, 3))
    return LatticeCluster(clusterbasis, eltlattice, clusterelements)
end
export hex9

function hex12(scale::T =  1.0) where{T<:Real}
    clusterelements = SVector(
        (-1, 1),(0, 1),
        (-1, 0),(0, 0),(1, 0),
        (-1, -1),(0, -1),(1, -1),(2, -1),
        (0, -2),(1, -2),(2, -2)
    )

    eltlattice = HexBasis3(scale)
    clusterbasis = LatticeBasis((2, 2), (-3, 2))
    return LatticeCluster(clusterbasis, eltlattice, clusterelements)
end
export hex12

#TODO fix cluster basis vectors.
"""WARNING: the cluster basis vectors are incorrect. Still need to figure out what these are. This will not tile properly"""
function hex18(scale::T =  1.0) where{T<:Real}
    throw(ErrorException("not yet working"))
    clusterelements = SVector(
        (-2,1),(-1, 1),(0, 1),(1,1),
        (-2,0),(-1, 0),(0, 0),(1, 0),(2,0),
        (-1, -1),(0, -1),(1, -1),(2, -1),
        (0, -2),(1, -2),(2, -2),
        (1,-3),(2,-3)
    )

    eltlattice = HexBasis3(scale)
    clusterbasis = LatticeBasis((2, 2), (-3, 2))
    return LatticeCluster(clusterbasis, eltlattice, clusterelements)
end
export hex18

hex19(scale::T = 1.0) where{T<:Real} = hexn(2,scale)
export hex19

hex37(scale::T = 1.0) where{T<:Real} = hexn(3,scale)
export hex37

""" function for creating symmetric clusters consisting of all lattice cells within regionsize of the origin. This includes hex1 (use regionsize = 0),hex7 (regionsize = 1),hex19 (region size = 2),hex37 (regionsize = 3), etc."""
function hexn(regionsize::Int64,scale::T = 1.0) where{T<:Real}
    eltlattice = HexBasis1(scale)
    clusterbasis = LatticeBasis((2*regionsize+1, -regionsize), (regionsize, regionsize+1))
    return LatticeCluster(clusterbasis, eltlattice, region(HexBasis1,(0,0),regionsize))
end
export hexn

