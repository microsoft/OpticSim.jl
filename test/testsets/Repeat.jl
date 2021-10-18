# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

@testset "Repeat" begin
    function hex3RGB()
        clusterelements = SVector((0,0),(-1,0),(-1,1))
        colors = [colorant"red",colorant"green",colorant"blue"]
        names = ["R","G","B"]
        eltlattice = Repeat.HexBasis1()
        clusterbasis = Repeat.LatticeBasis(( -1,2),(2,-1))
        lattice = Repeat.LatticeCluster(clusterbasis,eltlattice,clusterelements)
        properties =  DataFrame(Color = colors, Name = names)
        return Repeat.ClusterWithProperties(lattice,properties)
    end

    @test [-1 2;2 -1] == Repeat.basismatrix(Repeat.clusterbasis(hex3RGB()))
end