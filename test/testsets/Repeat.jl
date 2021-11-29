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

    #spherepoint tests
    @test isapprox(Repeat.Lenslets.spherepoint(1,π/2,0.0), [0.0,1.0,0.0])
    @test isapprox(Repeat.Lenslets.spherepoint(1,0.0,π/2), [1.0,0.0,0.0])
    @test isapprox(Repeat.Lenslets.spherepoint(1,0,0.0), [0.0,0.0,1.0])
    @test isapprox(Repeat.Lenslets.spherepoint(1,0.0,π/4), [sqrt(2)/2,0.0,sqrt(2)/2])


    """ Create a LatticeCluser with three elements at (0,0),(-1,0),(-1,1) coordinates in the HexBasis1 lattice"""
    function hex3cluster()
        clusterelts = SVector((0,0),(-1,0),(-1,1))
        eltlattice = Repeat.HexBasis1()
        clusterbasis = Repeat.LatticeBasis(( -1,2),(2,-1))
        return Repeat.LatticeCluster(clusterbasis,eltlattice,clusterelts)
    end

    @test [-1 2;2 -1] == Repeat.basismatrix(Repeat.clusterbasis(hex3RGB()))

    function basistest(a::Repeat.AbstractLatticeCluster)
        return Repeat.clusterbasis(a)
    end

    @test basistest(hex3cluster()) == basistest(hex3RGB())  

    #LatticeCluster testset
    cluster = Repeat.Lenslets.hex9()

    for iter in 1:100
        (i,j) = rand.((1:1000,1:1000))
        coords,tileindex = Repeat.cluster_coordinates_from_tile_coordinates(cluster,i,j)
        reconstructed = Repeat.tilecoordinates(cluster,coords...,tileindex)
        @test all((i,j) .== reconstructed)
    end
 
    #verify that the 0,0 cluster is correct
    for (index,element) in pairs(Repeat.clusterelements(cluster))
        coords, tileindex = Repeat.cluster_coordinates_from_tile_coordinates(cluster, element...)
        @test all(coords .== 0)
        @test tileindex == index
    end
end