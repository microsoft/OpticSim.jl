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
    @test isapprox(Repeat.Multilens.spherepoint(1,π/2,0.0), [0.0,1.0,0.0])
    @test isapprox(Repeat.Multilens.spherepoint(1,0.0,π/2), [1.0,0.0,0.0])
    @test isapprox(Repeat.Multilens.spherepoint(1,0,0.0), [0.0,0.0,1.0])
    @test isapprox(Repeat.Multilens.spherepoint(1,0.0,π/4), [sqrt(2)/2,0.0,sqrt(2)/2])


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
    cluster = Repeat.Multilens.hex9()

    #generate many tile coordinates. Compute the cluster index and tile index in that cluster for each tile coordinate.
    for iter in 1:100
        (i,j) = rand.((1:1000,1:1000))
        coords,tileindex = Repeat.cluster_coordinates_from_tile_coordinates(cluster,i,j)
        reconstructed = Repeat.tilecoordinates(cluster,coords...,tileindex)
        @test all((i,j) .== reconstructed)
    end

    function testassignment()
        #test assignment of eyebox numbers to RGB clusters
        rgb_cluster = Repeat.Multilens.hex12RGB()
        cluster_coords = map(x->Tuple(x),eachcol(Repeat.clustercoordinates(rgb_cluster,0,0))) #create the cluster coordinates corresponding to each of the tiles in the cluster
        eyeboxnumbers = (1,1,2,1,2,2,3,3,3,4,4,4) #correct eyebox number assignment for the tiles in the cluster
        for (index,coord) in enumerate(cluster_coords)
            boxnum = eyeboxnumbers[index]
            num = Repeat.Multilens.eyebox_number(coord,rgb_cluster,4)
            if num != boxnum
                return false
            end
        end
        return true
    end

    @test testassignment()

    """This verifies that the eyebox number assignment is invariant to the offset of the cluster"""
    function test_eyebox_assignment()
        subdivisions = (3,1)
        cluster = Repeat.Multilens.hex9RGB()
        eye_box_frame = Geometry.identitytransform()
        eye_box = (10*Unitful.DefaultSymbols.mm,9*Unitful.DefaultSymbols.mm)
        lattice_coordinates = Repeat.clusterelements(cluster)
    
       #compute subdivided eyebox polygons and assign to appropriate lenslets
       eyeboxpoly::SMatrix{3,4}  =  Unitful.DefaultSymbols.mm * (eye_box_frame * Repeat.Multilens.eyeboxpolygon(ustrip.(Unitful.DefaultSymbols.mm,eye_box)...)) #four corners of the eyebox frame which is assumed centered around the positive Z axis. Transformed to the eyeballframe. Have to switch back and forth between Unitful and unitless quantities because Transform doesn't work with Unitful values.
        subdivided_eyeboxpolys::Vector{SMatrix{3,4}} = Repeat.Multilens.compute_lenslet_eyebox_data(eye_box_frame,eyeboxpoly,subdivisions)
    
         
        # println("lenslet eyeboxes")
        # println(lenslet_eye_boxes)
    
        #assign each lenslet the index into subdivided_eyeboxpolys that corresponds to the subdivided eyebox polygon the lenslet is supposed to cover.
        lenslet_eyebox_numbers = Repeat.Multilens.eyebox_number.(lattice_coordinates,Ref(cluster),Ref(length(subdivided_eyeboxpolys)))
        # println(lattice_coordinates)
        # println("lenslet eybox nums $lenslet_eyebox_numbers")
        offset_latticecoords = [tuple(x...) for x in eachcol(Repeat.clustercoordinates(cluster,1,1))]
        # println("offset coords $offset_latticecoords")
        offset_eyebox_numbers = Repeat.Multilens.eyebox_number.(offset_latticecoords,Ref(cluster),Ref(length(subdivided_eyeboxpolys)))
        # println("eyebox numbers of offset cluster $offset_eyebox_numbers")
        
        lenslet_eyebox_numbers == offset_eyebox_numbers
    end

    @test test_eyebox_assignment()

  

    #verify that the 0,0 cluster is correct
    for (index,element) in pairs(Repeat.clusterelements(cluster))
        coords, tileindex = Repeat.cluster_coordinates_from_tile_coordinates(cluster, element...)
        @test all(coords .== 0)
        @test tileindex == index
    end

    #verify that choosecluster assertion doesn't fire incorrectly
    using Unitful,Unitful.DefaultSymbols
    function generate_clusters()
        freq = Vector{Int64}(undef,0)
        subdivs = Vector{Tuple{Int64,Int64}}(undef,0)
        areas = Vector(undef,0)

        try
            for cycles in 15:30
                OpticSim.Repeat.Multilens.system_properties(15mm,(10mm,9mm),(100°,70°),3.5mm,.1,cycles)
            end
        catch err #if any errors then failure
            return false
        end
        return true
    end

    @test generate_clusters() #shouldn't get assertion failures for any of the frequencies between 15:30
end