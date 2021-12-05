# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

"""A Cluster is a repeating pattern of indices defined in a lattice called the element lattice. The cluster of elements has its own lattice, which by definition is different from the element lattice unless each cluster consists of a single element. A LatticeCluster subtype must implement these methods:

clusterbasis(a::AbstractLatticeCluster) returns the AbstractBasis object for the cluster. This is not generally the same as the basis for the underlying lattice.

elementbasis(a::AbstractLatticeCluster) returns the AbstractBasis object the underlying lattice of the cluster.

clusterelements(a::AbstractLatticeCluster) returns the lattice indices, represented in the underlying lattice basis, for each of the elements in a unit cluster cell.

clustersize(a::AbstractLatticeCluster) returns the number of elements in a unit cluster cell.

For example this defines a Cluster of three hexagonal elements:

```
julia> function hex3cluster()
    clusterelts = SVector((0,0),(-1,0),(-1,1))
    eltlattice = HexBasis1()
    clusterbasis = LatticeBasis(( -1,2),(2,-1))
    return LatticeCluster(clusterbasis,eltlattice,clusterelts)
end

julia> lattice = hex3cluster()
LatticeCluster{3, 2, Float64, LatticeBasis{2, Int64}, HexBasis1{2, Float64}}(LatticeBasis{2, Int64}([-1 2; 2 -1]), HexBasis1{2, Float64}(), [(0, 0), (-1, 0), (-1, 1)])

julia> lattice[1,1]  #returns the locations of the 3 elements in the cluster at the cluster offset of 1,1
2×3 SMatrix{2, 3, Float64, 6} with indices SOneTo(2)×SOneTo(3):
 1.0  -0.5        1.0
 1.0   0.133975  -0.732051

```
"""
abstract type AbstractLatticeCluster end



"""Basic lattic cluster type"""
struct LatticeCluster{N1, N, T<:Real, B1<:AbstractBasis{N,Int},B2<:AbstractBasis{N,T}} <: AbstractLatticeCluster
    clusterbasis::B1 #this basis defines the offsets of the clusters
    elementbasis::B2 #this basis defines the underlying lattice

    clusterelements::SVector{N1,NTuple{N,Int64}} #vector containing the lattice indices of the elements in the cluster. Each column is one lattice coordinate. These indices are assumed to be offsets from the origin of the lattice.

    LatticeCluster(clusterbasis::B1,eltlattice::B2,clusterelements::SVector{N1,NTuple{N,Int64}}) where{N1,N,T<:Real,B1<:AbstractBasis{N,Int64},B2<:AbstractBasis{N,T}} = new{N1,N,T,B1,B2}(clusterbasis,eltlattice,clusterelements)
end
export LatticeCluster

LatticeCluster(clusterbasis::B1,eltlattice::B2,clusterelements::SMatrix{N,N1,Int64}) where{N1,N,T<:Real,B1<:AbstractBasis{N,Int64},B2<:AbstractBasis{N,T}} = LatticeCluster(clusterbasis,eltlattice,SVector{N1}(reinterpret(reshape,NTuple{N,Int64},clusterelements)))


clusterelements(a::LatticeCluster) = a.clusterelements
export clusterelements
elementbasis(a::LatticeCluster) = a.elementbasis
clustersize(a::LatticeCluster) = length(a.clusterelements)
clusterbasis(a::LatticeCluster) = a.clusterbasis
export clusterbasis

""" Lattice clusters have integer basis vectors. These represent the lattice in unit steps of the underlying element basis, not the Euclidean distance between lattice cluster centers. For example, the lattice cluster basis vectors for a hex3 lattice are (-1,2),(2,-1), where the units of the basis vectors represent steps in the underlying element basis. To get a Euclidean distance measure between cluster centers you need to multiply by the size of the element basis diameter. In the case of lenslet layout you need to multiply the cluster diameter by the lenslet diameter."""
latticediameter(a::Repeat.AbstractLatticeCluster) =   latticediameter(Repeat.basismatrix(Repeat.clusterbasis(a)))

"""returns the positions of every element in a cluster given the cluster indices"""
function Base.getindex(A::LatticeCluster{N1,N,T,B1,B2}, indices::Vararg{Int, N}) where{N1,N,T,B1<:AbstractBasis{N,Int},B2<:AbstractBasis{N,T}} 
    clusteroffset = A.clusterbasis[indices...]
    temp = MMatrix{N,N1,T}(undef)
    for i in 1:N1
        temp[:,i] = A.elementbasis[A.clusterelements[i]...] + clusteroffset
    end
    return SMatrix{N,N1,T}(temp) #positions of the cluster elements, not the lattice coordinates.
    
end

"""Can access any cluster in a cluster lattice so the range of indices is unlimited"""
function Base.size(::LatticeCluster{N1,N}) where{N1,N}
    return ntuple((i)->Base.IsInfinite(),N)
end

"""returns the integer coordinates of the tile at tileindex in the cluster at for a cluster with location cᵢ,cⱼ"""
function tilecoordinates(cluster::LatticeCluster{N1,N},cᵢ,cⱼ,tileindex) where{N1,N}
    @assert tileindex <= N1
    clustercenter = basismatrix(clusterbasis(cluster))*SVector(cᵢ,cⱼ)
    tileoffset = cluster.clusterelements[tileindex]
    return clustercenter .+ tileoffset
end
export tilecoordinates

Base.setindex!(A::AbstractBasis{N}, v, I::Vararg{Int, N}) where{N} = nothing #can't set lattice points. Might want to throw an exception instead.

""" returns the lattice indices of the elements in the cluster. These are generally not the positions of the elements"""
function clustercoordinates(a::LatticeCluster{N1,N},indices::Vararg{Int,N}) where{N1,N}
    clusteroffset = a.clusterbasis[indices...]
    temp = MMatrix{N,N1,Int}(undef)
    for i in 1:N1
        temp[:,i] = SVector{N,Int}(a.clusterelements[i]) + clusteroffset
    end
    return SMatrix{N,N1,Int}(temp) #the lattice coordinates in the underlying element basis of the cluster elements
end
export clustercoordinates

"""Given the (i,j) coordinates of a tile defined in the the underlying lattice basis of elements of the cluster compute the coordinates (cᵢ,cⱼ) of the cluster containing the tile, and the tile number of the tile in that cluster"""
function cluster_coordinates_from_tile_coordinates(cluster::LatticeCluster{N1,N},i::Int,j::Int) where{N1,N}
    found = false
    bmatrix = Rational.(basismatrix(clusterbasis(cluster)))

    local clustercoords::SVector{}
    tileindex = 0 #initialize to illegal value
    for coords in eachcol(clustercoordinates(cluster, 0,0)) #don't know which tile in the cluster the i,j coords come from so try all of them until find one that yields an integer value for the cluster coordinate.
        tileindex += 1
        possiblecenter = SVector{N,Rational}((i,j)) .- Rational.(coords)
        clustercoords = bmatrix \  possiblecenter
        if all(1 .== denominator.(clustercoords)) #If coordinates are integer this is the correct cluster value. If coordinates are not integer keep trying.
            found = true    #every tile position i,j corresponds to some cluster (cᵢ,cⱼ) so will always find an answer
            break
        end
    end
    @assert found == true #should always find a cluster corresponding to any i,j. If not then something is seriously wrong.

    return Int64.(clustercoords),tileindex  #return coordinates of the cluster and the index of the tile within that cluster
end
export cluster_coordinates_from_tile_coordinates



"""
May want to have many properties associated with the elements in a cluster, which is why properties is represented as a DataFrame. The DataFrame in the properties field should have as many rows as there are elements in a cluster. At a minimum it must have a :Color and a :Name column.

Example: a lattice cluster of three hexagonal elements, labeled R,G,B with colors red,green,blue.

```
function hex3RGB()
    clusterelements = SVector((0,0),(-1,0),(-1,1)) #lattice indices of cluster elements represent in the element lattice basis
    colors = [colorant"red",colorant"green",colorant"blue"]
    names = ["R","G","B"]
    eltlattice = HexBasis1()
    clusterbasis = LatticeBasis(( -1,2),(2,-1)) #lattice indices, in the element lattice basis, representing the repeating pattern of the cluster
    lattice = LatticeCluster(clusterbasis,eltlattice,clusterelements)
    properties =  DataFrame(Color = colors, Name = names)
    return ClusterWithProperties(lattice,properties)
end
```
"""
struct ClusterWithProperties{N1,N,T} <: AbstractLatticeCluster
    cluster::LatticeCluster{N1,N,T}
    properties::DataFrame
end
export ClusterWithProperties

Base.getindex(A::ClusterWithProperties{N1,N,T}, indices::Vararg{Int,N}) where{N1,N,T} = cluster(A)[indices...]
Base.setindex!(A::ClusterWithProperties, v, I::Vararg{Int, N}) where{N} = nothing #can't set lattice points. Might want to throw an exception instead.

properties(a::ClusterWithProperties) = a.properties
export properties
cluster(a::ClusterWithProperties) = a.cluster
export cluster
clustercoordinates(a::ClusterWithProperties,indices...) = clustercoordinates(cluster(a),indices...)
elementbasis(a::ClusterWithProperties) = elementbasis(cluster(a))
export elementbasis
clustersize(a::ClusterWithProperties) = clustersize(a.cluster)
export clustersize
clusterbasis(a::ClusterWithProperties) = clusterbasis(a.cluster)
export clusterbasis
