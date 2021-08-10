"""A Cluster is a repeating pattern of indices defined in a lattice called the element lattice. The cluster of elements has its own lattice, which by definition is different from the element lattice unless each cluster consists of a single element. For example this defines a Cluster of three hexagonal elements:

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
struct LatticeCluster{N1, N, T<:Real, B1<:Basis{N,Int},B2<:Basis{N,T}}
    clusterbasis::B1 #this basis defines the offsets of the clusters
    elementbasis::B2 #this basis defines the underlying lattice

    clusterelements::SVector{N1,NTuple{N,Int}} #vector containing the lattice indices of the elements in the cluster. Each column is one lattice coordinate. These indices are assumed to be offsets from the origin of the lattice.

    LatticeCluster(clusterbasis::B1,eltlattice::B2,clusterelements::SVector{N1,NTuple{N,Int}}) where{N1,N,T<:Real,B1<:Basis{N,Int},B2<:Basis{N,T}} = new{N1,N,T,B1,B2}(clusterbasis,eltlattice,clusterelements)
end
export LatticeCluster

clusterelements(a::LatticeCluster) = a.clusterelements
export clusterelements

"""returns the positions of every element in a cluster given the cluster indices"""
function Base.getindex(A::LatticeCluster{N1,N,T,B1,B2}, indices::Vararg{Int, N}) where{N1,N,T,B1<:Basis{N,Int},B2<:Basis{N,T}} 
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

Base.setindex!(A::Basis{N}, v, I::Vararg{Int, N}) where{N} = nothing #can't set lattice points. Might want to throw an exception instead.

""" returns the lattice indices of the elements in the cluster. These are generally not the positions of the elements"""
function clustercoordinates(a::LatticeCluster{N1,N},indices::Vararg{Int,N}) where{N1,N}
    clusteroffset = a.clusterbasis[indices...]
    temp = MMatrix{N,N1,Int}(undef)
    for i in 1:N1
        temp[:,i] = SVector{N,Int}(a.clusterelements[i]) + clusteroffset
    end
    return SMatrix{N,N1,Int}(temp) #positions of the cluster elements, not the lattice coordinates.
end
export clustercoordinates


"""
May want to have many properties associated with the elements in a cluster, which is why properties is represented as a DataFrame. The DataFrame in the properties field should have as many rows as there are elements in a cluster. At a minimum it must have a :Color and a :Name column.

Example:
lenslets = DataFrame()
"""
struct ClusterWithProperties{N1,N,T}
    cluster::LatticeCluster{N1,N,T}
    clusterindices::SVector{NTuple{N,T}}
    properties::DataFrame
end

Base.getindex(A::ClusterWithProperties{N1,N,T}, indices::Vararg{Int,N}) where{N1,N,T} = cluster(A)[indices...]
Base.setindex!(A::ClusterWithProperties, v, I::Vararg{Int, N}) where{N} = nothing #can't set lattice points. Might want to throw an exception instead.

properties(a::ClusterWithProperties) = a.properties
export properties
cluster(a::ClusterWithProperties) = a.cluster
export cluster
clustercoordinates(a::ClusterWithProperties,indices...) = clustercoordinates(lattice(a),indices...)

function hex3cluster()
    clusterelts = SVector((0,0),(-1,0),(-1,1))
    eltlattice = HexBasis1()
    clusterbasis = LatticeBasis(( -1,2),(2,-1))
    return LatticeCluster(clusterbasis,eltlattice,clusterelts)
end
export hex3cluster

""" Create a cluster with properties with three types of elements, R,G,B """
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

""" Create a cluster with properties with four types of elements, R,G,B,W """
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


