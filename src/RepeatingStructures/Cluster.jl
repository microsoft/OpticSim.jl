"""A Cluster is a repeating pattern of indices defined in a lattice called the element lattice. The cluster of elements has its own lattice, which by definition is different from the element lattice unless each cluster consists of a single element. For example this defines a Cluster of three hexagonal elements:

Cluster()"""
struct LatticeCluster{N1, N, B1<:Basis,B2<:Basis}
    clusterbasis::B1
    elementbasis::B2

    clusterelements::SVector{N1,NTuple{N,Int}} #vector containing the lattice indices of the elements in the cluster. Each column is one lattice coordinate. These indices are assumed to be offsets from the origin of the lattice.
end

"""returns the positions of every element in a cluster given the cluster indices"""
function Base.getindex(A::LatticeCluster{N1,N,B1,B2}, indices::Vararg{Int, N}) where{N1,N,T,B1<:Basis{N,T},B2<:Basis{N,T}} 
    clusteroffset = A.clusterbasis[indices]
    temp = MMatrix{N,N1,T}(undef)
    for i in 1:N1
        temp[:,i] = A.elementbasis[A.clusterelements[i]] + clusteroffset
    end
    return SMatrix{N,N1,T}(temp) #positions of the cluster elements, not the lattice coordinates.
end

"""Can access any cluster in a cluster lattice so the range of indices is unlimited"""
function Base.size(::LatticeCluster{N1,N}) where{N1,N}
    return ntuple((i)->Base.IsInfinite(),N)
end

Base.setindex!(A::Basis{N}, v, I::Vararg{Int, N}) where{N} = nothing #can't set lattice points. Might want to throw an exception instead.

"""
May want to have many properties associated with the elements in a cluster, which is why properties is represented as a DataFrame. The DataFrame in the properties field should have as many rows as there are elements in a cluster. At a minimum it must have a :Color, :Name, and :Lenslet column.

Example:
lenslets = DataFrame()
"""
struct LensletCluster
    lattice::LatticeCluster
    properties::DataFrame
end

function hex3RGB()
    clusterelements = SVector((0,0),(-1,0),(-1,-1))
    colors = [color("red"),color("green"),color("blue")]
    names = ["R","G","B"]
    eltlattice = HexBasis1()
    clusterbasis = LatticeBasis((0,2),(2,-1))
    lattice = LatticeCluster{3,2,HexBasis1,LatticeBasis}(clusterbasis,eltlattice,clusterelements)
    properties =  DataFrame(Colors = colors, Names = names)
    return LensletCluster(lattice,properties)
end
export hex3RGB