"""A Cluster is a repeating pattern of indices defined in a lattice called the element lattice. The cluster of elements has its own lattice, which by definition is different from the element lattice unless each cluster consists of a single element. For example this defines a Cluster of three hexagonal elements:

Cluster()"""
struct LatticeCluster{N1, B1<:Basis{N,T},B2<:Basis{N,T}}
    clusterbasis::B1
    elementbasis::B2

    clusterelements::SVector{N1,NTuple{N,Int}} #vector containing the lattice indices of the elements in the cluster. Each column is one lattice coordinate. These indices are assumed to be offsets from the origin of the lattice.
end

"""returns the positions of every element in a cluster given the cluster indices"""
function Base.getindex(A::LatticeCluster{N1,B1{N,T},B2{N,T}}, indices::Vararg{Int, N}) where{N} 
    clusteroffset = A.clusterbasis[indices]
    temp = MMatrix{N,N1,T}(undef)
    for i in 1:N1
        temp[:,i] = A.elementbasis[A.clusterelements[i]] + clusteroffset
    end
    return SMatrix{N,N1,T}(temp) #positions of the cluster elements, not the lattice coordinates.
end

"""Can access any cluster in a cluster lattice so the range of indices is unlimited"""
function Base.size(::LatticeCluster{N,B1{N,T},B2{N,T}}) where{N,T}
    return ntuple((i)->Base.IsInfinite(),N)
end

Base.setindex!(A::Basis{N,T}, v, I::Vararg{Int, N}) where{T,N} = nothing #can't set lattice points. Might want to throw an exception instead.

"""
Example:
lenslets = DataFrame()
"""
struct LensletCluster
    lattice::LatticeCluster
    properties::DataFrame
end