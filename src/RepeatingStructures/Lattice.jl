# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

""" Base class for defining points on a lattice. Lattice points are defined by basis vectors eᵢ and lattice indices indexᵢ:

point = Σᵢ indexᵢ*eᵢ

Example: define a 2D basis. 

```julia
julia> a = LatticeBasis([1,2],[3,4])
LatticeBasis{2, Int64}(SVector{2, Int64}[[1, 2], [3, 4]])
```

Array indexing is used to generate a lattice point: a[1,2] = 1*[1,2] + 2*[3,4] 

```julia
julia> a[1,2]
2-element SVector{2, Int64} with indices SOneTo(2):
  7
 10
```

Subtypes supporting the Basis interface should implement these functions:

Returns the neighbors in ring n surrounding centerpoint, excluding centerpoint
```
neighbors(::Type{B},centerpoint::Tuple{T,T},neighborhoodsize::Int) where{T<:Real,B<:Basis}
```
Returns the lattice basis vectors that define the lattice
```
basis(a::S) where{S<:Basis}
```
Returns the vertices of the unit polygon that tiles the plane for the basis
```
tilevertices(a::S) where{S<:Basis}
```
"""
abstract type Basis{N,T<:Real} end
export Basis


"""wrote this function because type inference couldn't handle argument types of Basis{N,T}. Not sure why"""
numbertype(a::Basis{N,T}) where{N,T} = T

function Base.getindex(A::Basis, indices::Vararg{Int, N}) where{N}
    T = numbertype(A)
    sum = MVector{N,T}(zeros(T,N))
    basisvecs = basis(A)
    for (i,index) in pairs(indices)
        sum += basisvecs[i]*index
    end
    return SVector{N,T}(sum)
end

Base.setindex!(A::Basis{N,T}, v, I::Vararg{Int, N}) where{T,N} = nothing #can't set lattice points. Might want to throw an exception instead.

struct LatticeBasis{N,T<:Real} <: Basis{N,T}
    basisvectors::SVector{N,SVector{N,T}}

    """ Convenience constructor that lets you use Vector arguments to describe the basis instead of SVector 
    
    Example:
    ```
    basis = LatticeBasis{2}([1.0,0.0],[0.0,1.0]) #basis for rectangular lattice
    ```
    """
    function LatticeBasis(vectors::Vararg{Vector{T},N}) where{T,N}
        dim = length(vectors[1])
        @assert N == dim
        for i in 2:length(vectors)
            @assert length(vectors[i])==dim "Vectors for the lattice basis were not all the same dimension"
        end
        
        temp = MVector{N,SVector{N,T}}(undef) #MVector is slightly faster and has slightly fewer allocations than Vector
 
        for (i,val) in pairs(vectors)
            temp[i] = SVector{N,T}(ntuple((j)->val[j],N)) #using ntuple is significantly faster than val... No idea why this should be true.
        end
        
        # return new{N,T}(SVector{N,SVector{N,T}}(temp...))
        return new{N,T}(SVector{N,SVector{N,T}}(ntuple((j)->temp[j],N)))
    end

    LatticeBasis(vectors::Vararg{SVector{N,T},N}) where{N,T} = new{N,T}(SVector{N,SVector{N,T}}(vectors)) #this function is considerably faster than the other constructor
end
export LatticeBasis

basis(a::LatticeBasis) = a.basisvectors

function Base.size(a::LatticeBasis{N,T}) where{N,T}
    return ntuple((i)->Base.IsInfinite(),N)
end