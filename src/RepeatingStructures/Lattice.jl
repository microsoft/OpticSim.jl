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

Base.getindex(A::Basis, indices::Vararg{Int, N}) where{N} = basis(A)*SVector{N,Int}(indices)

Base.setindex!(A::Basis{N,T}, v, I::Vararg{Int, N}) where{T,N} = nothing #can't set lattice points. Might want to throw an exception instead.

struct LatticeBasis{N,T<:Real} <: Basis{N,T}
    basisvectors::SMatrix{N,N,T}

    """ Convenience constructor that lets you use tuple arguments to describe the basis instead of SVector 
    
    Example:
    ```
    basis = LatticeBasis{2}((1.0,0.0),(0.0,1.0)) #basis for rectangular lattice
    ```
    """
    function LatticeBasis(vectors::Vararg{SVector{N,T}}) where{T,N}  
        temp = MMatrix{N,N,T}(undef)
 
        for (j,val) in pairs(vectors)
            for i in 1:N
                temp[i,j] = val[i]
            end
        end
        
        return new{N,T}(SMatrix{N,N,T}(temp))
    end

     LatticeBasis(vectors::SMatrix{N,N,T}) where{N,T<:Real} = new{N,T}(vectors)
end
export LatticeBasis

basis(a::LatticeBasis) = a.basisvectors

"""Can access any point in a lattice so the range of indices is unlimited"""
function Base.size(a::LatticeBasis{N,T}) where{N,T}
    return ntuple((i)->Base.IsInfinite(),N)
end