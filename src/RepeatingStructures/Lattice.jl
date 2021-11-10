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

Subtypes supporting the AbstractBasis interface should implement these functions:

Returns the neighbors in ring n surrounding centerpoint, excluding centerpoint
```
neighbors(::Type{B},centerpoint::Tuple{T,T},neighborhoodsize::Int) where{T<:Real,B<:AbstractBasis}
```
Returns the lattice basis vectors that define the lattice
```
basis(a::S) where{S<:AbstractBasis}
```
Returns the vertices of the unit polygon that tiles the plane for the basis
```
tilevertices(a::S) where{S<:AbstractBasis}
```
"""
abstract type AbstractBasis{N,T <: Real} end
export AbstractBasis

function Base.getindex(A::B1, indices::Vararg{Int,N}) where {N,T,B1 <: AbstractBasis{N,T}}
    return basismatrix(A) * SVector{N,Int}(indices)
end

#     temp::SVector{N,T} = (basis(A)*SVector{N,Int}(indices))::SVector{N,T}
#     return temp
# end

Base.setindex!(A::AbstractBasis{N,T}, v, I::Vararg{Int,N}) where {T,N} = nothing # can't set lattice points. Might want to throw an exception instead.

"""computes the tile vertices at the offset lattice coordinate"""
tilevertices(lattice::AbstractBasis, latticecoords::Union{AbstractVector,NTuple}) = lattice[latticecoords...] .+ tilevertices(lattice)

struct LatticeBasis{N,T <: Real} <: AbstractBasis{N,T}
    basisvectors::SMatrix{N,N,T}

    """ Convenience constructor that lets you use tuple arguments to describe the basis instead of SVector 
    
    Example:
    ```
    basis = LatticeBasis{2}((1.0,0.0),(0.0,1.0)) #basis for rectangular lattice
    ```
    This allocates 48 bytes for a 2D basis in Julia 1.6.2. It  shouldn't allocate anything but not performance critical.
    """
    function LatticeBasis(vectors::Vararg{NTuple{N,T},N}) where {T,N}  
        temp = MMatrix{N,N,T}(undef)
        
        for (j, val) in pairs(vectors)
            for i in 1:N
                temp[i,j] = val[i]
            end
        end
        
        return new{N,T}(SMatrix{N,N,T}(temp))
    end

     LatticeBasis(vectors::SMatrix{N,N,T}) where {N,T <: Real} = new{N,T}(vectors)
end
export LatticeBasis

function basismatrix(a::LatticeBasis{N,T})::SMatrix{N,N,T,N * N} where {N,T} 
    return a.basisvectors
end

"""Can access any point in a lattice so the range of indices is unlimited"""
function Base.size(a::LatticeBasis{N,T}) where {N,T}
    return ntuple((i) -> Base.IsInfinite(), N)
end

# These functions are used to find the lattice tiles that lie inside a polygonal shape

"""computes the matrix that transforms the basis set into the canonical basic vectors, eᵢ"""
basisnormalization(a::Repeat.AbstractBasis) = Matrix(inv(Repeat.basismatrix(a))) # unfortunately LazySets doesn't work with StaticArrays. If you return a static Array is messes with their functions.

"""warp the matrix into a coordinate frame where the basis vectors are canonical, eᵢ"""
function normalizedshape(a::Repeat.AbstractBasis, shape::LazySets.VPolygon) 
    normalizer = basisnormalization(a)
    return normalizer * shape
end
export normalizedshape

""" compute a conservative range of lattice indices that might intersect containingshape"""
function latticebox(containingshape::LazySets.VPolygon, lattice::Repeat.AbstractBasis)
    warpedshape = normalizedshape(lattice, containingshape)
    box = LazySets.interval_hull(warpedshape)
    return LazySets.Hyperrectangle(round.(box.center), (1, 1) .+ ceil.(box.radius)) # center the hyperrectangle on an integer point and enlarge the radius to take this into account
end
export latticebox


"""Returns the lattice coordinates of all tiles that lie at least partly inside containingshape. 

Compute the lattice tiles with non-zero intersection with containingshape. First compute a transformation that maps the lattice basis vectors to canonical unit basis vectors eᵢ (this is the inverse of the lattice basis matrix). Then transform containingshape into this coordinate frame and compute a bounding box. Unit steps along the coordinate axes in this space represent unit *lattice* steps in the original space. This makes it simple to determine coordinate bounds in the original space. Then test for intersection in the original space."""
function tilesinside(containingshape::LazySets.VPolygon, lattice::Repeat.AbstractBasis{2,T}) where {T}
    box = latticebox(containingshape, lattice)
    
    coords = Int64.(box.radius)
    tilevertices = LazySets.VPolygon(Matrix(Repeat.tilevertices(lattice))) # VPolygon will accept StaticArrays but other LazySets function will barf.
    result = Vector{NTuple{2,Int64}}(undef, 0)
    
    for i in -coords[1]:coords[1]
        for j in -coords[2]:coords[2]       
            center = Vector(lattice[i,j])
            offsethex = LazySets.translate(tilevertices, center)

            if !isempty(offsethex ∩ containingshape)
                push!(result, (i, j))
            end
        end
    end
    return result
end
export tilesinside

"""The vertices of the containing shape are the columns of the matrix containingshape"""
tilesinside(containingshape::AbstractMatrix,lattice::Repeat.AbstractBasis) = tilesinside(LazySets.VPolygon(containingshape), lattice)

tilesinside(xmin::T,ymin::T,xmax::T,ymax::T,lattice::Repeat.AbstractBasis) where {T <: Real} = tilesinside(SMatrix{2,4,T}(xmin, ymin, xmin, ymax, xmax, ymax, xmax, ymin), lattice)

using Plots

""" to see what the objects look like in the warped coordinate frame use inv(basismatrix(lattice)) as the transform"""
function plotall(containingshape, lattice, transform=[1.0 0.0;0.0 1.0])

     garb = tilevertices.(lattice, tilesinside(containingshape, lattice))
    #  garb = tilesinside(containingshape, lattice)
    for g in garb println(g)
    end

    tiles = LazySets.VPolygon.([x -> tilevertices(lattice, x) for x in tilesinside(containingshape, lattice)])
    for tile in tiles
        plot!(transform * tile, aspectratio=1)
    end
    plot!(containingshape, aspectratio=1)
end
export plotall


function testtilesinside()
    # poly = LazySets.VPolygon([-3.0 3.0 3.0; -3.0 3.0 -3.0])
    hex = HexBasis1()
    poly = LazySets.VPolygon(3 * tilevertices(hex))
    tilesinside(poly, hex)
    plotall(poly, hex)
end
export testtilesinside
