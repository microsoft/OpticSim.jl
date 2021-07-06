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

Returns the n neighbors surrounding centerpoint, excluding centerpoint
```
neighbors(::Type{B},centerpoint::Tuple{T,T},neighborhoodsize::Int) where{T<:Real,B<:Basis}
```
Returns the lattice basis vectors that define the lattice
```
basis(a::S) where{S<:Basis}
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




# """returns a Tuple = (2D coordinates of emitter in the lattice,emitter)"""
# emitter(a::L, i::Int, j::Int) where {L<:InfiniteLattice} = (a.origin + i * a.e1 + j * a.e2, a.emitter)

# emitterpitch_e₁(a::L) where {L<:Lattice} = norm(a.e₁)
# emitterpitch_e₂(a::L) where {L<:Lattice} = norm(a.e₂)

# struct SizedLattice{T<:Real,L<:InfiniteLattice{T}} <: FiniteLattice{T,L}
#     lattice::L
#     e₁::Int
#     e₂::Int

#     function SizedLattice(lattice::L, e₁::Int, e₂::Int) where {T<:Real,L<:InfiniteLattice{T}}
#         @assert e₁ >= T(0)
#         @assert e₂ >= T(0)
#         return new{T,L}(lattice, e₁, e₂)
#     end
# end

# emitterpitchx(a::SizedLattice) = emitterpitchx(norm(a.lattice.e₁))
# emitterpitchy(a::SizedLattice) = emitterpitchy(norm(a.lattice.e₂))

# function emitter(a::FiniteLattice{T,L}, i::Int, j::Int) where {T<:Real,L<:InfiniteLattice{T}}
#     @assert T(0) <= i <= a.e₁emitters
#     @assert T(0) <= j <= a.e₂emitters
#     return emitter(a.lattice, i, j) #because L is of type InfiteLattice and not FiniteLattice this call will not cause infinite recursion.
# end

# struct PrimitiveLattice{S<:AbstractSpectrum,P<:AbstractAngularPowerDistribution,T<:Real} <: InfiniteLattice{T}
#     origin::SVector{2,T}
#     e₁::SVector{2,T}
#     e2::SVector{2,T}
#     emitter::PlanarEmitter{S,P,T}

#     PrimitiveLattice(e₁::SVector{2,T}, e₂::SVector{2,T}, emitter::PlanarEmitter{S,P,T}; origin = SVector{2,T}(T(0), T(0))) where {S<:AbstractSpectrum,P<:AbstractAngularPowerDistribution,T<:Real} = new{S,P,T}(e₁, e₂, emitter, origin)
# end

# hexagonalemitterlattice(pitch::T, emitter::PlanarEmitter{S,P,T}; origin = SVector{2,T}(T(0), T(0))) where {S<:AbstractSpectrum,P<:AbstractAngularPowerDistribution,T<:Real} = PrimitiveLattice{S,P,T}(origin, SVector{2,T}(emitterpitch, T(0)), SVector{2,T}(emitterpitch / T(2), emitterpitch * sqrt(T(3)) / T(2)), emitter)

# squareemitterlattice(emitterpitch::T, emitter::PlanarEmitter{S,P,T}) where {S<:AbstractSpectrum,P<:AbstractAngularPowerDistribution,T<:Real} = PrimitiveLattice{S,P,T}(SVector{2,T}(emitterpitch, T(0)), SVector{2,T}(T(0), emitterpitch), emitter)

# #This needs work. We probably want some kind of hierarchical lattice with superemitter and subemitter functions.
# struct RGBPixelLattice{R<:AbstractSpectrum,G<:AbstractSpectrum,B<:AbstractSpectrum,P<:AbstractAngularPowerDistribution,T<:Real} <: InfiniteLattice{T}
#     red::PrimitiveLattice{R,P,T}
#     green::PrimitiveLattice{G,P,T}
#     blue::PrimitiveLattice{B,P,T}

#     function RGBPixelLattice(::Type{R}, ::Type{G}, ::Type{B}, ::Type{P}, subemitterwidth::T, rgbemitterpitch::T) where {R<:AbstractSpectrum,G<:AbstractSpectrum,B<:AbstractSpectrum,P<:AbstractAngularPowerDistribution,T<:Real}
#         @assert 3 * subemitterwidth <= rgbemitterpitch

#         subemitterpitch = rgbemitterpitch / T(3.0)
#         e₁ = SVector{2,T}(rgbemitterpitch, T(0))
#         e2 = SVector{2,T}(T(0), rgbemitterpitch)

#         return new{R,G,B,P,T}(PrimitiveLattice(e1, e2, PlanarEmitter{R,P,T}(subemitterwidth, rgbemitterpitch)), PrimitiveLattice(e1, e2, PlanarEmitter{G,P,T}(subemitterwidth, rgbemitterpitch), origin = SVector{2,T}(subemitterpitch, T(0.0))), PrimitiveLattice(e1, e2, PlanarEmitter{B,P,T}(subemitterwidth, rgbemitterpitch), origin = SVector{2,T}(2 * subemitterpitch, T(0.0))))
#     end
# end
