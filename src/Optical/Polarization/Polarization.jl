module Polarization

using StaticArrays
using LinearAlgebra
using ..OpticSim

abstract type AbstractPolarization{T<:Real} end

struct Chipman{T<:Real} <: AbstractPolarization{T}
    P::SMatrix{3,3,Complex{T},9} #3D polarization matrix
    Chipman{T}() = new{T}(SMatrix{3,3,Complex{T},9}(0.0 + 0.0im,0.0 + 0.0im,0.0 + 0.0im,0.0 + 0.0im,0.0 + 0.0im,0.0 + 0.0im,0.0 + 0.0im,0.0 + 0.0im,0.0 + 0.0im)) where{T<:Real}
end

"""Computes the P matrix (Polarized Light and Optical Systems by Chimpan et al, chapter 9) which represents the s-p coordinate frame for an interface. In OpticSim surface normals are defined positive pointing outward from a solid. In Chipman the surface normal is defined positive on the exiting side of an interface. Need to take this account when passing `surfacenormal` argument."""
function localtoworld(surfacenormal::AbstractVector{T},incidentvector::AbstractVector{T}) where{T<:Real}
    k̂ = normalize(incidentvector)
    ŝ = normalize(k̂ × surfacenormal)
    p̂ = incidentvector × surfacenormal
    return SMatrix{3,3,T,9}(ŝ[1],ŝ[2],ŝ[3],p̂[1],p̂[2],p̂[3],k̂[1],k̂[2],k̂[3])
end

worldtolocal(loctoworld::SMatrix{3,3,T,9}) = loctoworld'

"""For Fresnel reflection need to compute reflected and/or refracted P matrix. This looks like:
see if normal is on the """

function composepmatrix(pinputmatrix::Chipman{T},interface::OpticalInterface{T}, ray::OpticalRay{T}) where{T<:Real}
    #compute worldtolocal for incident ray
    plocal = pinputmatrix*worldtolocal(localtoworld(normal(interface),direction(ray)))
    #compute s,p components for reflected and transmitted values these are Jones matrix values, potentially complex.
    #compute reflected,refracted transformation

end

end #module