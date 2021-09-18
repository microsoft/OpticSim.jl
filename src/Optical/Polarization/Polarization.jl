module Polarization

using StaticArrays
using LinearAlgebra
using ..OpticSim
import Unitful
import Unitful.DefaultSymbols

abstract type AbstractPolarization{T<:Real} end
export AbstractPolarization

struct Chipman{T<:Real} <: AbstractPolarization{T}
    electricfieldvector::SVector{3,T}
    pmatrix::SMatrix{3,3,Complex{T},9} #3D polarization matrix
    
    Chipman{T}() where{T<:Real}= new{T}(SVector{3,T}(T(1),T(0),T(0)), SMatrix{3,3,Complex{T},9}(T(0) + T(0)im,T(0) + T(0)im,T(0) + T(0)im,T(0) + T(0)im,T(0) + T(0)im,T(0) + T(0)im,T(0) + T(0)im,T(0) + T(0)im,T(0) + T(0)im))
end
export Chipman

pmatrix(a::Chipman) = a.pmatrix
electricfieldvector(a::Chipman) = a.electricfieldvector

struct NoPolarization{T<:Real} <: AbstractPolarization{T}
end

"""Computes the O matrix (Polarized Light and Optical Systems by Chimpan et al, chapter 9) which represents the s-p coordinate frame for an interface. In OpticSim surface normals are defined positive pointing outward from a solid. In Chipman the surface normal is defined positive on the exiting side of an interface. Need to take this account when passing `surfacenormal` argument."""
function localtoworld(surfacenormal::AbstractVector{T},incidentvector::AbstractVector{T}) where{T<:Real}
    k̂ = normalize(incidentvector)
    ŝ = normalize(k̂ × surfacenormal)
    p̂ = incidentvector × surfacenormal
    return SMatrix{3,3,T,9}(ŝ[1],ŝ[2],ŝ[3],p̂[1],p̂[2],p̂[3],k̂[1],k̂[2],k̂[3])
end

worldtolocal(surfacenormal::AbstractVector{T},incidentvector::AbstractVector{T}) where{T<:Real} = localtoworld(surfacenormal,incidentvector)'

"""For Fresnel reflection need to compute reflected and/or refracted P matrix. This looks like:
see if normal is on the """

function composepmatrix(interface::OpticalInterface{T}, ray::OpticalRay{T,N,Polarization.Chipman{T}}) where{T<:Real,N}
    polarizationinput = polarization(ray)
    pinput = pmatrix(polarizationinput)
    evector = electricfieldvector(polarizationinput)

    #compute worldtolocal for incident ray
    plocal = pinputmatrix*worldtolocal(localtoworld(normal(interface),direction(ray)))
    #compute s,p components for reflected and transmitted values these are Jones matrix values, potentially complex.
    #compute reflected,refracted transformation

end



end #module
export Polarization