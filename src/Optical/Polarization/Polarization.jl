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
       
    Chipman{T}() where{T<:Real}= new{T}(SVector{3,Complex{T}}(T(1),T(0),T(0)), SMatrix{3,3,Complex{T},9}(I))
    Chipman{T}(evec::SVector{3,Complex{T}},pmatrix::SMatrix{3,3,Complex{T},9}) where{T<:Real}= new{T}(evec,pmatrix)
end
export Chipman

pmatrix(a::Chipman) = a.pmatrix
export pmatrix
electricfieldvector(a::Chipman) = a.electricfieldvector
export electricfieldvector

struct NoPolarization{T<:Real} <: AbstractPolarization{T}
end

"""Computes the O matrix (Polarized Light and Optical Systems by Chimpan et al, chapter 9) which represents the s-p coordinate frame for an interface. In OpticSim surface normals are defined positive pointing outward from a solid. In Chipman the surface normal is defined positive on the exiting side of an interface. Need to take this account when passing `surfacenormal` argument."""
function localtoworld(surfacenormal::AbstractVector{T},incidentvector::AbstractVector{T}) where{T<:Real}
    k̂ = normalize(incidentvector)
    ŝ = normalize(k̂ × surfacenormal)
    p̂ = incidentvector × ŝ
    return SMatrix{3,3,T,9}(ŝ[1],ŝ[2],ŝ[3],p̂[1],p̂[2],p̂[3],k̂[1],k̂[2],k̂[3])
end
export localtoworld

worldtolocal(surfacenormal::AbstractVector{T},incidentvector::AbstractVector{T}) where{T<:Real} = localtoworld(surfacenormal,incidentvector)'
export worldtolocal

"""Create jones matrix for dielectric interfaces"""
jonesmatrix(s::Complex{T},p::Complex{T}) where{T<:Real} = SMatrix{3,3,Complex{T},9}(s,0,0,0,p,0,0,0,1)
jonesmatrix(s::T,p::Complex{T}) where{T<:Real} = jonesmatrix(Complex{T}(s),p)
jonesmatrix(s::Complex{T},p::T) where{T<:Real} = jonesmatrix(s,Complex(p))
jonesmatrix(s::T,p::T) where{T<:Real} = jonesmatrix(Complex(s),Complex(p))
export jonesmatrix

""" specialized function for no polarization case"""
function composepolarization(s::Complex{T},p::Complex{T}, Tₐ::T, normal::SVector{N,T}, incidentray::SVector{3,T},exitray::SVector{3,T},polarizationinput::P) where{T<:Real,N,P<:Polarization.NoPolarization{T}}
    return NoPolarization{T}()
end

""" computes the polarization information for an interface given an input polarization of type Chipman"""
function composepolarization(s::Complex{T},p::Complex{T}, Tₐ::T, normal::SVector{N,T}, incidentray::SVector{3,T},exitray::SVector{3,T},polarizationinput::P) where{T<:Real,N,P<:Polarization.Chipman{T}}
    pinput = Polarization.pmatrix(polarizationinput)
    evector = Polarization.electricfieldvector(polarizationinput)

    #compute worldtolocal for incident ray
    plocal = Polarization.worldtolocal(normal,incidentray)*pinput
    temp = Polarization.jonesmatrix(s,p) * plocal
    poutput = Polarization.localtoworld(normal,exitray)*temp
    evectorout = poutput*evector*sqrt(Tₐ)
    return Polarization.Chipman{T}(evectorout,poutput)
end
export composepolarization


end #module
export Polarization