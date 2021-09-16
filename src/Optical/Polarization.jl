abstract type AbstractPolarization end

struct Chipman <: AbstractPolarization
    P::SMatrix{3,3,Complex,9} #3D polarization matrix
    Chipman() = new(SMatrix{3,3,Complex,9}(0.0 + 0.0im,0.0 + 0.0im,0.0 + 0.0im,0.0 + 0.0im,0.0 + 0.0im,0.0 + 0.0im,0.0 + 0.0im,0.0 + 0.0im,0.0 + 0.0im))
end