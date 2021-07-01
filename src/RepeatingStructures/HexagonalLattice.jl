abstract type HexagonalLattice{N,T} <: Basis{N,T} end

const sin60 = .5*sqrt(3)
const cos60 = .5

const hexcoords = [
    1 0;
    cos60 -sin60;
    -cos60 -sin60;
    -1 0;
    -cos60 sin60;
    cos60 sin60
    ]


hexe₁(::Type{T}=Float64) where{T<:Real} = SVector{2,T}(T(1.5),T(.5)*sqrt(T(3)))
hexe₂(::Type{T}=Float64) where{T<:Real} = SVector{2,T}(T(1.5),T(-.5)*(sqrt(T(3))))


#Alternative basis which makes it a bit easier to figure out how coordinates correspond with movements in the hex lattice
const hex2e₁(::Type{T}=Float64) where{T<:Real} = SVector{2,T}(T(-1.5),T(.5)*sqrt(T(3)))
const hex2e₂(::Type{T}=Float64) where{T<:Real} = SVector{2,T}(T(0),sqrt(T(3)))

struct HexBasis1{N,T} <: HexagonalLattice{N,T}
    HexBasis1(::Type{T} = Float64) where{T<:Real} = new{2,T}()
end
export HexBasis1

basis(a::HexBasis1{2,T}) where{T} = SVector{2,SVector{2,T}}(hexe₁(T),hexe₂(T))

# coordinate offsets to move up, down, etc., one hex cell in a hex lattice defined in HexBasis2
hexup(a::HexBasis1) = (1,-1)
hexdown(a::HexBasis1) = (-1,1)
hexupright(a::HexBasis1) = (1,0)
hexdownleft(a::HexBasis1) = (-1,0)
hexdownright(a::HexBasis1) = (0,1)
hexupleft(a::HexBasis1) = (0,-1)
export hexup,hexdown,hexupright,hexdownleft,hexdownright,hexupleft

struct HexBasis2{N,T} <: HexagonalLattice{N,T}
    HexBasis2(::Type{T} = Float64) where{T<:Real} = new{2,T}()
end
export HexBasis2

basis(a::HexBasis2{2,T}) where{T} = SVector{2,SVector{2,T}}(hex2e₁(T),hex2e₂(T))

hexagonallattice(pitch::T = 1.0) where{T<:Real} = LatticeBasis(pitch*SVector{2,T}(T(1.5),T(.5)*sqrt(T(3))),pitch*SVector{2,T}(T(1.5),T(-.5)*sqrt(T((3)))))
export hexagonallattice

#offset constants used for generating 1 and 2 hexagonal rings. 
const ringoffsets = (
    (a::HexagonalLattice) -> (hexdown(a),hexupright(a),hexup(a),hexupleft(a),hexdownleft(a),hexdown(a)),
    (a::HexagonalLattice) -> (hexdown(a),hexupright(a),hexupright(a),hexup(a),hexup(a),hexupleft(a),hexupleft(a),hexdownleft(a),hexdownleft(a),hexdown(a),hexdown(a),hexdownright(a))
)


export ring1offsets,ring2offsets

"""returns i,j coordinates of tiles in the ring"""
function ring(a::HexagonalLattice, startingcoordinates::NTuple{N,T}, ringnum) where{T<:Int,N}
    offsets = ringoffsets[ringnum](a)
    ringlength = length(offsets)
    result = MVector{ringlength,NTuple{N,T}}(undef)
    result[1] = startingcoordinates
    for i in 1:ringnum
        result[1] = result[1] .+ hexdown(a)
    end

    for i in 2:ringlength
        result[i] = result[i-1] .+ offsets[i]
    end
    return SVector{ringlength,NTuple{N,T}}(result)
end
export ring

""" returns lattices positions rather than lattice indices in a hex7 pattern about latticepoint"""
hex7points(latticepoint::SVector{N,T},basis::HexagonalLattice) where{N,T<:Real} = map(offset -> basis[offset[1],offset[2]] + latticepoint,hex7basiscoordinates(basis))
export hex7points

"""returns lattice offsets (i,j) in a hex13 pattern about latticepoint"""
hex13() = SVector{13,NTuple{2,Int64}}((0,0),hexring(2)...)
export hex13

""" returns lattice positions rather than lattice indices in a hex13 pattern about latticepoint"""
hex13points(latticepoint::SVector{N,T}) where{N,T<:Real} = map(offset -> hexbasis1[offset[1],offset[2]] + latticepoint,hex13())
export hex13points
	
	

# const hexagon = hexsize*[Luxor.Point(hexcoords[i,:]...) for i in 1:6]

# function drawhex(i,j,color)
#     Luxor.arrow(Luxor.Point(0.0,0.0),hexsize*Luxor.Point(e₁...))
#     Luxor.sethue("blue")
#     Luxor.arrow(Luxor.Point(0.0,0.0),hexsize*Luxor.Point(e₂...))
#     Luxor.sethue("black")
#     offset = Luxor.Point(hexsize*(i*e₁ + j*e₂)...)
#     Luxor.translate(offset)
    
#     Luxor.sethue(color)
#     Luxor.poly(hexagon, :fill, close=true)
#     Luxor.sethue("black")
#     Luxor.poly(hexagon, :stroke, close=true)
#     Luxor.text("$i, $j")
#     Luxor.translate(-offset)
# end

function xbounds(numi)
    numevens = div(numi,2)
    numodds = numevens + mod(numi,2)
    basewidth = numodds + 2*numevens + 1
    
    if iseven(numi)
        maxx = basewidth
    else
        maxx = basewidth + sin(deg2rad(30))
    end
    
    minx  = -maxx
    return (minx,maxx)
end

function ybounds(numj)
    maxy = 2*sin60*(numj + 1)
    return (-maxy,maxy-sin60)
end

"""only works for [-1.5,sin60],[0.0],2.0*sin60] basis"""
function bbox(numi,numj)
    return (xbounds(numi),ybounds(numj))
    
end

"""computes the basis coordinates of hexagonal cells in a rectangular pattern 2*numi+1 by 2*numj+1 wide and high respectively."""
function hexcells(numi,numj)
    result = Array{Tuple{Int64,Int64},2,}(undef,2*numi+1,2*numj+1)
    for i in -numi:numi
        let offsetj
            if i < 0
            offsetj = -abs(div(i,2))
            else
                offsetj = div(i+1,2)
            end
            for j in -numj:numj
                result[i+numi+1,j+numj+1] = (i, j - offsetj)
            end
        end
    end
    return result
end
