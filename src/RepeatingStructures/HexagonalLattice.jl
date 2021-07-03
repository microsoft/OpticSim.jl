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

# coordinate offsets to move up, down, etc., one hex cell in a hex lattice defined in HexBasis1
hexup() = (1,-1)
hexdown() = (-1,1)
hexupright() = (1,0)
hexdownleft() = (-1,0)
hexdownright() = (0,1)
hexupleft() = (0,-1)
export hexup,hexdown,hexupright,hexdownleft,hexdownright,hexupleft

"""HexBasis2 is defined because certainly conditional tests, such as determing which hex coordinates lie inside a rectangula box, are more easily performed in this basis. For all other purposes use HexBasis1"""
struct HexBasis2{N,T} <: HexagonalLattice{N,T}
    HexBasis2(::Type{T} = Float64) where{T<:Real} = new{2,T}()
end


basis(a::HexBasis2{2,T}) where{T} = SVector{2,SVector{2,T}}(hex2e₁(T),hex2e₂(T))


export ring1offsets,ring2offsets

"""returns i,j coordinates of tiles in the ring"""
function ring(a::HexagonalLattice, startingcoordinates::NTuple{N,T}, ringnum) where{T<:Int,N}
    offsets = ringoffsets[ringnum]()
    ringlength = length(offsets)
    result = MVector{ringlength,NTuple{N,T}}(undef)
    result[1] = startingcoordinates
    for i in 1:ringnum
        result[1] = result[1] .+ hexdown()
    end

    for i in 2:ringlength
        result[i] = result[i-1] .+ offsets[i]
    end
    return SVector{ringlength,NTuple{N,T}}(result)
end
export ring

hexcycle = SVector{6,NTuple{2,Int64}}(hexupright(),hexup(),hexupleft(),hexdownleft(),hexdown(),hexdownright())

ring(n::Int64) = ring(Val{n})
function ring(::Type{Val{N}}) where N
    temp = MVector{N*6,NTuple{2,Int64}}(undef)

    for i in 1:6
        for j in 1:N
            temp[(i-1)*N + j] = hexcycle[i]
        end
    end
    
    return SVector{N*6,NTuple{2,Int64}}(temp)
end
export ring

function allhexcells(latticecoord, n::Int64) 
    f(i) = i==0 ? () : (hexpoints(latticecoord,i)...,f(i-1)...)
    return ((0,0),f(n)...)
end
export allhexcells

# hexoffsets(::Type{Val{1}}) = SVector{6,NTuple{2,Int64}}((hexdown(),hexupright(),hexup(),hexupleft(),hexdownleft(),hexdown()))
# hexoffsets(::Type{Val{2}}) = SVector{12,NTuple{2,Int64}}(hexdown() .+ hexdown(),hexupright(),hexupright(),hexup(),hexup(),hexupleft(),hexupleft(),hexdownleft(),hexdownleft(),hexdown(),hexdown(),hexdownright())

function hexpoints(latticepoint::Tuple{Int64,Int64},n) where{T}
    temp = MVector{n*6,Tuple{Int64,Int64}}(undef)
    hoffsets = ring(Val{n})
    latticepoint = latticepoint .+ n .* (hexdown())
    println(hoffsets)
    println()
    for i in 1:length(temp)
        println(latticepoint)
        temp[i] = latticepoint
        latticepoint = latticepoint .+ hoffsets[i] #last computed value of latticepoint won't be used
    end
    return SVector{n*6,Tuple{Int64,Int64}}(temp)
end

# ring1(latticepoint) = ((0,0),hexpoints(latticepoint,1)...)
# ring2(latticepoint) = (ring1(latticepoint)...,hexpoints(latticepoint,2)...)

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

hexbasis2to1(i,j) = (i+j,-j)

"""computes the basis coordinates of hexagonal cells defined in the Hex1Basis in a rectangular pattern 2*numi+1 by 2*numj+1 wide and high respectively."""
function hexcells(numi,numj)
    #i and j are coordinates in the Hex2Basis. The conditional tests are simpler in this basis. Convert back to
    #Hex1Basis coordinates before returning.
    result = Array{Tuple{Int64,Int64},2,}(undef,2*numi+1,2*numj+1)
    for i2 in -numi:numi
        let offsetj
            if i2 < 0
            offsetj = -abs(div(i2,2))
            else
                offsetj = div(i2+1,2)
            end
            for j2 in -numj:numj
                # conversion from Hex2Basis i2,j2 to Hex1Basis i,j is 
                # (i,j) = (i2,0) + j2*(1,-1)
                jtemp = j2-offsetj
                # i1 = i2 + jtemp
                # j1 = -jtemp
                result[i2+numi+1,j2+numj+1] = hexbasis2to1(i2,jtemp)
            end
        end
    end
    return result
end
