abstract type HexagonalLattice{N,T} <: Basis{N,T} end

hexe₁(::Type{T}=Float64) where{T<:Real} = SVector{2,T}(T(1.5),T(.5)*sqrt(T(3)))
hexe₂(::Type{T}=Float64) where{T<:Real} = SVector{2,T}(T(1.5),T(-.5)*(sqrt(T(3))))

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

""" (i,j) offsets of a ring represented in the HexBasis1 lattice basis. The hexagonal grid is assumed to start from the bottom center hex cell. To generate larger rings repeat this pattern
 ring 1 = (1, 0)      (1, -1)       (0, -1)       (-1, 0)       (-1, 1)       (0, 1)
 ring 2 = (1, 0)(1, 0)(1, -1)(1, -1)(0, -1)(0, -1)(-1, 0)(-1, 0)(-1, 1)(-1, 1)(0, 1)(0, 1)
"""
const hexcycle = SVector{6,NTuple{2,Int64}}(hexupright(),hexup(),hexupleft(),hexdownleft(),hexdown(),hexdownright())


"""returns (i,j) offsets of ring n defined in the HexBasis1 lattice basis."""
ringoffsets(n::Int64) = ringoffsets(Val{n})
function ringoffsets(::Type{Val{N}}) where N
    temp = MVector{N*6,NTuple{2,Int64}}(undef)

    for i in 1:6
        for j in 1:N
            temp[(i-1)*N + j] = hexcycle[i]
        end
    end
    
    return SVector{N*6,NTuple{2,Int64}}(temp)
end
export ringoffsets

"""Returns the positions of all the cells in ring n centered at ```centerpoint``` plus all cells enclosed by this ring. This is a hexagonal region of the lattice n times as large as the unit hexagon cell.

Example:

```
Vis.@wrapluxor Vis.drawhexcells(50,filledhexagon((0,0),2))
```
"""
function cellsenclosedbyring(centerpoint::Tuple{Int64,Int64}, n::Int64) 
    f(i) = i==0 ? () : (ringcells(centerpoint,i)...,f(i-1)...)
    return (centerpoint,f(n)...)
end
export cellsenclosedbyring

"""Returns all hex cells contained in the ring of size n centered around centerpoint"""
function ringcells(centerpoint::Tuple{Int64,Int64},n) where{T}
    temp = MVector{n*6,Tuple{Int64,Int64}}(undef)
    hoffsets = ringoffsets(Val{n})
    latticepoint = centerpoint .+ n .* (hexdown())

    #convert hexoffsets to absolute coordinate positions
    for i in 1:length(temp)
        temp[i] = latticepoint
        latticepoint = latticepoint .+ hoffsets[i] #last computed value of latticepoint won't be used
    end
    return SVector{n*6,Tuple{Int64,Int64}}(temp)
end
export ringcells

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

"""computes bounding box of hex tiles arranged in an approximately rectangular layout, with 2*numj+1 hexagons in a row, and 2*numi+1 hexagons in a column"""
function bbox(numi,numj)
    return (xbounds(numi),ybounds(numj))
    
end
export bbox


"""computes the basis coordinates of hexagonal cells defined in the Hex1Basis in a rectangular pattern 2*numi+1 by 2*numj+1 wide and high respectively."""
function hexcellsinbox(numi,numj)
    #i2 and j2 are coordinates in the basis (1.5,.5√3),(0,-√3). The conditional tests are simpler in this basis. Convert to Hex1Basis1 coordinates before returning.
    hexbasis2to1(i,j) = (i+j,-j)

    result = Array{Tuple{Int64,Int64},1}(undef,(2*numi+1)*(2*numj+1))
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
                result[(i2+numi)*(2*numj+1) + j2+numj+1] = hexbasis2to1(i2,jtemp)
            end
        end
    end
    return result
end
export Repeat
