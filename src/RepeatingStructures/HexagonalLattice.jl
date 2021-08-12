# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

const sin60 = .5*sqrt(3)
const cos60 = .5
"""coordinates of a unit hexagon tile polygon with unit length sides centered about the point (0,0)"""
const hexcoords = [
		1 0;
		cos60 -sin60;
		-cos60 -sin60;
		-1 0;
		-cos60 sin60;
		cos60 sin60
		]

struct HexBasis1{N,T} <: Basis{N,T}
    HexBasis1(::Type{T} = Float64) where{T<:Real} = new{2,T}()
end
export HexBasis1

basis(::HexBasis1{2,T}) where{T} = SMatrix{2,2,T}(T(1.5),T(.5)*sqrt(T(3)),T(1.5),T(-.5)*(sqrt(T(3))))

# SVector{2,SVector{2,T}}(hexe₁(T),hexe₂(T))

"""Returns the vertices of the unit tile polygon for the basis"""
function tilevertices(::HexBasis1{2,T}) where{T}
sin60 = T(.5)*sqrt(T(3))
cos60 = T(.5)
return T.([
		1 0;
		cos60 -sin60;
		-cos60 -sin60;
		-1 0;
		-cos60 sin60;
		cos60 sin60
		])
end

# coordinate offsets to move up, down, etc., numsteps hex cells in a hex lattice defined in HexBasis1
hexup(numsteps = 1) = numsteps .* (1,-1)
hexdown(numsteps = 1) = numsteps .* (-1,1)
hexupright(numsteps = 1) = numsteps .* (1,0)
hexdownleft(numsteps = 1) = numsteps .* (-1,0)
hexdownright(numsteps = 1) = numsteps .* (0,1)
hexupleft(numsteps = 1) = numsteps .* (0,-1)

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

"""Returns ```centerpoint``` plus all its neighbors in a ring of size n. This is a hexagonal region of the lattice n times as large as the unit hexagon cell.

Example:

```
Vis.@wrapluxor Vis.drawhexcells(50, hexregion((0,0),2))
```
"""
function region(::Type{HexBasis1},centerpoint::Tuple{Int64,Int64}, n::Int64) 
    f(i) = i==0 ? () : (neighbors(HexBasis1,centerpoint,i)...,f(i-1)...)
    cells = (centerpoint...,f(n)...)
    numcells = length(cells) ÷ 2
    println(numcells)
    # convert to matrix form since that is what the drawing code expects
    temp = MMatrix{2,numcells,Int64}(undef)
    for i in 1:numcells
        println(cells[i])
        temp[1,i] = cells[2*i-1]
        temp[2,i] = cells[2*i]
    end

    return SMatrix{2,numcells,Int64}(temp)
end
export region

"""Returns all hex cells contained in the ring of size n centered around centerpoint. Use region if you want all the cells contained in the ring of size n."""
function neighbors(::Type{HexBasis1}, centerpoint::Tuple{Int64,Int64},n::Int) where{T}
    temp = MMatrix{2,n*6,Int64}(undef)
    hoffsets = ringoffsets(Val{n})
    latticepoint = centerpoint .+ n .* (hexdown())

    #convert hexoffsets to absolute coordinate positions
    for i in 1:size(temp)[2]
        temp[1,i] = latticepoint[1]
        temp[2,i] = latticepoint[2]
        latticepoint = latticepoint .+ hoffsets[i] #last computed value of latticepoint won't be used
    end
    return SMatrix{2,n*6,Int64}(temp)
end
export neighbors

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
    hexbasis2to1(i,j) = [i+j,-j]

    result = Matrix{Int}(undef,2,(2*numi+1)*(2*numj+1))
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
                result[:,(i2+numi)*(2*numj+1) + j2+numj+1] = hexbasis2to1(i2,jtemp)
            end
        end
    end
    return result
end
export Repeat

struct HexBasis3{N,T}<:Basis{N,T}
    HexBasis3(::Type{T} = Float64) where{T} = new{2,T}()
end
export HexBasis3

function tilevertices(::HexBasis3{2,T}) where{T}
    sin60 = T(.5)*sqrt(T(3))
    cos60 = T(.5)
    return T.([
        0 1;
        -sin60 cos60;
        -sin60 -cos60;
        0 -1;
        sin60 -cos60;
        sin60 cos60;
            ])
end

basis(::HexBasis3{2,T}) where{T} = SMatrix{2,2,T}(T(2*sin60),T(0),T(sin60),T(1.5))

