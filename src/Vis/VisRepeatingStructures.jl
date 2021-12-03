# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.



#############################################################################



# lattice visualizations are drawn with Luxor because it is easier to do 2D drawings with Luxor than with Makie.

function draw(tilebasis::AbstractBasis, tilesize, i, j, color, name)
    vertices = Repeat.tilevertices(tilebasis)
    tile = tilesize * [Luxor.Point(vertices[:,i]...) for i in 1:size(vertices)[2]]
    pt = tilesize * tilebasis[i,j]
    offset = Luxor.Point(pt[1], -pt[2]) # flip y so indices show up correctly
    
    Luxor.translate(offset)
    
    Luxor.sethue(color)
   
    Luxor.poly(tile, :fill, close=true)
    Luxor.sethue("lightgrey")
    Luxor.poly(tile, :stroke, close=true)
    Luxor.sethue("black")

    # scale and offset text so coordinates are readable
    Luxor.fontsize(tilesize / 3)
    Luxor.text("$i, $j", Luxor.Point(-tilesize / 3, tilesize / 2.5))
    if name !== nothing
        Luxor.text(name, Luxor.Point(-tilesize / 4, -tilesize / 4))
    end

    Luxor.translate(-offset)
end

function drawcells(clstr::Repeat.ClusterWithProperties,scale,points) 
    _,npts = size(points)
    repeats = npts ÷ clustersize(clstr)
    props = repeat(Repeat.properties(clstr), repeats)
    drawcells(elementbasis(clstr), scale, points, color=props[:,:Color], name=props[:,:Name])
end

drawcells(clstr::Repeat.LatticeCluster,scale,points) = drawcells(elementbasis(clstr),scale,points)

"""Draws a list of hexagonal cells, represented by their lattice coordinates, which are represented as a 2D matrix, with each column being one lattice coordinate."""
function drawcells(tilebasis::AbstractBasis, tilesize, cells::AbstractMatrix; color::Union{AbstractArray,Nothing}=nothing, name::Union{AbstractArray{String},Nothing}=nothing, format=:png, resolution=(1000, 1000))

    numcells = size(cells)[2]
    Luxor.Drawing(resolution[1], resolution[2], format)
    Luxor.origin()
    Luxor.background(Colors.RGBA(0, 1, 1, 0.0))
    if color === nothing
        color = Colors.distinguishable_colors(numcells, lchoices=30:250) # type unstable but not performance critical code
    end

    for i in 1:numcells
        cell = cells[:,i]
        cellname = name === nothing ? nothing : name[i]
        draw(tilebasis, tilesize, cell[1], cell[2], color[i], cellname)
    end
 
    Luxor.finish()
    if (format == :svg)
        res = Luxor.svgstring()
        return res
    end
    if (format == :png)
        Luxor.preview()
    end
end

""" draw the ClusterWithProperties at coordinates specified by lattice_coordinate_offset """
function draw(clstr::Repeat.AbstractLatticeCluster, cluster_coordinate_offset::AbstractMatrix{T} , scale=50.0) where{T}
    dims = size(cluster_coordinate_offset)
    clstrsize = clustersize(clstr)
    points = Matrix(undef, dims[1], dims[2] * clstrsize)
    for i in 1:dims[2]
        points[:,(i - 1) * clstrsize + 1:i * clstrsize] = clustercoordinates(clstr, cluster_coordinate_offset[:,i]...)
    end
   
    drawcells(clstr,scale,points)
end
    


