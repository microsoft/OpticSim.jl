# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.



#############################################################################



# lattice visualizations are drawn with Luxor because it is easier to do 2D drawings with Luxor than with Makie.

function draw(tilebasis::Basis,tilesize,i,j,color,name)
    tile = tilesize*[Luxor.Point(Repeat.tilevertices(tilebasis)[i,:]...) for i in 1:6]
    pt = tilesize*tilebasis[i,j]
    offset = Luxor.Point(pt[1],-pt[2]) #flip y so indices show up correctly
    
    Luxor.translate(offset)
    
    Luxor.sethue(color)
   
    Luxor.poly(tile, :fill, close=true)
    Luxor.sethue("lightgrey")
    Luxor.poly(tile, :stroke, close=true)
    Luxor.sethue("black")

    #scale and offset text so coordinates are readable
    Luxor.fontsize(tilesize/3)
    Luxor.text("$i, $j",Luxor.Point(-tilesize/3,tilesize/2.5))
    if name !== nothing
        Luxor.text(name,Luxor.Point(-tilesize/4,-tilesize/4))
    end

    Luxor.translate(-offset)
end

"""Draws a list of hexagonal cells, represented by their lattice coordinates"""
function drawcells(tilebasis::Basis, tilesize,cells; color::Union{AbstractArray,Nothing} = nothing, name::Union{AbstractArray{String},Nothing} = nothing, format=:png, resolution=(500,500))

    Luxor.Drawing(resolution[1], resolution[2], format)
    Luxor.origin()
    Luxor.background(Colors.RGBA(0, 1, 1, 0.0))
    if color === nothing
        color = Colors.distinguishable_colors(length(cells),lchoices = 30:250) #type unstable but not performance critical code
    end

    for i in 1:size(cells)[2]
        cell = cells[:,i]
        cellname = name === nothing ? nothing : name[i]
        draw(tilebasis,tilesize,cell[1],cell[2],color[i],cellname)
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

""" draw the LatticeCluster offset to (0,0) """
draw(clstr::LatticeCluster,tilesize = 50.0, cells = hcat([0,0])) = drawcells(clstr.elementbasis,tilesize,clustercoordinates(clstr,0,0))

""" draw the ClusterWithProperties at coordinates specified by cells """
function draw(clstr::Repeat.ClusterWithProperties, cells = hcat([0,0]), scale = 50.0) #this is how you have to make a 2x1 matrix. One would expect [0;0] to work but it doesn't.
    dims = size(cells)
    @assert dims[1] == 2 "dims[1] must be 2 instead was $(dims[1])"
    clstrsize = clustersize(clstr)
    points = Matrix(undef,dims[1],dims[2]*clstrsize)
    for i in 1:dims[2]
        points[:,(i-1)*clstrsize+1:i*clstrsize] = clustercoordinates(clstr,cells[1,i],cells[2,i])
    end
    props = repeat(Repeat.properties(clstr),dims[2])
    drawcells(elementbasis(cluster(clstr)),scale,points,color = props[:,:Color],name = props[:,:Name])
end
    


