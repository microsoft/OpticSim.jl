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
    Luxor.sethue("grey")
    Luxor.poly(tile, :stroke, close=true)
    Luxor.sethue("black")
    Luxor.circle(Luxor.Point(0,0),2.0,:fill)
    #scale and offset text so coordinates are readable
    Luxor.fontsize(tilesize/3)
    Luxor.text("$i, $j",Luxor.Point(-tilesize/3,tilesize/2.5))
    if name !== nothing
        Luxor.text(name,Luxor.Point(-tilesize/4,-tilesize/4))
    end

    Luxor.translate(-offset)
end

"""Draws a list of hexagonal cells, represented by their lattice coordinates"""
function drawcells(tilebasis::Basis, tilesize,cells; color::Union{AbstractArray,Nothing} = nothing, name::Union{AbstractArray,Nothing} = nothing, format=:png, resolution=(500,500))
    Luxor.Drawing(resolution[1], resolution[2], format)
    Luxor.origin()
    Luxor.background(Colors.RGBA(0, 1, 1, 0.0))
    if color === nothing
        color = Colors.distinguishable_colors(length(cells),lchoices = 30:250) #type unstable but not performance critical code
    end

    for (i,cell) in pairs(cells)
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

function draw(clstr::Repeat.ClusterWithProperties,scale = 50.0)
    points = clustercoordinates(clstr,0,0)
    println(points)
    ptvecs = [points[:,i] for i in 1:size(points)[2]]
    props = Repeat.properties(clstr)
    drawcells(elementbasis(cluster(clstr)),scale,ptvecs,color = props[:,:Color],name = props[:,:Name])
end
    


