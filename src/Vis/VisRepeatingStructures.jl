# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.



#############################################################################



# lattice visualizations are drawn with Luxor because it is easier to do 2D drawings with Luxor than with Makie.

function drawhex(hexbasis::Repeat.Basis,hexsize,i,j,color)
    hexagon = hexsize*[Luxor.Point(Repeat.tilevertices(hexbasis)[i,:]...) for i in 1:6]
    pt = hexsize*hexbasis[i,j]
    offset = Luxor.Point(pt[1],-pt[2]) #flip y so indices show up correctly
    Luxor.translate(offset)
    
    Luxor.sethue(color)
   
    Luxor.poly(hexagon, :fill, close=true)
    Luxor.sethue("grey")
    # Luxor.setdash("dash")
    Luxor.poly(hexagon, :stroke, close=true)
    Luxor.sethue("black")
    Luxor.circle(Luxor.Point(0,0),2.0,:fill)
    #scale and offset text so coordinates are readable
    Luxor.fontsize(hexsize/3)
    Luxor.text("$i, $j",Luxor.Point(-hexsize/3,hexsize/2.5))

    # arrowlength = hexsize*.5*sqrt(3)/norm(e₁)
    # Luxor.arrow(Luxor.Point(0.0,0.0),arrowlength*Luxor.Point(e₁...))
    # Luxor.sethue("blue")
    # Luxor.arrow(Luxor.Point(0.0,0.0),arrowlength*Luxor.Point(e₂...))
    # Luxor.sethue("black")
    Luxor.translate(-offset)
end

function drawhexcells(hexsize,cells, color::Union{AbstractArray,Nothing} = nothing; format=:png, resolution=(500,500))
    Luxor.Drawing(resolution[1], resolution[2], format)
    Luxor.origin()
    Luxor.background(Colors.RGBA(0, 1, 1, 0.0))
    if color === nothing
        distcolors = Colors.distinguishable_colors(length(cells),lchoices = range(40,stop=100,length = 15))
            for (i,cell) in pairs(cells)
            drawhex(Repeat.HexBasis1(),hexsize,cell[1],cell[2],distcolors[i])
        end
    else
        for (i,cell) in pairs(cells)
            drawhex(Repeat.HexBasis1(),hexsize,cell[1],cell[2],color[i])
        end
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

function draw(lattice::Repeat.Basis, scale = 50.0)
    Luxor.sethue("black")
    pt = scale*lattice[i,j]
    offset = Luxor.Point(pt[1],-pt[2]) #flip y so indices show up correctly
    luxor.circle(offset,scale*.1,:fill)
     #scale and offset text so coordinates are readable
     Luxor.fontsize(scale/3)
     Luxor.text("$i, $j",Luxor.Point(-scale/3,scale/9))
end
export draw

