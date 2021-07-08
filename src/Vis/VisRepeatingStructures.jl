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
    Luxor.sethue("black")
    Luxor.poly(hexagon, :stroke, close=true)
    #scale and offset text so coordinates are readable
    Luxor.fontsize(hexsize/3)
    Luxor.text("$i, $j",Luxor.Point(-hexsize/3,hexsize/9))

    # arrowlength = hexsize*.5*sqrt(3)/norm(e₁)
    # Luxor.arrow(Luxor.Point(0.0,0.0),arrowlength*Luxor.Point(e₁...))
    # Luxor.sethue("blue")
    # Luxor.arrow(Luxor.Point(0.0,0.0),arrowlength*Luxor.Point(e₂...))
    # Luxor.sethue("black")
    Luxor.translate(-offset)
end

function drawhexcells(hexsize,cells, color = nothing)
    if color === nothing
        colors = Colors.distinguishable_colors(length(cells),lchoices = range(40,stop=100,length = 15))
        for (i,cell) in pairs(cells)
            drawhex(Repeat.HexBasis1(),hexsize,cell[1],cell[2],colors[i])
        end
    else
        for cell in cells
            drawhex(Repeat.HexBasis1(),hexsize,cell[1],cell[2],color)
        end
    end
end

macro wrapluxor(f)
    return :(Luxor.@draw begin
        $f
    end 1000 1000)
end
export @wrapluxor


function drawhexring() end
