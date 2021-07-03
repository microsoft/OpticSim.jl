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


function drawhex(hexbasis::Repeat.Basis,hexsize,i,j,color)
    hexagon = hexsize*[Luxor.Point(hexcoords[i,:]...) for i in 1:6]
   
    offset = Luxor.Point(hexsize*(hexbasis[i,j])...)
    Luxor.translate(offset)
    Luxor.sethue(color)
    Luxor.poly(hexagon, :fill, close=true)
    Luxor.sethue("black")
    Luxor.poly(hexagon, :stroke, close=true)
    Luxor.text("$i, $j")

    # arrowlength = hexsize*.5*sqrt(3)/norm(e₁)
    # Luxor.arrow(Luxor.Point(0.0,0.0),arrowlength*Luxor.Point(e₁...))
    # Luxor.sethue("blue")
    # Luxor.arrow(Luxor.Point(0.0,0.0),arrowlength*Luxor.Point(e₂...))
    # Luxor.sethue("black")
    Luxor.translate(-offset)
end

function drawhexcells(hexsize,cells)
    colors = Colors.distinguishable_colors(length(cells))

    for (cell, color) in zip(cells,colors)
        
        drawhex(Repeat.HexBasis1(),hexsize,cell[1],cell[2],color)
    end
end

function drawhexcellsinbox(hexsize, numi, numj)
        cells = Repeat.hexcells(numi,numj)
        for cell in cells
            drawhex(Repeat.HexBasis1(), hexsize,cell[1],cell[2],"yellow")
        end
         
        (minx,maxx),(miny,maxy) = bbox(numi,numj)
        Luxor.sethue("black")
        
        Luxor.box(Luxor.Point(hexsize*minx,hexsize*miny),Luxor.Point(hexsize*maxx,hexsize*maxy), :stroke)
end

macro wrapluxor(f...)
    :(Luxor.@draw begin
        $(f...)
    end 1000 1000)
end
export @wrapluxor

function drawhexring() end
