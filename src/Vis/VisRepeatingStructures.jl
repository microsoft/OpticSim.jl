function drawhex(hexsize,i,j,color)
    Luxor.arrow(Luxor.Point(0.0,0.0),hexsize*Luxor.Point(e₁...))
    Luxor.sethue("blue")
    Luxor.arrow(Luxor.Point(0.0,0.0),hexsize*Luxor.Point(e₂...))
    Luxor.sethue("black")
    offset = Luxor.Point(hexsize*(i*e₁ + j*e₂)...)
    Luxor.translate(offset)
    
    Luxor.sethue(color)
    Luxor.poly(hexagon, :fill, close=true)
    Luxor.sethue("black")
    Luxor.poly(hexagon, :stroke, close=true)
    Luxor.text("$i, $j")
    Luxor.translate(-offset)
end

function drawhexcells(hexsize, numi, numj)
    Luxor.@draw begin
        cells = Repeat.hexcells
        
                    drawhex(hexsize,i,j,"yellow")
                end
            end
        end
        (minx,maxx),(miny,maxy) = bbox(numi,numj)
        println("$minx $maxx")
        Luxor.sethue("black")
        
        Luxor.box(Luxor.Point(hexsize*minx,hexsize*miny),Luxor.Point(hexsize*maxx,hexsize*maxy), :stroke)
    end 1000 1000
end