function timeparaxiallens()
    lens = ParaxialLensRect(10.0,5.0,5.0,[0.0,0.0,-1.0],[0.0,0.0,0.0])
    detector = Rectangle(5.0,5.0,[0.0,0.0,1.0],[0.0,0.0,-10.0],interface = opaqueinterface(Float64))
    sys = CSGOpticalSystem(LensAssembly(lens), detector)
    # raygenerator = Emitters.Sources.Source(transform = translation(0.0,0.0,10.0), origins = Emitters.Origins.RectGrid(10.0,10.0,250,250),directions = Emitters.Directions.Constant(0.0,0.0,-1.0)) 
    
    raygenerator = Emitters.Sources.Source(transform = translation(0.0,0.0,10.0), origins = Emitters.Origins.Point(10.0,10.0,10.0),directions = Emitters.Directions.UniformCone(π/4,62500))
    # ray =  raygenerator[1]
    # return ()->  traceMT(sys,raygenerator,printprog=false)
    println("numrays $(length(raygenerator))")

    return ()-> for ray in raygenerator end
    # for i in 1:2000 
        # trace(sys,ray)
     # end
    # Vis.drawtracerays(sys, test = true, trackallrays = true, colorbynhits = true)
end
export timeparaxiallens
