function timeparaxiallens()
    lens = ParaxialLensRect(10.0,5.0,5.0,[0.0,0.0,-1.0],[0.0,0.0,0.0])
    detector = Rectangle(5.0,5.0,[0.0,0.0,1.0],[0.0,0.0,-10.0],interface = opaqueinterface(Float64))
    sys = CSGOpticalSystem(LensAssembly(lens), detector)
    raygenerator = Emitters.Sources.Source(transform = translation(0.0,0.0,10.0), origins = Emitters.Origins.RectGrid(10.0,10.0,2500,2500),directions = Emitters.Directions.Constant(0.0,0.0,-1.0))
    return ()-> traceMT(sys,raygenerator)
    # Vis.drawtracerays(sys, test = true, trackallrays = true, colorbynhits = true)
end
export timeparaxiallens
