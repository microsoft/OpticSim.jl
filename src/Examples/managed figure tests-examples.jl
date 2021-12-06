
using OpticSim
using OpticSim.GlassCat
using OpticSim.Geometry
using OpticSim.Vis
using OpticSim.Emitters
using OpticSim.Examples
using GLMakie

using DataFrames
using LinearAlgebra

using StaticArrays
using FileIO

##
Vis.activateFigureView(true)

lens = SphericalLens(OpticSim.GlassCat.SCHOTT.N_BK7, 0.0, 10.0, 10.0, 5.0, 5.0) 
Vis.draw(lens)
##
topsurf = Plane(
    SVector(0.0, 0.0, 1.0),
    SVector(0.0, 0.0, 0.0),
    interface = FresnelInterface{Float64}(SCHOTT.N_BK7, Air),
    vishalfsizeu = 12.00075,
    vishalfsizev = 12.00075)
botsurf = AcceleratedParametricSurface(ZernikeSurface(
    12.00075,
    radius = -1.14659768e+4,
    aspherics = [(4, 3.68090959e-7), (6, 2.73643352e-11), (7, 1.0e-18), (8, 3.20036892e-14)]),
    17,
    interface = FresnelInterface{Float64}(SCHOTT.N_BK7, Air))
coverlens = Cylinder(12.00075, 1.4) ∩ topsurf ∩ leaf(botsurf, Transform(rotmatd(0, 180, 0), Makie.Vec3(0.0, 0.0, -0.65)))
Vis.draw(coverlens)

 # big mirror with a hole in it
 frontsurfacereflectance = 1.0
 bigmirror = (
     ConicLens(SCHOTT.N_BK7, -72.65, -95.2773500000134, 0.077235, Inf, 0.0, 0.2, 12.18263; frontsurfacereflectance) -
     leaf(Cylinder(4.0, 0.3, interface = opaqueinterface()), Geometry.translation(0.0, 0.0, -72.75))
 )
Vis.draw!(bigmirror)
##

g1, g2 = SCHOTT.N_SK16, GlassCat.SCHOTT.N_SF2
Vis.draw(AxisymmetricOpticalSystem{Float64}(DataFrame(
        SurfaceType  = ["Object", "Standard", "Standard", "Standard", "Stop", "Standard", "Standard", "Image"],
        Radius       = [Inf,      26.777,     66.604,     -35.571,    35.571, 35.571,     -26.777,    Inf    ],
        Thickness    = [Inf,      4.0,        2.0,        4.0,        2.0,    4.0,        44.748,     missing],
        Material     = [Air,      g1,         Air,        g2,         Air,    g1,         Air,        missing],
        SemiDiameter = [Inf,      8.580,      7.513,      7.054,      6.033,  7.003,      7.506,      15.0   ],
    )))

##

function draw_cooketriplet(filename::Union{Nothing,AbstractString} = nothing)
    g1, g2 = SCHOTT.N_SK16, GlassCat.SCHOTT.N_SF2
    sys = AxisymmetricOpticalSystem{Float64}(DataFrame(
        SurfaceType  = ["Object", "Standard", "Standard", "Standard", "Stop", "Standard", "Standard", "Image"],
        Radius       = [Inf,      26.777,     66.604,     -35.571,    35.571, 35.571,     -26.777,    Inf    ],
        Thickness    = [Inf,      4.0,        2.0,        4.0,        2.0,    4.0,        44.748,     missing],
        Material     = [Air,      g1,         Air,        g2,         Air,    g1,         Air,        missing],
        SemiDiameter = [Inf,      8.580,      7.513,      7.054,      6.033,  7.003,      7.506,      15.0   ],
    ))

    origins = Origins.Hexapolar(8, 15.0, 15.0)
    directions = Directions.Constant(0.0, 0.0, -1.0)
    s1 = Sources.Source(; origins, directions, sourcenum=1)

    transform = Transform(rotmatd(10, 0, 0), unitZ3())
    s2 = Sources.Source(; transform, origins, directions, sourcenum=2)

    raygenerator = Sources.CompositeSource(Transform(), [s1, s2])

    trackallrays = test = colorbysourcenum = true; resolution = (1000, 700)
    Vis.drawtracerays(sys; raygenerator, trackallrays, test, colorbysourcenum, resolution)
    return sys
end

sys = draw_cooketriplet()

sys

##

function draw_zoomlenses(filenames::Vector{<:Union{Nothing,AbstractString}} = repeat([nothing], 3))
    stops = [2.89, 3.99, 4.90]
    zooms = [9.48, 4.48, 2.00]
    dists = [4.46970613, 21.21, 43.81]

    transform = Geometry.translation(0.0, 0.0, 10.0)
    origins = Origins.Hexapolar(8, 10.0, 10.0)
    directions = Directions.Constant(0.0, 0.0, -1.0)
    raygenerator = Sources.Source(; transform, origins, directions)

    aspherics = [
        ["4" => 1.0386E-04, "6" => 1.4209E-07, "8" => -8.8495E-09, "10" => 1.2477E-10, "12" => -1.0367E-12, "14" => 3.6556E-15],
        ["4" => 4.2721E-05, "6" => 1.2484E-07, "8" => 9.7079E-09, "10" => -1.8444E-10, "12" => 1.8644E-12, "14" => -7.7975E-15],
        ["4" => 1.1339E-04, "6" => 4.8165E-07, "8" => 1.8778E-08, "10" => -5.7571E-10, "12" => 8.9994E-12, "14" => -4.6768E-14],
    ]
    syss = [
        AxisymmetricOpticalSystem{Float64}(DataFrame(
            SurfaceType = ["Object", "Stop", "Standard", "Standard", "Standard", "Aspheric", "Standard", "Standard", "Aspheric", "Aspheric", "Standard", "Standard", "Standard", "Standard", "Standard", "Standard", "Image"],
            Radius = [Inf64, Inf64, -1.6202203499676E+01, -4.8875855327468E+01, 1.5666614444619E+01, -4.2955326460481E+01, 1.0869565217391E+02, 2.3623907394283E+01, -1.6059097478722E+01, -4.2553191489362E+02, -3.5435861091425E+01, -1.4146272457208E+01, -2.5125628140704E+02, -2.2502250225023E+01, -1.0583130489999E+01, -4.4444444444444E+01, Inf64],
            Parameters = [missing, missing, missing, missing, missing, aspherics[1], missing, missing, aspherics[2], aspherics[3], missing, missing, missing, missing, missing, missing, missing],
            Thickness = [Inf64, 0.0, 5.18, 0.10, 4.40, 0.16, 1.0, 4.96, zoom, 4.04, 1.35, 1.0, 2.80, 3.0, 1.22, dist, missing],
            Material = [Air, Air, OHARA.S_LAH66, Air, NIKON.LLF6, Air, OHARA.S_TIH6, OHARA.S_FSL5, Air, OHARA.S_FSL5, Air, OHARA.S_LAL8, OHARA.S_FSL5, Air, OHARA.S_LAH66, Air, missing],
            SemiDiameter = [Inf64, stop, 3.85433218451, 3.85433218451, 4.36304692871, 4.36304692871, 4.72505505439, 4.72505505439, 4.72505505439, 4.45240784026, 4.45240784026, 4.50974054117, 4.50974054117, 4.50974054117, 4.76271114409, 4.76271114409, 15.0]))
        for (stop, zoom, dist) in zip(stops, zooms, dists)]

    for (sys, filename) in zip(syss, filenames)
        fV = FigureView(figureName = filename)
        Vis.drawtracerays!(sys; raygenerator, trackallrays=true, test=true, numdivisions=50)
    end
    return syss
end

draw_zoomlenses()

##

function draw_schmidtcassegraintelescope(filename::Union{Nothing,AbstractString} = nothing)
    # glass entrance lens on telescope
    topsurf = Plane(
        SVector(0.0, 0.0, 1.0),
        SVector(0.0, 0.0, 0.0),
        interface = FresnelInterface{Float64}(SCHOTT.N_BK7, Air),
        vishalfsizeu = 12.00075,
        vishalfsizev = 12.00075)
    botsurf = AcceleratedParametricSurface(ZernikeSurface(
        12.00075,
        radius = -1.14659768e+4,
        aspherics = [(4, 3.68090959e-7), (6, 2.73643352e-11), (7, 1.0e-18), (8, 3.20036892e-14)]),
        17,
        interface = FresnelInterface{Float64}(SCHOTT.N_BK7, Air))
    coverlens = Cylinder(12.00075, 1.4) ∩ topsurf ∩ leaf(botsurf, Transform(rotmatd(0, 180, 0), Vec3(0.0, 0.0, -0.65)))

    # big mirror with a hole in it
    frontsurfacereflectance = 1.0
    bigmirror = (
        ConicLens(SCHOTT.N_BK7, -72.65, -95.2773500000134, 0.077235, Inf, 0.0, 0.2, 12.18263; frontsurfacereflectance) -
        leaf(Cylinder(4.0, 0.3, interface = opaqueinterface()), translation(0.0, 0.0, -72.75))
    )

    # small mirror supported on a spider
    backsurfacereflectance = 1.0
    smallmirror = SphericalLens(SCHOTT.N_BK7, -40.65, Inf, -49.6845, 1.13365, 4.3223859; backsurfacereflectance)

    obscuration1 = Circle(4.5, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, -40.649), interface = opaqueinterface())
    obscurations2 = Spider(3, 0.5, 12.0, SVector(0.0, 0.0, -40.65))

    # put it together with the detector
    la = LensAssembly(coverlens(), bigmirror(), smallmirror(), obscuration1, obscurations2...)
    det = Circle(3.0, SVector(0.0, 0.0, 1.0), SVector(0.0, 0.0, -92.4542988), interface = opaqueinterface())
    sys = CSGOpticalSystem(la, det)

    # define ray generator
    transform = translation(0.0, 0.0, 10.0)
    origins = Origins.Hexapolar(8, 20.0, 20.0)
    directions = Directions.Constant(0.0, 0.0, -1.0)
    raygenerator = Sources.Source(; transform, origins, directions)

    # draw and output
    Vis.drawtracerays(sys; raygenerator, trackallrays = true, colorbynhits = true, test = true, numdivisions = 100, drawgen = false)
    return nothing
end

draw_schmidtcassegraintelescope()

##

function draw_stackedbeamsplitters(filenames::Vector{<:Union{Nothing,AbstractString}} = repeat([nothing], 3))
    # ReflectOrTransmit: nondeterministic
    # Transmit: deterministic, all beamsplitters transmissive
    # Reflect: deterministic, all beamsplitters reflective
    interfacemodes = [ReflectOrTransmit, Transmit, Reflect]

    for (interfacemode, filename) in zip(interfacemodes, filenames)
        interface = FresnelInterface{Float64}(SCHOTT.N_BK7, Air; reflectance=0.5, transmission=0.5, interfacemode)
        bs_1 = OpticSim.transform(
                Cuboid(10.0, 20.0, 2.0, interface=interface),
                translation(0.0, 0.0, -30.0-2*sqrt(2))*rotationX(π/4))

        l1 = OpticSim.transform(
            SphericalLens(SCHOTT.N_BK7, -70.0, 30.0, Inf, 5.0, 10.0),
            translation(0.0, -1.34, 0.0))

        bs_2 = OpticSim.transform(
            Cuboid(10.0, 20.0, 2.0, interface=interface),
            translation(0.0, 40.0, -30.0+2*sqrt(2))*rotationX(π/4))
            
        l2 = OpticSim.transform(
            SphericalLens(SCHOTT.N_BK7, -70.0, 30.0, Inf, 5.0, 10.0),
            translation(0.0, 40.0, 0.0))

        la = LensAssembly(bs_1(), l1(), bs_2(), l2())
        
        detector = Rectangle(20.0, 40.0, SVector(0.0, 0.0, 1.0), SVector(0.0, 20.0, -130.0); interface = opaqueinterface())
        sys = CSGOpticalSystem(la, detector)

        Vis.drawtracerays(sys; trackallrays=true, rayfilter=nothing, colorbynhits=true)
    end
    return nothing
end

draw_stackedbeamsplitters()


##

book = [
    -.240408 0.0        0.0;
    0       .130825    .501818;
    0       -.501818    -.710275
    ]





    # surface, radius, thickness, index
    patent1 = [0 Inf  -0.20   1
    1   1.962   1.19    1.471   
    2   33.398  0.93    1
    3   -2.182   0.75   1.603  
    4   -6.367  0.10    1
    5   5.694   0.89    1.510   
    6   9.192   0.16    1
    7   1.674   0.85    1.510  
    8   1.509   0.70    1
    9   Inf   0.40    1.516   
    10 Inf 0.64   1]

    # after transpose A3, A4, A5, A6, A7, A8, A9, A10
    patent2 = transpose([
        -1.895E-02	-4.966E-03	-4.388E-02 -1.131E-01	-7.876E-02	9.694E-03	7.429E-02	1.767E-03
        2.426E-02	-1.434E-02	-2.555E-02 -7.863E-02 7.020E-02	-2.516E-03	-6.933E-02 	-4.652E-02
        -5.123E-02	-6.139E-03	5.160E-02	1.094E-01	1.575E-03	-3.606E-03	-5.811E-03	1.625E-02
        8.371E-04	-9.284E-05 -4.307E-02	6.228E-03	-9.958E-03 -2.497E-04 2.396E-03 -3.522E-03
        7.850E-03	6.438E-03	-2.831E-02	-2.916E-02	-7.322E-03	-6.840E-04	2.100E-03 -7.106E-04
        4.091E-03	-5.720E-03	3.162E-02	-5.890E-03	6.914E-04 -1.414E-04 -3.119E-04 3.825E-04
        -7.732E-03	-2.385E-02	4.630E-02	4.123E-03	2.540E-03 2.932E-04	-5.552E-05 6.271E-05
    -4.265E-03	1.108E-02	-4.877E-02	1.041E-03  -7.650E-04 -7.284E-05 7.969E-06 -2.631E-05
    ])
x,y = size(patent2)
aspherics = [[string(j+2)=>patent2[i,j] for j in 1:y] for i in 1:x]


    K=Array([	2.153E+00	4.018E+01	2.105E+00	3.382E+00	-2.211E+02	9.331E-01 -7.617E+00	-2.707E+00][:])

    names =["stop", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "image"]
    apers = [0.5 * 5.57/2.8, 1.2, 1.2, 1.2, 1.2, 2.1, 2.1, 2.5, 2.5, 2.7, 2.7, 3.]



function draw_cellPhoneCameraLens(filenames::Vector{<:Union{Nothing,AbstractString}} = repeat([nothing], 3))

    # surface, radius, thickness, index
    patentTable1 = [
    -1 Inf64 Inf64 1 
    0 Inf64  -0.20   1
    1   1.962   1.19    1.471   
    2   33.398  0.93    1
    3   -2.182   0.75   1.603  
    4   -6.367  0.10    1
    5   5.694   0.89    1.510   
    6   9.192   0.16    1
    7   1.674   0.85    1.510  
    8   1.509   0.70    1
    9   Inf   0.40    1.516   
    10 Inf 0.64   1
    11 Inf missing 1]

    # after transpose A3, A4, A5, A6, A7, A8, A9, A10
    patentTable2 = transpose([
        -1.895E-02	-4.966E-03	-4.388E-02 -1.131E-01	-7.876E-02	9.694E-03	7.429E-02	1.767E-03
        2.426E-02	-1.434E-02	-2.555E-02 -7.863E-02 7.020E-02	-2.516E-03	-6.933E-02 	-4.652E-02
        -5.123E-02	-6.139E-03	5.160E-02	1.094E-01	1.575E-03	-3.606E-03	-5.811E-03	1.625E-02
        8.371E-04	-9.284E-05 -4.307E-02	6.228E-03	-9.958E-03 -2.497E-04 2.396E-03 -3.522E-03
        7.850E-03	6.438E-03	-2.831E-02	-2.916E-02	-7.322E-03	-6.840E-04	2.100E-03 -7.106E-04
        4.091E-03	-5.720E-03	3.162E-02	-5.890E-03	6.914E-04 -1.414E-04 -3.119E-04 3.825E-04
        -7.732E-03	-2.385E-02	4.630E-02	4.123E-03	2.540E-03 2.932E-04	-5.552E-05 6.271E-05
    -4.265E-03	1.108E-02	-4.877E-02	1.041E-03  -7.650E-04 -7.284E-05 7.969E-06 -2.631E-05
    ])


    K=Array([	2.153E+00	4.018E+01	2.105E+00	3.382E+00	-2.211E+02	9.331E-01 -7.617E+00	-2.707E+00][:])
    conicConstants = K .- 1.0
    names =["image", "stop", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "image"]
    apers = [Inf, 0.5 * 5.57/2.8, 1.2, 1.2, 1.2, 1.2, 2.1, 2.1, 2.5, 2.5, 2.7, 2.7, 3.]

    
    transform = translation(0.0, 0.0, 10.0)
    origins = Origins.Hexapolar(8, 10.0, 10.0)
    directions = Directions.Constant(0.0, 0.0, -1.0)
    raygenerator = Sources.Source(; transform, origins, directions)
    ix,jy = size(patent2)
    aspherics = [[string(j+2)=>patentTable2[i, j] for j in 1:jy] for i in 1:ix]
    materials = map(refIndex -> modelglass(refIndex, 1000.0, 0.0), patentTable1[:, 4])
    
    sys = AxisymmetricOpticalSystem{Float64}(DataFrame(
        SurfaceType = ["Object", "Stop", "Aspheric", "Aspheric","Aspheric","Aspheric","Aspheric","Aspheric","Aspheric","Aspheric","Standard", "Standard","Image"],
        Radius = patentTable1[:,2],
        Conic = [missing, missing, conicConstants..., missing, missing, missing],
        Parameters = [missing, missing, aspherics..., missing, missing,missing],
        Thickness = patentTable1[:,3],
        Material = materials,
        SemiDiameter = apers
        ))

    visplane = Plane(SVector(1.0, 0., 0.), SVector(0.0, 0.0,0.0))
    Vis.draw(sys ;resolution=(1600, 1000))  

    #Vis.drawtracerays(sys; raygenerator, trackallrays=true, test=true, numdivisions=50, resolution=(1600, 1000))
    #Vis.make2dy()
    #Vis.save(filename)

    return sys
end

draw_cellPhoneCameraLens()

##

brain = load(assetpath("brain.stl"))

b = Meshes.Mesh(brain)
Meshes.boundingbox(brain)

m = mesh(
    brain,
    color = [tri[1][2] for tri in brain for i in 1:3],
    colormap = Reverse(:magma)
)


r = Rect(0.0, 0.0, 0.0, 60.0, 60., 60.)

c = Cuboid(100.0, 10.0, 100.0)

r ∩ c

m1 = mesh(
    brain ∩ c,
    color = [tri[1][2] for tri in brain for i in 1:3],
    colormap = Reverse(:magma)
)

display(m)

