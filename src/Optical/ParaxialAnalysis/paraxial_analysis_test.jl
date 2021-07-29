function testvirtualpoint()
    lens = ParaxialLensRect(10.0,100.0,100.0,[0.0,0.0,1.0],[0.0,0.0,0.0])
    displaypoint = [0.0,0.0,-8.0]
    r1 = Ray(displaypoint,[1.0,0.0,1.0])
    intsct = OpticSim.surfaceintersection(lens,r1)
    intsctpt = point(intsct)
    # compute the refracted ray and project this backwards to its intersection with the optical axis. This should be the same position as returned by virtualpoint()
    refrac,_,_ = OpticSim.processintersection(OpticSim.interface(lens),OpticSim.point(intsct),OpticSim.normal(lens),OpticalRay(r1,1.0,.55),OpticSim.TEMP_REF,OpticSim.PRESSURE_REF,false)
    # compute intersection of ray with optical axis. this should match the position of the virtual point
    slope = refrac[1]/refrac[3]
    virtptfromslope = [0.0,0.0,-intsctpt[1]/slope]
    # insert test virtptfromslope == point(virtualpoint(lens,displaypoint))
    println(virtptfromslope)
    println(point(virtualpoint(lens,displaypoint)))
end
export testvirtualpoint

function testproject()
    focallength = 10.0
    lens = ParaxialLensRect(focallength,100.0,100.0,[0.0,0.0,1.0],[0.0,0.0,0.0])
    display = Display(1000,1000,1.0μm,1.0μm,translation(0.0,0.0,-focallength))
    lenslet = LensletAssembly(lens,identitytransform(),display)
    displaypoint = SVector(0.0,0.0,-8.0)
    pupilpoints = SMatrix{3,2}(10.0,10.0,10.0,-10.0,-10.0,20.0)
    project(lenslet,displaypoint,pupilpoints)
end
export testproject

function testintersection()
    lens =  ParaxialLensRect(10.0,1.0,1.0,[0.0,0.0,1.0],[0.0,0.0,0.0])
    lensverts = vertices(lens)
    # projectedverts = SVector{4,SVector{3,Float64}}(SVector(0.0, 0.0, 0.0),SVector(0.0,2.0, 0.0),SVector(2.0,2.0, 0.0),SVector(2.0,0.0, 0.0))
    projectedverts = SMatrix{2,4}(0.0,0.0,0.0,2.0,2.0,2.0,2.0,0.0)
    
    LazySets.VPolygon(projectedverts) ∩ LazySets.VPolygon(lensverts)
end
export testintersection

