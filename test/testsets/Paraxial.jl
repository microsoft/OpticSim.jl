@testset "ParaxialLens" begin
    @testset "Virtual point" begin
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
            @test isapprox(virtptfromslope,point(virtualpoint(lens,displaypoint)))
        end
    end
end