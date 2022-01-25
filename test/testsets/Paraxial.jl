# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.


@testset "ParaxialLens" begin
    @testset "Virtual point" begin
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

    function rayintersection(lens,incomingray)
        intsct = OpticSim.surfaceintersection(lens,incomingray)
        refrac,_,_ = OpticSim.processintersection(OpticSim.interface(lens),OpticSim.point(intsct),OpticSim.normal(lens),OpticalRay(incomingray,1.0,.55),OpticSim.TEMP_REF,OpticSim.PRESSURE_REF,false)
        return refrac
    end

    @testset "Refracted rays" begin
         #test all 4 combinations: n⋅r > 0, n⋅r < 0, focal length > 0, focal length < 0
        #n⋅r > 0 focal length > 0
        lens = ParaxialLensRect(1.0,100.0,100.0,[0.0,0.0,1.0],[0.0,0.0,0.0])
        r = Ray([1.0,0.0,-1.0],[0.0,0.0,1.0])
        refrac = rayintersection(lens,r)
        @test isapprox([-sqrt(2)/2,0.0,sqrt(2)/2],refrac)
        
        #n⋅r < 0 focal length > 0
        r = Ray([1.0,0.0,1.0],[0.0,0.0,-1.0])
        refrac = rayintersection(lens,r)
        @test isapprox([-sqrt(2)/2,0.0,-sqrt(2)/2],refrac)

        #n⋅r > 0 focal length < 0
        lens = ParaxialLensRect(-1.0,100.0,100.0,[0.0,0.0,1.0],[0.0,0.0,0.0])
        r = Ray([1.0,0.0,-1.0],[0.0,0.0,1.0])
        refrac = rayintersection(lens,r)
        @test isapprox([sqrt(2)/2,0.0,sqrt(2)/2],refrac)
        
        #n⋅r < 0 focal length < 0
        r = Ray([1.0,0.0,1.0],[0.0,0.0,-1.0])
        refrac = rayintersection(lens,r)
        @test isapprox([sqrt(2)/2,0.0,-sqrt(2)/2],refrac)
    end

    @testset "Reversibility of rays for paraxial lenses" begin
        lens = ParaxialLensRect(10.0,100.0,100.0,[0.0,0.0,1.0],[0.0,0.0,0.0])
        displaypoint = [0.0,0.0,-8.0]
        direction = [0.0,0.0,1.0]
        r1 = Ray(displaypoint,direction)
        intsct = OpticSim.surfaceintersection(lens,r1)
        intsctpt = point(intsct)
        # compute the refracted ray and project this backwards to its intersection with the optical axis. This should be the same position as returned by virtualpoint()
        refrac1,_,_ = OpticSim.processintersection(OpticSim.interface(lens),OpticSim.point(intsct),OpticSim.normal(lens),OpticalRay(r1,1.0,.55),OpticSim.TEMP_REF,OpticSim.PRESSURE_REF,false)
        r2 = Ray(-displaypoint,-direction)

        intsct = OpticSim.surfaceintersection(lens,r2)
        intsctpt = point(intsct)
        # compute the refracted ray and project this backwards to its intersection with the optical axis. This should be the same position as returned by virtualpoint()
        refrac2,_,_ = OpticSim.processintersection(OpticSim.interface(lens),OpticSim.point(intsct),OpticSim.normal(lens),OpticalRay(r2,1.0,.55),OpticSim.TEMP_REF,OpticSim.PRESSURE_REF,false)

        @test isapprox(-refrac2,refrac1)
    end
end