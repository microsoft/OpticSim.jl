@testset "ParaxialAnalysis" begin
    @testset "Lens" begin
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
        @test isapprox(virtptfromslope, point(virtualpoint(lens,displaypoint)))
    end

    @testset "Projection" begin
        focallength = 10.0
        lens = ParaxialLensRect(focallength,100.0,100.0,[0.0,0.0,1.0],[0.0,0.0,0.0])
        display = Display(1000,1000,1.0μm,1.0μm,translation(0.0,0.0,-focallength))
        lenslet = LensletAssembly(lens,identitytransform(),display)
        displaypoint = SVector(0.0,0.0,-8.0)
        pupilpoints = SMatrix{3,2}(10.0,10.0,10.0,-10.0,-10.0,20.0)
        project(lenslet,displaypoint,pupilpoints)
    end

    @testset "BeamEnergy" begin
        focallength = 10.0
        lens = ParaxialLensRect(focallength,1.0,1.0,[0.0,0.0,1.0],[0.0,0.0,0.0])
        display = Display(1000,1000,1.0μm,1.0μm,translation(0.0,0.0,-focallength))
        lenslet = LensletAssembly(lens,identitytransform(),display)
        displaypoint = SVector(0.0,0.0,-8.0)
        pupil = Rectangle(1.0,1.0,SVector(0.0,0.0,-1.0),SVector(2.0,2.0,40.0))
    
        beamenergy(lenslet,displaypoint,Geometry.vertices3d(pupil))
    end

    @testset "SphericalPolygon" begin
        """creates a circular polygon that subtends a half angle of θ. If you double θ the spherical area should double"""
        function sphericalcircle(θ, nsides = 10)
            temp = MMatrix{3,nsides,Float64}(undef)
            for i in 0:1:(nsides-1)
                ϕ = i*2π/nsides
                temp[1,i+1] = sin(θ)*cos(ϕ)
                temp[2,i+1] = cos(θ)
                temp[3,i+1] = sin(θ)*sin(ϕ)
            end
            return SphericalPolygon(SMatrix{3,nsides,Float64}(temp),SVector(0.0,0.0,0.0),1.0)
        end

        # function circlepoly(nsides; offset = [0.0,0.0,1.0])
        #     step = 2π/nsides
        #     result = MMatrix{3,nsides,Float64}(undef)
        #     for i in 0:nsides-1
        #         y,x = sincos(step*i)
        #         result[:,i+1] = [x,y,0.0] .+ offset
        #     end
        #     return SMatrix{3,nsides,Float64}(result)
        # end
        
        oneeigthsphere() = SphericalTriangle(SMatrix{3,3,Float64}(
            0.0,1.0,0.0,
            1.0,0.0,0.0,
            0.0,0.0,1.0),
            SVector(0.0,0.0,0.0),
            1.0)
        
        onesixteenthphere() = SphericalTriangle(SMatrix{3,3,Float64}(
            0.0,1.0,0.0,
            1.0,1.0,0.0,
            0.0,0.0,1.0),
            SVector(0.0,0.0,0.0),
            1.0)
         
        threesidedpoly() = SphericalPolygon(SMatrix{3,3,Float64}(
            0.0,1.0,0.0,
            1.0,0.0,0.0,
            0.0,0.0,1.0),
            SVector(0.0,0.0,0.0),
            1.0)
         
        
        foursidedpoly() = SphericalPolygon(SMatrix{3,4,Float64}(
            0.0,1.0,0.0,
            1.0,1.0,-.9,
            1.0,0.0,0.0,
            0.0,0.0,1.0),
            SVector(0.0,0.0,0.0),
            1.0)
    
        @test isapprox(area(oneeigthsphere()), area(threesidedpoly()))
        @test isapprox(area(onesixteenthphere()), area(oneeigthsphere())/2.0)
        @test isapprox(area(circlepoly(π/8.0,1000)),area(circlepoly(π/16.0,1000))/4.0)
        
        
        end
    end
end

