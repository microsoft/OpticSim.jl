@testset "Lenses" begin
    test_n = 5000

    """Creates a 3D vector uniformly distributed on the sphere by rejection sampling, i.e., discarding all points with norm > 1.0"""
    function randunit()
        let v = rand(3)
            while (norm(v) > 1.0)
                v = rand(3)
            end
            return normalize(v)
        end
    end

    @testset "Refraction" begin
        Random.seed!(SEED)
        ηᵢ = 1.4
        ηₜ = 1.2
        for i in 1:test_n
            nₛ = randunit()
            rᵢ = randunit()
            rₜ = refractedray(ηᵢ, ηₜ, nₛ, rᵢ)
            if !(rₜ === nothing)
                sinθₜ = norm(cross(-nₛ, rₜ))
                sinθᵢ = norm(cross(nₛ, rᵢ))
                a = isapprox(ηₜ * sinθₜ, ηᵢ * sinθᵢ, rtol = RTOLERANCE, atol = ATOLERANCE) #verify snell's law for the the transmitted and incident ray
                b = isapprox(1.0, norm(rₜ), rtol = RTOLERANCE, atol = ATOLERANCE) #verify transmitted ray is unit
                perp = normalize(cross(nₛ, rᵢ))
                c = isapprox(0.0, rₜ ⋅ perp, rtol = RTOLERANCE, atol = ATOLERANCE) #verify incident and transmitted ray are in the same plane
                @test a && b && c
            end
        end
    end # testset refraction

    # snell
    @testset "Snell" begin
        Random.seed!(SEED)
        for i in 1:test_n
            r = randunit()
            nₛ = randunit()
            n1 = 1.0
            n2 = 1.4
            sθ1, sθ2 = snell(nₛ, r, n1, n2)
            sθ3, sθ4 = snell(nₛ, r, n2, n1)
            @test isapprox(sθ1 * n1, sθ2 * n2, rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(sθ3 * n2, sθ4 * n1, rtol = RTOLERANCE, atol = ATOLERANCE)
        end
    end # testset snell

    @testset "Reflection" begin
        Random.seed!(SEED)
        # reflection
        for i in 1:test_n
            r = randunit()
            nₛ = randunit()
            if abs(r ⋅ nₛ) < 0.9999 # if r and n are exactly aligned then their sum will be zero, not a vector pointing in the direction of the normal
                reflected = reflectedray(nₛ, r)
                #ensure reflected and refracted rays are in the same plane
                a = isapprox(norm(reflected ⋅ cross(r, nₛ)), 0.0, rtol = RTOLERANCE, atol = ATOLERANCE)
                b = isapprox(norm(reflected), 1.0, rtol = RTOLERANCE, atol = ATOLERANCE)
                c = isapprox(0.0, reflected ⋅ nₛ + r ⋅ nₛ, rtol = RTOLERANCE, atol = ATOLERANCE)
                #ensure sum of reflected and origin ray are in the direction of the normal
                d = isapprox(0.0, norm(nₛ - sign(reflected ⋅ nₛ) * (normalize(reflected - r))), rtol = RTOLERANCE, atol = ATOLERANCE)
                @test a & b & c & d
            end
        end
    end # testset reflection

    @testset "Paraxial" begin
        # check that normally incident rays are focussed across the lens
        l = ParaxialLensEllipse(100.0, 10.0, 10.0, [1.0, 1.0, 0.0], [3.0, 3.0, 3.0])
        r = OpticalRay([2.0, 2.0, 3.0], [1.0, 1.0, 0.0], 1.0, 0.55)
        intsct = halfspaceintersection(surfaceintersection(l, r))
        ref, _, _ = OpticSim.processintersection(OpticSim.interface(intsct), OpticSim.point(intsct), OpticSim.normal(intsct), r, 20.0, 1.0, true, true)
        @test isapprox(ref, normalize([1.0, 1.0, 0.0]), rtol = RTOLERANCE, atol = ATOLERANCE)
        l = ParaxialLensRect(100.0, 10.0, 10.0, [1.0, 1.0, 0.0], [3.0, 3.0, 3.0])
        r = OpticalRay([2.3, 2.4, 3.1], [1.0, 1.0, 0.0], 1.0, 0.55)
        intsct = halfspaceintersection(surfaceintersection(l, r))
        fp = [3.0, 3.0, 3.0] + 100 * normalize([1.0, 1.0, 0.0])
        ref, _, _ = OpticSim.processintersection(OpticSim.interface(intsct), OpticSim.point(intsct), OpticSim.normal(intsct), r, 20.0, 1.0, true, true)
        @test isapprox(ref, normalize(fp - OpticSim.point(intsct)), rtol = RTOLERANCE, atol = ATOLERANCE)
    end # testset Paraxial
end # testset lenses
