@testset "SurfaceDefs" begin
    @testset "Knots" begin
        # find span
        knots = KnotVector{Int64}([0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5])
        function apply(knots::KnotVector, curveorder, u, correctindex)
            knotnum = findspan(knots, curveorder, u)
            return knotnum == correctindex
        end
        @test apply(knots, 2, 2.5, 5) && apply(knots, 2, 0, 3) && apply(knots, 2, 4, 6) && apply(knots, 2, 5, 8)

        # insert knots
        knots = [0, 0, 0, 0, 1, 2, 3, 3, 4, 4, 5, 5, 6, 7, 7, 7, 7]
        insertions = knotstoinsert(knots, 3)
        @test insertions == [(5, 2), (6, 2), (7, 1), (9, 1), (11, 1), (13, 2)]
    end

    @testset "BSpline" begin
        # bspline conversion
        correctverts = [[0.0, 0.0, 0.0], [0.0, 0.49625, 0.0], [0.6625, 0.49625, 0.65625], [0.0, 0.0, 0.0], [0.6625, 0.49625, 0.65625], [0.6625, 0.0, 0.0], [0.6625, 0.0, 0.0], [0.6625, 0.49625, 0.65625], [1.33, 0.49625, 0.0], [0.6625, 0.0, 0.0], [1.33, 0.49625, 0.0], [1.33, 0.0, 0.0], [0.0, 0.49625, 0.0], [0.0, 1.0, 0.0], [0.6625, 1.0, 0.0], [0.0, 0.49625, 0.0], [0.6625, 1.0, 0.0], [0.6625, 0.49625, 0.65625], [0.6625, 0.49625, 0.65625], [0.6625, 1.0, 0.0], [1.33, 1.0, 0.0], [0.6625, 0.49625, 0.65625], [1.33, 1.0, 0.0], [1.33, 0.49625, 0.0]]
        correcttris = [
            0x00000001 0x00000002 0x00000003
            0x00000004 0x00000005 0x00000006
            0x00000007 0x00000008 0x00000009
            0x0000000a 0x0000000b 0x0000000c
            0x0000000d 0x0000000e 0x0000000f
            0x00000010 0x00000011 0x00000012
            0x00000013 0x00000014 0x00000015
            0x00000016 0x00000017 0x00000018
        ]
        surf = TestData.bsplinesurface()
        points, triangles = makiemesh(makemesh(surf, 2))
        segments = tobeziersegments(surf) # TODO just testing that this evaluates for now
        @test triangles == correcttris
        @test all(isapprox.(points, correctverts, rtol = RTOLERANCE, atol = ATOLERANCE))
    end # testset BSpline

    @testset "PowerBasis" begin
        # polynomial is 3x^2 + 2x + 1
        coeff = [1.0 2.0 3.0]
        curve = PowerBasisCurve{OpticSim.Euclidean,Float64,1,2}(coeff)
        @test !all(isapprox.(coeff, coefficients(curve, 1), rtol = RTOLERANCE, atol = ATOLERANCE)) # TODO not sure if this is correct/what this is testing?
        correct = true
        for x in -10.0:0.1:10
            exact = 3 * x^2 + 2 * x + 1
            calculated = point(curve, x)[1]
            @test isapprox(exact, calculated, rtol = RTOLERANCE, atol = ATOLERANCE)
        end
    end # testset PowerBasis

    @testset "BezierSurface" begin
        surf = TestData.beziersurface()
        fdm = central_fdm(10, 1)
        for u in 0:0.1:1, v in 0:0.1:1
            (du, dv) = partials(surf, u, v)
            fu(u) = point(surf, u, v)
            fv(v) = point(surf, u, v)
            # compute accurate finite difference approximation to derivative
            fdu, fdv = (fdm(fu, u), fdm(fv, v))
            @test isapprox(du, fdu, rtol = 1e-12, atol = 2 * eps(Float64))
            @test isapprox(dv, fdv, rtol = 1e-12, atol = 2 * eps(Float64))
        end
    end # testset BezierSurface

    @testset "ZernikeSurface" begin
        surf = TestData.zernikesurface1a() # with normradius
        fdm = central_fdm(10, 1)
        for ρ in 0.05:0.1:0.95, ϕ in 0:(π / 10):(2π)
            (dρ, dϕ) = partials(surf, ρ, ϕ)
            fρ(ρ) = point(surf, ρ, ϕ)
            fϕ(ϕ) = point(surf, ρ, ϕ)
            # compute accurate finite difference approximation to derivative
            fdρ, fdϕ = (fdm(fρ, ρ), fdm(fϕ, ϕ))
            @test isapprox(dρ, fdρ, rtol = 1e-12, atol = 2 * eps(Float64))
            @test isapprox(dϕ, fdϕ, rtol = 1e-12, atol = 2 * eps(Float64))
        end

        @test !any(isnan.(normal(TestData.conicsurface(), 0.0, 0.0)))
    end # testset ZernikeSurface

    @testset "QTypeSurface" begin
        # test predefined values against papers
        # Forbes, G. W. "Robust, efficient computational methods for axially symmetric optical aspheres." OpticSim express 18.19 (2010): 19700-19712.
        # Forbes, G. W. "Characterizing the shape of freeform optics." OpticSim express 20.3 (2012): 2483-2499.

        f0_true = [2, sqrt(19 / 4), 4 * sqrt(10 / 19), 1 / 2 * sqrt(509 / 10), 6 * sqrt(259 / 509), 1 / 2 * sqrt(25607 / 259)]
        for i in 1:length(f0_true)
            @test isapprox(OpticSim.QType.f0(i - 1), f0_true[i], rtol = 2 * eps(Float64))
        end
        g0_true = [-1 / 2, -5 / (2 * sqrt(19)), -17 / (2 * sqrt(190)), -91 / (2 * sqrt(5090)), -473 / (2 * sqrt(131831))]
        for i in 1:length(g0_true)
            @test isapprox(OpticSim.QType.g0(i - 1), g0_true[i], rtol = 2 * eps(Float64))
        end
        h0_true = [-1 / 2, -6 / sqrt(19), -3 / 2 * sqrt(19 / 10), -20 * sqrt(10 / 509)]
        for i in 1:length(h0_true)
            @test isapprox(OpticSim.QType.h0(i - 1), h0_true[i], rtol = 2 * eps(Float64))
        end

        F_true = [1/4 1/2 27/32 5/4; 15/32 7/8 35/16 511/128; 17/72 35/36 35/16 23/6; 29/40 67/40 243/80 12287/2560]
        N, M = size(F_true)
        for m in 1:M
            for n in 0:(N - 1)
                @test isapprox(OpticSim.QType.F(m, n), F_true[n + 1, m], rtol = 2 * eps(Float64))
            end
        end
        G_true = [1/4 3/8 15/32 35/64; -1/24 -5/48 -7/32 -21/64; -7/40 -7/16 -117/160 -33/32]
        N, M = size(G_true)
        for m in 1:M
            for n in 0:(N - 1)
                @test isapprox(OpticSim.QType.G(m, n), G_true[n + 1, m], rtol = 2 * eps(Float64))
            end
        end
        f_true = [1/2 1/sqrt(2) 3 / 4*sqrt(3 / 2) sqrt(5)/2; 1 / 4*sqrt(7 / 2) 1 / 4*sqrt(19 / 2) 1 / 4*sqrt(185 / 6) 3 / 32*sqrt(427); 1 / 6*sqrt(115 / 14) 1 / 2*sqrt(145 / 38) 1 / 4*sqrt(12803 / 370) 1 / 2*sqrt(2785 / 183); 1 / 5*sqrt(3397 / 230) 1 / 4*sqrt(6841 / 290) 9 / 4*sqrt(14113 / 25606) 1 / 16*sqrt(1289057 / 1114)]
        N, M = size(f_true)
        for m in 1:M
            for n in 0:(N - 1)
                @test isapprox(OpticSim.QType.f(m, n), f_true[n + 1, m], rtol = 2 * eps(Float64))
            end
        end
        g_true = [1/2 3/(4 * sqrt(2)) 5/(4 * sqrt(6)) 7 * sqrt(5)/32; -1/(3 * sqrt(14)) -5/(6 * sqrt(38)) -7 / 4*sqrt(3 / 370) -1 / 2*sqrt(7 / 61); -21 / 10*sqrt(7 / 230) -7 / 4*sqrt(19 / 290) -117 / 4*sqrt(37 / 128030) -33 / 16*sqrt(183 / 2785)]
        N, M = size(g_true)
        for m in 1:M
            for n in 0:(N - 1)
                @test isapprox(OpticSim.QType.g(m, n), g_true[n + 1, m], rtol = 2 * eps(Float64))
            end
        end
        A_true = [2 3 5 7 9 11; -4/3 2/3 2 76/27 85/24 106/25; 9/5 26/15 2 117/50 203/75 108/35; 55/28 66/35 2 536/245 135/56 130/49; 161/81 122/63 2 515/243 143/63 161/66]
        N, M = size(A_true)
        for m in 1:M
            for n in 0:(N - 1)
                @test isapprox(OpticSim.QType.A(m, n), A_true[n + 1, m], rtol = 2 * eps(Float64))
            end
        end
        B_true = [-1 -2 -4 -6 -8 -10; -8/3 -4 -4 -40/9 -5 -28/5; -24/5 -4 -4 -21/5 -112/25 -24/5; -30/7 -4 -4 -144/35 -30/7 -220/49; -112/27 -4 -4 -110/27 -88/21 -13/3]
        N, M = size(B_true)
        for m in 1:M
            for n in 0:(N - 1)
                @test isapprox(OpticSim.QType.B(m, n), B_true[n + 1, m], rtol = 2 * eps(Float64))
            end
        end
        C_true = [NaN NaN NaN NaN NaN NaN; -11/3 -3 -5/3 -35/27 -9/8 -77/75; 0 5/9 7/15 21/50 88/225 13/35; 27/28 21/25 27/35 891/1225 39/56 33/49; 80/81 45/49 55/63 1430/1701 40/49 1105/1386]
        N, M = size(C_true)
        for m in 1:M
            for n in 0:(N - 1)
                @test isapprox(OpticSim.QType.C(m, n), C_true[n + 1, m], rtol = 2 * eps(Float64)) || (isnan(OpticSim.QType.C(m, n)) && isnan(C_true[n + 1, m]))
            end
        end

        # test deriv with finite differences
        surf = TestData.qtypesurface1()
        fdm = central_fdm(10, 1)
        for ρ in 0.05:0.1:0.95, ϕ in 0:(π / 10):(2π)
            (dρ, dϕ) = partials(surf, ρ, ϕ)
            fρ(ρ) = point(surf, ρ, ϕ)
            fϕ(ϕ) = point(surf, ρ, ϕ)
            # compute accurate finite difference approximation to derivative
            fdρ, fdϕ = (fdm(fρ, ρ), fdm(fϕ, ϕ))
            @test isapprox(dρ, fdρ, rtol = 1e-12, atol = 2 * eps(Float64))
            @test isapprox(dϕ, fdϕ, rtol = 1e-12, atol = 2 * eps(Float64))
        end
    end # testset QTypeSurface

    @testset "GridSag" begin
        surf = TestData.gridsagsurfacelinear()
        fdm = central_fdm(10, 1)
        # doesn't enforce C1 across patch boundaries meaning that finite differences won't match at all
        # just test within a patch to make sure partials() is working
        for ρ in 0.02:0.01:0.18, ϕ in 0:(π / 30):(2π)
            (dρ, dϕ) = partials(surf, ρ, ϕ)
            fρ(ρ) = point(surf, ρ, ϕ)
            fϕ(ϕ) = point(surf, ρ, ϕ)
            # compute accurate finite difference approximation to derivative
            fdρ, fdϕ = (fdm(fρ, ρ), fdm(fϕ, ϕ))
            @test isapprox(dρ, fdρ, rtol = 1e-12, atol = 2 * eps(Float64))
            @test isapprox(dϕ, fdϕ, rtol = 1e-12, atol = 2 * eps(Float64))
        end

        surf = TestData.gridsagsurfacebicubic()
        fdm = central_fdm(10, 1)
        # doesn't enforce C2 across patch boundaries meaning that finite differences won't match exactly
        # just test within a patch to make sure partials() is working
        for ρ in 0.02:0.01:0.18, ϕ in 0:(π / 30):(2π)
            (dρ, dϕ) = partials(surf, ρ, ϕ)
            fρ(ρ) = point(surf, ρ, ϕ)
            fϕ(ϕ) = point(surf, ρ, ϕ)
            # compute accurate finite difference approximation to derivative
            fdρ, fdϕ = (fdm(fρ, ρ), fdm(fϕ, ϕ))
            @test isapprox(dρ, fdρ, rtol = 1e-12, atol = 2 * eps(Float64))
            @test isapprox(dϕ, fdϕ, rtol = 1e-12, atol = 2 * eps(Float64))
        end

        surf = TestData.gridsagsurfacebicubiccheby()
        fdm = central_fdm(10, 1)
        # not valid at the very boundary of the surface as stuff gets weird outside |u| < 1 and |v| < 1
        for u in -0.3:0.01:0.3, v in -0.15:0.01:0.15
            (du, dv) = partials(surf, u, v)
            fu(u) = point(surf, u, v)
            fv(v) = point(surf, u, v)
            # compute accurate finite difference approximation to derivative
            fdu, fdv = (fdm(fu, u), fdm(fv, v))
            @test isapprox(du, fdu, rtol = 1e-12, atol = 2 * eps(Float64))
            @test isapprox(dv, fdv, rtol = 1e-12, atol = 2 * eps(Float64))
        end
    end # testset GridSag

    @testset "ChebyshevSurface" begin
        surf = TestData.chebyshevsurface1() # with normradius
        fdm = central_fdm(10, 1)
        # not valid at the very boundary of the surface as stuff gets weird outside |u| < 1 and |v| < 1
        for u in -0.99:0.1:0.99, v in -0.99:0.1:0.99
            (du, dv) = partials(surf, u, v)
            fu(u) = point(surf, u, v)
            fv(v) = point(surf, u, v)
            # compute accurate finite difference approximation to derivative
            fdu, fdv = (fdm(fu, u), fdm(fv, v))
            @test isapprox(du, fdu, rtol = 1e-11, atol = 2 * eps(Float64)) #changed rtol from 1e-12 to 1e-11. FiniteDifferences approximation to the derivative had larger than expected error.
            @test isapprox(dv, fdv, rtol = 1e-11, atol = 2 * eps(Float64)) #changed rtol from 1e-12 to 1e-11. FiniteDifferences approximation to the derivative had larger than expected error.
        end
    end # testset ChebyshevSurface

end # testset SurfaceDefs
