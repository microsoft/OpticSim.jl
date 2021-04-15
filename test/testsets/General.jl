@testset "General" begin
    @testset "QuadraticRoots" begin
        Random.seed!(SEED)
        similarroots(r1, r2, (x1, x2)) = isapprox([r1, r2], [x1, x2], rtol = 1e-9) || isapprox([r2, r1], [x1, x2], rtol = 1e-9)

        @testcase "Fuzz Test" begin
            for i in 1:10000
                r1, r2, scale = rand(3) .- 0.5
                a, b, c = scale .* (1, r1 + r2, r1 * r2)
                @test similarroots(-r1, -r2, quadraticroots(a, b, c))
            end
        end

        @testcase "Unit Test" begin
            @test similarroots(1, -1, quadraticroots(1, 0, -1))
        end

        @testcase "Repeated Roots" begin
            # (x+1)(x+1)= x^2 + 2x + 1, double root at one
            @test similarroots(-1, -1, quadraticroots(1, 2, 1))
        end
    end # testset QuadraticRoots

    @testset "Transform" begin
        Random.seed!(SEED)

        @testcase begin
            @test isapprox(rotmatd(180, 0, 0), [1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 -1.0], rtol = RTOLERANCE, atol = ATOLERANCE)
            @test isapprox(rotmatd(0.0, 180.0, 0.0), [-1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 -1.0], rtol = RTOLERANCE, atol = ATOLERANCE)
            @test isapprox(rotmatd(0, 0, 180), [-1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 1.0], rtol = RTOLERANCE, atol = ATOLERANCE)
            @test isapprox(rotmatd(0, 90, 0), [0.0 0.0 1.0; 0.0 1.0 0.0; -1.0 0.0 0.0], rtol = RTOLERANCE, atol = ATOLERANCE)
            @test isapprox(rotmatd(45, -45, 45), [0.5 -0.8535533905932737 0.1464466094067261; 0.5 0.14644660940672644 -0.8535533905932737; 0.7071067811865475 0.5 0.5], rtol = RTOLERANCE, atol = ATOLERANCE)
        end

        @testcase begin
            x, y, z = rand(3)
            @test isapprox(rotmatd(x * 180 / π, y * 180 / π, z * 180 / π), rotmat(x, y, z), rtol = RTOLERANCE, atol = ATOLERANCE)
        end

        @testcase begin
            @test isapprox(Transform(rotmatd(0, 90, 0), SVector(0.0, 0.0, 1.0)) * SVector(1.0, 0.0, 0.0), [0.0, 0.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE)
        end

        @testcase begin
            ta = Transform(rotmatd(0, 90, 0), SVector(0.0, 0.0, 1.0))
            tb = Transform(rotmatd(90, 0, 0), SVector(1.0, 0.0, 0.0))
            @test isapprox(collect(ta * tb), collect(Transform(rotmatd(90, 90, 0), SVector(0.0, 0.0, 0.0))), rtol = RTOLERANCE, atol = ATOLERANCE)
            @test isapprox(collect(ta * inv(ta)), collect(identitytransform()), rtol = RTOLERANCE, atol = ATOLERANCE)
            @test isapprox(collect(tb * inv(tb)), collect(identitytransform()), rtol = RTOLERANCE, atol = ATOLERANCE)
            @test isapprox(collect(inv(ta)), collect(Transform(rotmatd(0, -90, 0), SVector(1.0, 0.0, 0.0))), rtol = RTOLERANCE, atol = ATOLERANCE)
        end
    end # testset Transform

    @testset "Interval" begin
        intersectionat = TestData.intersectionat

        pt1 = Intersection(0.3, [4.0, 5.0, 6.0], normalize(rand(3)), 0.4, 0.5, NullInterface())
        pt2 = Intersection(0.5, [1.0, 2.0, 3.0], normalize(rand(3)), 0.2, 0.3, NullInterface())

        @testcase begin
            @test pt1 < pt2 && pt1 <= pt2 && pt2 > pt1 && pt2 >= pt1 && pt1 != pt2 && pt1 <= pt1
        end

        @testcase begin
            intvl2 = positivehalfspace(pt1)
            @test α(halfspaceintersection(intvl2)) == α(pt1)
        end

        @testset "interval intersection" begin
            @testcase "interval/interval" begin
                # one in another
                a = Interval(intersectionat(0.1), intersectionat(0.2))
                b = Interval(intersectionat(0.0), intersectionat(0.3))
                @test intervalintersection(a, b) == intervalintersection(b, a) == a
                # overlap
                a = Interval(intersectionat(0.1), intersectionat(0.3))
                b = Interval(intersectionat(0.2), intersectionat(0.4))
                @test intervalintersection(a, b) == intervalintersection(b, a) == Interval(intersectionat(0.2), intersectionat(0.3))
                # no overlap
                a = Interval(intersectionat(0.1), intersectionat(0.2))
                b = Interval(intersectionat(0.3), intersectionat(0.4))
                @test intervalintersection(a, b) == intervalintersection(b, a) isa EmptyInterval
                # start == end
                a = Interval(intersectionat(0.1), intersectionat(0.2))
                b = Interval(intersectionat(0.2), intersectionat(0.3))
                @test intervalintersection(a, b) == intervalintersection(b, a) isa EmptyInterval
            end

            @testcase "interval/du" begin
                # encompass one
                a = Interval(intersectionat(0.0), intersectionat(0.3))
                b = DisjointUnion(Interval(intersectionat(0.1), intersectionat(0.2)), Interval(intersectionat(0.5), intersectionat(0.6)))
                @test intervalintersection(a, b) == intervalintersection(b, a) == b[1]
                # encompass two
                a = Interval(intersectionat(0.0), intersectionat(0.8))
                b = DisjointUnion(Interval(intersectionat(0.1), intersectionat(0.2)), Interval(intersectionat(0.5), intersectionat(0.6)))
                @test intervalintersection(a, b) == intervalintersection(b, a) == b
                # overlap one
                a = Interval(intersectionat(0.0), intersectionat(0.2))
                b = DisjointUnion(Interval(intersectionat(0.1), intersectionat(0.3)), Interval(intersectionat(0.5), intersectionat(0.6)))
                @test intervalintersection(a, b) == intervalintersection(b, a) == Interval(intersectionat(0.1), intersectionat(0.2))
                # overlap two
                a = Interval(intersectionat(0.2), intersectionat(0.5))
                b = DisjointUnion(Interval(intersectionat(0.1), intersectionat(0.3)), Interval(intersectionat(0.4), intersectionat(0.6)))
                res = intervalintersection(a, b)
                @test res == intervalintersection(b, a) && res[1] == Interval(intersectionat(0.2), intersectionat(0.3)) && res[2] == Interval(intersectionat(0.4), intersectionat(0.5))
                # no overlap
                a = Interval(intersectionat(0.1), intersectionat(0.2))
                b = DisjointUnion(Interval(intersectionat(0.3), intersectionat(0.4)), Interval(intersectionat(0.5), intersectionat(0.6)))
                @test intervalintersection(a, b) == intervalintersection(b, a) isa EmptyInterval
            end

            @testcase "du/du" begin
                # no overlap
                a = DisjointUnion(Interval(intersectionat(0.0), intersectionat(0.1)), Interval(intersectionat(0.4), intersectionat(0.5)))
                b = DisjointUnion(Interval(intersectionat(0.2), intersectionat(0.3)), Interval(intersectionat(0.6), intersectionat(0.7)))
                @test intervalintersection(a, b) == intervalintersection(b, a) isa EmptyInterval
                # one in a overlaps one in b
                a = DisjointUnion(Interval(intersectionat(0.0), intersectionat(0.2)), Interval(intersectionat(0.4), intersectionat(0.5)))
                b = DisjointUnion(Interval(intersectionat(0.1), intersectionat(0.3)), Interval(intersectionat(0.6), intersectionat(0.7)))
                @test intervalintersection(a, b) == intervalintersection(b, a) == Interval(intersectionat(0.1), intersectionat(0.2))
                # one in a encompasses one in b
                a = DisjointUnion(Interval(intersectionat(0.0), intersectionat(0.3)), Interval(intersectionat(0.4), intersectionat(0.5)))
                b = DisjointUnion(Interval(intersectionat(0.1), intersectionat(0.2)), Interval(intersectionat(0.6), intersectionat(0.7)))
                @test intervalintersection(a, b) == intervalintersection(b, a) == b[1]
                # one in a overlaps two in b
                a = DisjointUnion(Interval(intersectionat(0.2), intersectionat(0.5)), Interval(intersectionat(0.7), intersectionat(0.8)))
                b = DisjointUnion(Interval(intersectionat(0.1), intersectionat(0.3)), Interval(intersectionat(0.4), intersectionat(0.6)))
                res = intervalintersection(a, b)
                @test res == intervalintersection(b, a) && res[1] == Interval(intersectionat(0.2), intersectionat(0.3)) && res[2] == Interval(intersectionat(0.4), intersectionat(0.5))
                # one in a encompasses two in b
                a = DisjointUnion(Interval(intersectionat(0.0), intersectionat(0.5)), Interval(intersectionat(0.7), intersectionat(0.8)))
                b = DisjointUnion(Interval(intersectionat(0.1), intersectionat(0.2)), Interval(intersectionat(0.3), intersectionat(0.4)))
                @test intervalintersection(a, b) == intervalintersection(b, a) == b
                # each in a overlap one in b
                a = DisjointUnion(Interval(intersectionat(0.0), intersectionat(0.2)), Interval(intersectionat(0.4), intersectionat(0.6)))
                b = DisjointUnion(Interval(intersectionat(0.1), intersectionat(0.3)), Interval(intersectionat(0.5), intersectionat(0.7)))
                res = intervalintersection(a, b)
                @test res == intervalintersection(b, a) && res[1] == Interval(intersectionat(0.1), intersectionat(0.2)) && res[2] == Interval(intersectionat(0.5), intersectionat(0.6))
                # each in a encompass one in b
                a = DisjointUnion(Interval(intersectionat(0.0), intersectionat(0.3)), Interval(intersectionat(0.4), intersectionat(0.7)))
                b = DisjointUnion(Interval(intersectionat(0.1), intersectionat(0.2)), Interval(intersectionat(0.5), intersectionat(0.6)))
                @test intervalintersection(a, b) == intervalintersection(b, a) == b
            end

            @testcase "empties" begin
                intvl = positivehalfspace(pt1)
                du = intervalcomplement(Interval(pt1, pt2))
                @test intervalintersection(EmptyInterval(), EmptyInterval()) isa EmptyInterval
                @test intervalintersection(EmptyInterval(), intvl) isa EmptyInterval
                @test intervalintersection(intvl, EmptyInterval()) isa EmptyInterval
                @test intervalintersection(EmptyInterval(), du) isa EmptyInterval
                @test intervalintersection(du, EmptyInterval()) isa EmptyInterval
            end
        end

        @testset "interval union" begin
            @testcase "interval/interval" begin
                # one in another
                a = Interval(intersectionat(0.1), intersectionat(0.2))
                b = Interval(intersectionat(0.0), intersectionat(0.3))
                @test intervalunion(a, b) == intervalunion(b, a) == b
                # overlap
                a = Interval(intersectionat(0.1), intersectionat(0.3))
                b = Interval(intersectionat(0.2), intersectionat(0.4))
                @test intervalunion(a, b) == intervalunion(b, a) == Interval(intersectionat(0.1), intersectionat(0.4))
                # no overlap
                a = Interval(intersectionat(0.1), intersectionat(0.2))
                b = Interval(intersectionat(0.3), intersectionat(0.4))
                res = intervalunion(a, b)
                @test res == intervalunion(b, a) && res[1] == a && res[2] == b
                # start == end
                a = Interval(intersectionat(0.1), intersectionat(0.2))
                b = Interval(intersectionat(0.2), intersectionat(0.3))
                @test intervalunion(a, b) == intervalunion(b, a) == Interval(intersectionat(0.1), intersectionat(0.3))
            end

            @testcase "interval/du" begin
                # encompass one
                a = Interval(intersectionat(0.0), intersectionat(0.3))
                b = DisjointUnion(Interval(intersectionat(0.1), intersectionat(0.2)), Interval(intersectionat(0.5), intersectionat(0.6)))
                res = intervalunion(a, b)
                @test res == intervalunion(b, a) && res[1] == a && res[2] == b[2]
                # encompass two
                a = Interval(intersectionat(0.0), intersectionat(0.8))
                b = DisjointUnion(Interval(intersectionat(0.1), intersectionat(0.2)), Interval(intersectionat(0.5), intersectionat(0.6)))
                @test intervalunion(a, b) == intervalunion(b, a) == a
                # overlap one
                a = Interval(intersectionat(0.0), intersectionat(0.2))
                b = DisjointUnion(Interval(intersectionat(0.1), intersectionat(0.3)), Interval(intersectionat(0.5), intersectionat(0.6)))
                res = intervalunion(a, b)
                @test res == intervalunion(b, a) && res[1] == Interval(intersectionat(0.0), intersectionat(0.3)) && res[2] == b[2]
                # overlap two
                a = Interval(intersectionat(0.2), intersectionat(0.5))
                b = DisjointUnion(Interval(intersectionat(0.1), intersectionat(0.3)), Interval(intersectionat(0.4), intersectionat(0.6)))
                @test intervalunion(a, b) == intervalunion(b, a) == Interval(intersectionat(0.1), intersectionat(0.6))
                # no overlap
                a = Interval(intersectionat(0.1), intersectionat(0.2))
                b = DisjointUnion(Interval(intersectionat(0.3), intersectionat(0.4)), Interval(intersectionat(0.5), intersectionat(0.6)))
                res = intervalunion(a, b)
                @test res == intervalunion(b, a) && res[1] == a && res[2] == b[1] && res[3] == b[2]
            end

            @testcase "du/du" begin
                # no overlap
                a = DisjointUnion(Interval(intersectionat(0.0), intersectionat(0.1)), Interval(intersectionat(0.4), intersectionat(0.5)))
                b = DisjointUnion(Interval(intersectionat(0.2), intersectionat(0.3)), Interval(intersectionat(0.6), intersectionat(0.7)))
                res = intervalunion(a, b)
                @test res == intervalunion(b, a) && res[1] == a[1] && res[2] == b[1] && res[3] == a[2] && res[4] == b[2]
                # one in a overlaps one in b
                a = DisjointUnion(Interval(intersectionat(0.0), intersectionat(0.2)), Interval(intersectionat(0.4), intersectionat(0.5)))
                b = DisjointUnion(Interval(intersectionat(0.1), intersectionat(0.3)), Interval(intersectionat(0.6), intersectionat(0.7)))
                res = intervalunion(a, b)
                @test res == intervalunion(b, a) && res[1] == Interval(intersectionat(0.0), intersectionat(0.3)) && res[2] == a[2] && res[3] == b[2]
                # one in a encompasses one in b
                a = DisjointUnion(Interval(intersectionat(0.0), intersectionat(0.3)), Interval(intersectionat(0.4), intersectionat(0.5)))
                b = DisjointUnion(Interval(intersectionat(0.1), intersectionat(0.2)), Interval(intersectionat(0.6), intersectionat(0.7)))
                res = intervalunion(a, b)
                @test res == intervalunion(b, a) && res[1] == a[1] && res[2] == a[2] && res[3] == b[2]
                # one in a overlaps two in b
                a = DisjointUnion(Interval(intersectionat(0.2), intersectionat(0.5)), Interval(intersectionat(0.7), intersectionat(0.8)))
                b = DisjointUnion(Interval(intersectionat(0.1), intersectionat(0.3)), Interval(intersectionat(0.4), intersectionat(0.6)))
                res = intervalunion(a, b)
                @test res == intervalunion(b, a) && res[1] == Interval(intersectionat(0.1), intersectionat(0.6)) && res[2] == a[2]
                # one in a encompasses two in b
                a = DisjointUnion(Interval(intersectionat(0.0), intersectionat(0.5)), Interval(intersectionat(0.7), intersectionat(0.8)))
                b = DisjointUnion(Interval(intersectionat(0.1), intersectionat(0.2)), Interval(intersectionat(0.3), intersectionat(0.4)))
                @test intervalunion(a, b) == intervalunion(b, a) == a
                # each in a overlap one in b
                a = DisjointUnion(Interval(intersectionat(0.0), intersectionat(0.2)), Interval(intersectionat(0.4), intersectionat(0.6)))
                b = DisjointUnion(Interval(intersectionat(0.1), intersectionat(0.3)), Interval(intersectionat(0.5), intersectionat(0.7)))
                res = intervalunion(a, b)
                @test res == intervalunion(b, a) && res[1] == Interval(intersectionat(0.0), intersectionat(0.3)) && res[2] == Interval(intersectionat(0.4), intersectionat(0.7))
                # each in a encompass one in b
                a = DisjointUnion(Interval(intersectionat(0.0), intersectionat(0.3)), Interval(intersectionat(0.4), intersectionat(0.7)))
                b = DisjointUnion(Interval(intersectionat(0.1), intersectionat(0.2)), Interval(intersectionat(0.5), intersectionat(0.6)))
                @test intervalunion(a, b) == intervalunion(b, a) == a
            end

            @testcase "empties" begin
                intvl = positivehalfspace(pt1)
                du = intervalcomplement(Interval(pt1, pt2))
                @test intervalunion(EmptyInterval(), EmptyInterval()) isa EmptyInterval
                @test intervalunion(EmptyInterval(), intvl) == intvl
                @test intervalunion(intvl, EmptyInterval()) == intvl
                @test intervalunion(EmptyInterval(), du) == du
                @test intervalunion(du, EmptyInterval()) == du
            end
        end

        @testset "interval complement" begin
            @testcase begin
                int = intervalcomplement(EmptyInterval())
                @test lower(int) isa RayOrigin && upper(int) isa Infinity
            end

            @testcase begin
                int = intervalcomplement(rayorigininterval(Infinity()))
                @test int isa EmptyInterval
            end

            @testcase begin
                int = intervalcomplement(rayorigininterval(pt1))
                @test lower(int) == pt1 && upper(int) isa Infinity
            end

            @testcase begin
                int = intervalcomplement(positivehalfspace(pt1))
                @test lower(int) isa RayOrigin && upper(int) == pt1
            end

            @testcase begin
                int = intervalcomplement(Interval(pt1, pt2))
                @test lower(int[1]) isa RayOrigin && upper(int[1]) == pt1 && lower(int[2]) == pt2 && upper(int[2]) isa Infinity
            end

            @testcase begin
                a = DisjointUnion(Interval(intersectionat(0.1), intersectionat(0.3)), Interval(intersectionat(0.4), intersectionat(0.7)))
                res = intervalcomplement(a)
                @test res[1] == rayorigininterval(intersectionat(0.1)) && res[2] == Interval(intersectionat(0.3), intersectionat(0.4)) && res[3] == positivehalfspace(intersectionat(0.7))
            end

            @testcase begin
                a = DisjointUnion(rayorigininterval(intersectionat(0.2)), Interval(intersectionat(0.4), intersectionat(0.7)))
                res = intervalcomplement(a)
                @test res[1] == Interval(intersectionat(0.2), intersectionat(0.4)) && res[2] == positivehalfspace(intersectionat(0.7))
            end
        end
    end # testset interval
end # testset General
