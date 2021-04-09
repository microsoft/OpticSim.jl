@testset "Intersection" begin
    function samplepoints(numsamples, lowu, highu, lowv, highv)
        samples = Array{Tuple{Float64,Float64},1}(undef, 0)

        for i in 0:numsamples
            u = lowu * i / numsamples + (1 - i / numsamples) * highu
            for j in 0:numsamples
                v = lowv * j / numsamples + (1 - j / numsamples) * highv
                push!(samples, (u, v))
            end
        end
        return samples
    end

    @testset "Cylinder" begin
        # Random.seed!(SEED)

        # # random point intersections
        # ur, vr = uvrange(Cylinder)
        # for i in 1:10000
        #     cyl = Cylinder(rand() * 3.0, rand() * 50.0)
        #     u, v = rand() * (ur[2] - ur[1]) + ur[1], rand() * (vr[2] - vr[1]) + vr[1]
        #     pointon = point(cyl, u, v)
        #     origin = SVector(4.0, 4.0, 4.0)
        #     r = Ray(origin, pointon .- origin)
        #     halfspace = surfaceintersection(cyl, r)
        #     @test halfspace !== nothing && !((lower(halfspace) isa RayOrigin) || (upper(halfspace) isa Infinity))
        #     @test isapprox(pointon, point(lower(halfspace)), rtol = RTOLERANCE, atol = ATOLERANCE) || isapprox(pointon, point(upper(halfspace)), rtol = RTOLERANCE, atol = ATOLERANCE)
        # end

        # on xy
        cyl = Cylinder(0.5)
        r = Ray([1.0, 1.0, 0.0], [-1.0, -1.0, 0.0])
        intscts = surfaceintersection(cyl, r)
        xy = sqrt(2) / 4.0
        pt1 = [xy, xy, 0.0]
        pt2 = [-xy, -xy, 0.0]
        cpt1 = point(lower(intscts))
        cpt2 = point(upper(intscts))
        @test isapprox(pt1, cpt1, rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(pt2, cpt2, rtol = RTOLERANCE, atol = ATOLERANCE)

        # starting inside
        r = Ray([0.0, 0.0, 0.0], [1.0, 1.0, 0.0])
        intscts = surfaceintersection(cyl, r)
        @test lower(intscts) isa RayOrigin && isapprox(pt1, point(upper(intscts)), rtol = RTOLERANCE, atol = ATOLERANCE)

        # inside infinite
        r = Ray([0.0, 0.0, 0.0], [0.0, 0.0, 1.0])
        intscts = surfaceintersection(cyl, r)
        @test lower(intscts) isa RayOrigin && upper(intscts) isa Infinity

        # on surface
        r = Ray([0.5, 0.0, 0.0], [0.0, 0.0, 1.0])
        intscts = surfaceintersection(cyl, r)
        @test lower(intscts) isa RayOrigin && upper(intscts) isa Infinity

        # on diagonal
        r = Ray([1.0, 1.0, 1.0], [-1.0, -1.0, -1.0])
        intscts = surfaceintersection(cyl, r)
        cpt1 = point(lower(intscts))
        cpt2 = point(upper(intscts))
        pt1 = [xy, xy, xy]
        pt2 = [-xy, -xy, -xy]
        @test isapprox(pt1, cpt1, rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(pt2, cpt2, rtol = RTOLERANCE, atol = ATOLERANCE)

        # missing
        r = Ray([5.0, 5.0, 5.0], [0.0, 0.0, -1.0])
        intscts = surfaceintersection(cyl, r)
        @test intscts isa EmptyInterval
    end # testset cylinder

    @testset "Sphere" begin
        # Random.seed!(SEED)

        # # random point intersections
        # ur, vr = uvrange(Sphere)
        # for i in 1:10000
        #     sph = Sphere(rand() * 3)
        #     u, v = rand() * (ur[2] - ur[1]) + ur[1], rand() * (vr[2] - vr[1]) + vr[1]
        #     pointon = point(sph, u, v)
        #     origin = SVector(4.0, 4.0, 4.0)
        #     r = Ray(origin, pointon .- origin)
        #     halfspace = surfaceintersection(sph, r)
        #     @test halfspace !== nothing && !((lower(halfspace) isa RayOrigin) || (upper(halfspace) isa Infinity))
        #     @test isapprox(pointon, point(lower(halfspace)), rtol = RTOLERANCE, atol = ATOLERANCE) || isapprox(pointon, point(upper(halfspace)), rtol = RTOLERANCE, atol = ATOLERANCE)
        # end

        # on xy plane
        sph = Sphere(0.5)
        r = Ray([1.0, 1.0, 0.0], [-1.0, -1.0, 0.0])
        intscts = surfaceintersection(sph, r)
        xy = sqrt(2) / 4.0
        pt1 = [xy, xy, 0.0]
        pt2 = [-xy, -xy, 0.0]
        cpt1 = point(lower(intscts))
        cpt2 = point(upper(intscts))
        @test isapprox(pt1, cpt1, rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(pt2, cpt2, rtol = RTOLERANCE, atol = ATOLERANCE)

        # starting inside
        r = Ray([0.0, 0.0, 0.0], [1.0, 1.0, 0.0])
        intscts = surfaceintersection(sph, r)
        @test lower(intscts) isa RayOrigin && isapprox(pt1, point(upper(intscts)), rtol = RTOLERANCE, atol = ATOLERANCE)

        # on diagonal
        r = Ray([1.0, 1.0, 1.0], [-1.0, -1.0, -1.0])
        intscts = surfaceintersection(sph, r)
        cpt1 = point(lower(intscts))
        cpt2 = point(upper(intscts))
        xyz = sqrt(3) / 6.0
        pt1 = [xyz, xyz, xyz]
        pt2 = [-xyz, -xyz, -xyz]
        @test isapprox(pt1, cpt1, rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(pt2, cpt2, rtol = RTOLERANCE, atol = ATOLERANCE)

        # missing
        r = Ray([5.0, 5.0, 5.0], [0.0, 0.0, -1.0])
        intscts = surfaceintersection(sph, r)
        @test intscts isa EmptyInterval

        # tangent
        r = Ray([0.5, 0.5, 0.0], [0.0, -1.0, 0.0])
        intscts = surfaceintersection(sph, r)
        @test intscts isa EmptyInterval
    end # testset sphere

    @testset "Spherical Cap" begin
        @test_throws AssertionError SphericalCap(0.0, 1.0)
        @test_throws AssertionError SphericalCap(1.0, 0.0)

        sph = SphericalCap(0.5, π / 3, SVector(1.0, 1.0, 0.0), SVector(1.0, 1.0, 1.0))

        r = Ray([10.0, 10.0, 1.0], [-1.0, -1.0, 0.0])
        int = surfaceintersection(sph, r)
        @test isapprox(point(lower(int)), [1.0, 1.0, 1.0], rtol = RTOLERANCE, atol = ATOLERANCE) && upper(int) isa Infinity

        r = Ray([0.8, 1.4, 1.2], [1.0, -1.0, 0.0])
        int = surfaceintersection(sph, r)
        @test isapprox(point(lower(int)), [0.898231126982741, 1.301768873017259, 1.2], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(int)), [1.301768873017259, 0.8982311269827411, 1.2], rtol = RTOLERANCE, atol = ATOLERANCE)

        r = Ray([2.0, 2.0, 1.5], [-1.0, -1.0, 0.0])
        @test surfaceintersection(sph, r) isa EmptyInterval
    end # testset spherical cap

    @testset "Triangle" begin
        Random.seed!(SEED)

        tri = Triangle(SVector{3}(2.0, 1.0, 1.0), SVector{3}(1.0, 2.0, 1.0), SVector{3}(1.0, 1.0, 2.0), SVector{2}(0.0, 0.0), SVector{2}(1.0, 0.0), SVector{2}(0.5, 0.5))

        # random point intersections
        for i in 1:1000
            barycentriccoords = rand(3)
            barycentriccoords /= sum(barycentriccoords)
            pointontri = point(tri, barycentriccoords...)
            origin = SVector(3.0, 3.0, 3.0)
            r = Ray(origin, pointontri .- origin)
            halfspace = surfaceintersection(tri, r)
            @test !(halfspace isa EmptyInterval) && upper(halfspace) isa Infinity
            t = α(lower(halfspace))
            @test isapprox(point(r, t), point(lower(halfspace)), rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(pointontri, point(lower(halfspace)), rtol = RTOLERANCE, atol = ATOLERANCE)
        end

        # front and back face intersections
        tri = Triangle(SVector{3}(0.0, 0.0, 0.0), SVector{3}(1.0, 0.0, 0.0), SVector{3}(0.0, 1.0, 0.0), SVector{2}(0.0, 0.0), SVector{2}(1.0, 0.0), SVector{2}(0.0, 1.0))
        r1 = Ray([0.1, 0.1, 1.0], [0.0, 0.0, -1.0])
        r2 = Ray([0.1, 0.1, -1.0], [0.0, 0.0, 1.0])

        intsct1 = halfspaceintersection(surfaceintersection(tri, r1))
        intsct2 = halfspaceintersection(surfaceintersection(tri, r2))

        @test isapprox(point(intsct1), point(intsct2), rtol = RTOLERANCE, atol = ATOLERANCE)

        # ray should miss
        r = Ray([1.0, 1.0, 1.0], [0.0, 0.0, -1.0])
        intsct = surfaceintersection(tri, r)
        @test intsct isa EmptyInterval
    end # testset triangle

    @testset "Plane" begin
        rinplane = Ray([0.0, 0.0, 1.0], [1.0, 1.0, 0.0])
        routside = Ray([0.0, 0.0, 2.0], [1.0, 1.0, 1.0])
        rinside = Ray([0.0, 0.0, -1.0], [1.0, 1.0, -1.0])
        routintersects = Ray([0.0, 0.0, 2.0], [0.0, 0.0, -1.0])
        rinintersects = Ray([0.0, 0.0, -1.0], [0.0, 0.0, 1.0])
        pln = Plane(0.0, 0.0, 1.0, 0.0, 0.0, 1.0)

        # parallel to plane and outside
        @test surfaceintersection(pln, routside) isa EmptyInterval

        # NOTE for coplanar faces to work with visualization, rays in the plane must count as being 'inside'
        # parallel and on plane
        res = surfaceintersection(pln, rinplane)
        @test lower(res) isa RayOrigin && upper(res) isa Infinity

        # inside but not hitting
        res = surfaceintersection(pln, rinside)
        @test lower(res) isa RayOrigin && upper(res) isa Infinity

        # starts outside and hits
        res = surfaceintersection(pln, routintersects)
        @test isapprox(point(lower(res)), [0.0, 0.0, 1.0], rtol = RTOLERANCE, atol = ATOLERANCE) && upper(res) isa Infinity

        # starts inside and hits
        res = surfaceintersection(pln, rinintersects)
        @test lower(res) isa RayOrigin && isapprox(point(upper(res)), [0.0, 0.0, 1.0], rtol = RTOLERANCE, atol = ATOLERANCE)
    end # testset plane

    @testset "Rectangle" begin
        @test_throws AssertionError Rectangle(0.0, 1.0)
        @test_throws AssertionError Rectangle(1.0, 0.0)

        rinplane = Ray([0.0, 0.0, 0.0], [1.0, 1.0, 0.0])
        routside = Ray([0.0, 0.0, 1.0], [1.0, 1.0, 1.0])
        rinside = Ray([0.0, 0.0, -1.0], [1.0, 1.0, -1.0])
        routintersects = Ray([0.0, 0.0, 1.0], [0.0, 0.0, -1.0])
        rinintersects = Ray([0.0, 0.0, -1.0], [0.0, 0.0, 1.0])
        routmiss = Ray([0.0, 2.0, 1.0], [0.0, 0.0, -1.0])
        rinmiss = Ray([0.0, 2.0, -1.0], [0.0, 0.0, 1.0])
        rinbounds = Ray([0.5, 0.5, -1.0], [0.0, 0.0, 1.0])
        routbounds = Ray([0.5, 0.5, 1.0], [0.0, 0.0, -1.0])

        rect = Rectangle(0.5, 0.5)

        # parallel to plane and outside
        @test surfaceintersection(rect, routside) isa EmptyInterval

        # parallel and on plane
        @test surfaceintersection(rect, rinplane) isa EmptyInterval

        # inside but not hitting
        @test surfaceintersection(rect, rinside) isa EmptyInterval

        # starts outside and hits
        res = surfaceintersection(rect, routintersects)
        @test isapprox(point(lower(res)), [0.0, 0.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE) && upper(res) isa Infinity

        # starts inside and hits
        res = surfaceintersection(rect, rinintersects)
        @test lower(res) isa RayOrigin && isapprox(point(upper(res)), [0.0, 0.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE)

        # starts inside and misses bounds
        @test surfaceintersection(rect, rinmiss) isa EmptyInterval

        # starts outside and misses bounds
        @test surfaceintersection(rect, routmiss) isa EmptyInterval

        # starts outside and hits bounds
        res = surfaceintersection(rect, routbounds)
        @test isapprox(point(lower(res)), [0.5, 0.5, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE) && upper(res) isa Infinity

        # starts inside and hits bounds
        res = surfaceintersection(rect, rinbounds)
        @test lower(res) isa RayOrigin && isapprox(point(upper(res)), [0.5, 0.5, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE)

        # test a rectangle with translation and rotation
        rect2 = Rectangle(0.3, 0.5, SVector(3.0, 1.0, 4.0), SVector(0.2, 0.3, 0.4))
        rayhit = Ray([0.6, 0.6, 0.6], [-0.3, -0.1, -0.3])
        raymiss = Ray([0.7, 0.4, 0.4], [-0.3, -0.1, -0.3])
        res = surfaceintersection(rect2, rayhit)
        @test isapprox(point(lower(res)), [0.2863636363636363, 0.4954545454545454, 0.2863636363636363], rtol = RTOLERANCE, atol = ATOLERANCE) && upper(res) isa Infinity
        @test surfaceintersection(rect2, raymiss) isa EmptyInterval
    end # testset Rectangle

    @testset "Ellipse" begin
        @test_throws AssertionError Ellipse(0.0, 1.0)
        @test_nowarn Circle(0.5)

        rinplane = Ray([0.0, 0.0, 0.0], [1.0, 1.0, 0.0])
        routside = Ray([0.0, 0.0, 1.0], [1.0, 1.0, 1.0])
        rinside = Ray([0.0, 0.0, -1.0], [1.0, 1.0, -1.0])
        routintersects = Ray([0.0, 0.0, 1.0], [0.0, 0.0, -1.0])
        rinintersects = Ray([0.0, 0.0, -1.0], [0.0, 0.0, 1.0])
        routmiss = Ray([0.0, 2.0, 1.0], [0.0, 0.0, -1.0])
        rinmiss = Ray([0.0, 2.0, -1.0], [0.0, 0.0, 1.0])
        rinbounds = Ray([0.5, 0.0, -1.0], [0.0, 0.0, 1.0])
        routbounds = Ray([0.5, 0.0, 1.0], [0.0, 0.0, -1.0])
        rasymhit = Ray([0.0, 0.9, 1.0], [0.0, 0.0, -1.0])
        rasymmiss = Ray([0.9, 0.0, 1.0], [0.0, 0.0, -1.0])

        ell = Ellipse(0.5, 1.0)

        # parallel to plane and outside
        @test surfaceintersection(ell, routside) isa EmptyInterval

        # parallel and on plane
        @test surfaceintersection(ell, rinplane) isa EmptyInterval

        # inside but not hitting
        @test surfaceintersection(ell, rinside) isa EmptyInterval

        # starts outside and hits
        res = surfaceintersection(ell, routintersects)
        @test isapprox(point(lower(res)), [0.0, 0.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE) && upper(res) isa Infinity

        # starts inside and hits
        res = surfaceintersection(ell, rinintersects)
        @test lower(res) isa RayOrigin && isapprox(point(upper(res)), [0.0, 0.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE)

        # starts inside and misses bounds
        @test surfaceintersection(ell, rinmiss) isa EmptyInterval

        # starts outside and misses bounds
        @test surfaceintersection(ell, routmiss) isa EmptyInterval

        # starts outside and hits bounds
        res = surfaceintersection(ell, routbounds)
        @test isapprox(point(lower(res)), [0.5, 0.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE) && upper(res) isa Infinity

        # starts inside and hits bounds
        res = surfaceintersection(ell, rinbounds)
        @test lower(res) isa RayOrigin && isapprox(point(upper(res)), [0.5, 0.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE)

        # asymmetric hit
        res = surfaceintersection(ell, rasymhit)
        @test isapprox(point(lower(res)), [0.0, 0.9, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE) && upper(res) isa Infinity

        # asymmetric miss
        @test surfaceintersection(ell, rasymmiss) isa EmptyInterval

        # test an ellipse with translation and rotation
        ell2 = Ellipse(0.3, 0.5, SVector(3.0, 1.0, 4.0), SVector(0.2, 0.3, 0.4))
        rayhit = Ray([0.6, 0.6, 0.6], [-0.3, -0.1, -0.3])
        raymiss = Ray([0.7, 0.4, 0.4], [-0.3, -0.1, -0.3])
        res = surfaceintersection(ell2, rayhit)
        @test isapprox(point(lower(res)), [0.2863636363636363, 0.4954545454545454, 0.2863636363636363], rtol = RTOLERANCE, atol = ATOLERANCE) && upper(res) isa Infinity
        @test surfaceintersection(ell2, raymiss) isa EmptyInterval
    end # testset Ellipse

    @testset "Hexagon" begin
        @test_nowarn Hexagon(0.5)

        rinplane = Ray([0.0, 0.0, 0.0], [1.0, 1.0, 0.0])
        routside = Ray([0.0, 0.0, 1.0], [1.0, 1.0, 1.0])
        rinside = Ray([0.0, 0.0, -1.0], [1.0, 1.0, -1.0])
        routintersects = Ray([0.0, 0.0, 1.0], [0.0, 0.0, -1.0])
        rinintersects = Ray([0.0, 0.0, -1.0], [0.0, 0.0, 1.0])
        routmiss = Ray([0.0, 2.0, 1.0], [0.0, 0.0, -1.0])
        rinmiss = Ray([0.0, 2.0, -1.0], [0.0, 0.0, 1.0])
        rinbounds = Ray([sqrt(3) / 2, 0.0, -1.0], [0.0, 0.0, 1.0])
        routbounds = Ray([sqrt(3) / 2, 0.0, 1.0], [0.0, 0.0, -1.0])
        rasymhit = Ray([0.0, 1.0, 1.0], [0.0, 0.0, -1.0])
        rasymmiss = Ray([1.0, 0.0, 1.0], [0.0, 0.0, -1.0])

        hex = Hexagon(1.0)

        # parallel to plane and outside
        @test surfaceintersection(hex, routside) isa EmptyInterval

        # parallel and on plane
        @test surfaceintersection(hex, rinplane) isa EmptyInterval

        # inside but not hitting
        @test surfaceintersection(hex, rinside) isa EmptyInterval

        # starts outside and hits
        res = surfaceintersection(hex, routintersects)
        @test isapprox(point(lower(res)), [0.0, 0.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE) && upper(res) isa Infinity

        # starts inside and hits
        res = surfaceintersection(hex, rinintersects)
        @test lower(res) isa RayOrigin && isapprox(point(upper(res)), [0.0, 0.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE)

        # starts inside and misses bounds
        @test surfaceintersection(hex, rinmiss) isa EmptyInterval

        # starts outside and misses bounds
        @test surfaceintersection(hex, routmiss) isa EmptyInterval

        # starts outside and hits bounds
        res = surfaceintersection(hex, routbounds)
        @test isapprox(point(lower(res)), [sqrt(3) / 2, 0.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE) && upper(res) isa Infinity

        # starts inside and hits bounds
        res = surfaceintersection(hex, rinbounds)
        @test lower(res) isa RayOrigin && isapprox(point(upper(res)), [sqrt(3) / 2, 0.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE)

        # asymmetric hit
        res = surfaceintersection(hex, rasymhit)
        @test isapprox(point(lower(res)), [0.0, 1.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE) && upper(res) isa Infinity

        # asymmetric miss
        @test surfaceintersection(hex, rasymmiss) isa EmptyInterval

        # test an ellipse with translation and rotation
        hex2 = Hexagon(0.4, SVector(3.0, 1.0, 4.0), SVector(0.2, 0.3, 0.4))
        rayhit = Ray([0.6, 0.6, 0.6], [-0.3, -0.1, -0.3])
        raymiss = Ray([0.7, 0.4, 0.4], [-0.3, -0.1, -0.3])
        res = surfaceintersection(hex2, rayhit)
        @test isapprox(point(lower(res)), [0.2863636363636363, 0.4954545454545454, 0.2863636363636363], rtol = RTOLERANCE, atol = ATOLERANCE) && upper(res) isa Infinity
        @test surfaceintersection(hex2, raymiss) isa EmptyInterval
    end # testset Hexagon

    @testset "Stops" begin
        infiniterect = RectangularAperture(0.4, 0.8, SVector(0.0, 1.0, 0.0), SVector(0.0, 0.0, 1.0))
        finiterect = RectangularAperture(0.4, 0.8, 1.0, 2.0, SVector(0.0, 1.0, 0.0), SVector(0.0, 0.0, 1.0))
        # through hole
        r = Ray([0.0, 1.0, 1.0], [0.0, -1.0, 0.0])
        @test surfaceintersection(infiniterect, r) isa EmptyInterval
        @test surfaceintersection(finiterect, r) isa EmptyInterval
        r = Ray([0.0, -1.0, 1.0], [0.0, 1.0, 0.0])
        @test surfaceintersection(infiniterect, r) isa EmptyInterval
        @test surfaceintersection(finiterect, r) isa EmptyInterval
        # on edge
        r = Ray([0.0, 1.0, 0.6], [0.0, -1.0, 0.0])
        res = surfaceintersection(infiniterect, r)
        @test isapprox(point(lower(res)), [0.0, 0.0, 0.6], rtol = RTOLERANCE, atol = ATOLERANCE) && upper(res) isa Infinity
        res = surfaceintersection(finiterect, r)
        @test isapprox(point(lower(res)), [0.0, 0.0, 0.6], rtol = RTOLERANCE, atol = ATOLERANCE) && upper(res) isa Infinity
        # through hole asym
        r = Ray([0.7, 1.0, 1.0], [0.0, -1.0, 0.0])
        @test surfaceintersection(infiniterect, r) isa EmptyInterval
        @test surfaceintersection(finiterect, r) isa EmptyInterval
        # on finite edge
        r = Ray([1.0, -1.0, 2.0], [0.0, 1.0, 0.0])
        res = surfaceintersection(infiniterect, r)
        @test lower(res) isa RayOrigin && isapprox(point(upper(res)), [1.0, 0.0, 2.0], rtol = RTOLERANCE, atol = ATOLERANCE)
        res = surfaceintersection(finiterect, r)
        @test lower(res) isa RayOrigin && isapprox(point(upper(res)), [1.0, 0.0, 2.0], rtol = RTOLERANCE, atol = ATOLERANCE)
        # outside finite
        r = Ray([2.0, 1.0, 3.0], [0.0, -1.0, 0.0])
        @test isapprox(point(surfaceintersection(infiniterect, r)), [2.0, 0.0, 3.0], rtol = RTOLERANCE, atol = ATOLERANCE)
        @test surfaceintersection(finiterect, r) isa EmptyInterval

        infinitecirc = CircularAperture(0.4, SVector(0.0, 1.0, 0.0), SVector(0.0, 0.0, 1.0))
        finitecirc = CircularAperture(0.4, 1.0, 2.0, SVector(0.0, 1.0, 0.0), SVector(0.0, 0.0, 1.0))
        # through hole
        r = Ray([0.0, 1.0, 1.0], [0.0, -1.0, 0.0])
        @test surfaceintersection(infinitecirc, r) isa EmptyInterval
        @test surfaceintersection(finitecirc, r) isa EmptyInterval
        r = Ray([0.0, -1.0, 1.0], [0.0, 1.0, 0.0])
        @test surfaceintersection(infinitecirc, r) isa EmptyInterval
        @test surfaceintersection(finitecirc, r) isa EmptyInterval
        # on edge
        r = Ray([0.0, 1.0, 0.6], [0.0, -1.0, 0.0])
        res = surfaceintersection(infinitecirc, r)
        @test isapprox(point(lower(res)), [0.0, 0.0, 0.6], rtol = RTOLERANCE, atol = ATOLERANCE) && upper(res) isa Infinity
        res = surfaceintersection(finitecirc, r)
        @test isapprox(point(lower(res)), [0.0, 0.0, 0.6], rtol = RTOLERANCE, atol = ATOLERANCE) && upper(res) isa Infinity
        # through hole 2
        r = Ray([0.3, 1.0, 1.0], [0.0, -1.0, 0.0])
        @test surfaceintersection(infinitecirc, r) isa EmptyInterval
        @test surfaceintersection(finitecirc, r) isa EmptyInterval
        # on finite edge
        r = Ray([1.0, -1.0, 2.0], [0.0, 1.0, 0.0])
        res = surfaceintersection(infinitecirc, r)
        @test lower(res) isa RayOrigin && isapprox(point(upper(res)), [1.0, 0.0, 2.0], rtol = RTOLERANCE, atol = ATOLERANCE)
        res = surfaceintersection(finitecirc, r)
        @test lower(res) isa RayOrigin && isapprox(point(upper(res)), [1.0, 0.0, 2.0], rtol = RTOLERANCE, atol = ATOLERANCE)
        # outside finite
        r = Ray([2.0, 1.0, 3.0], [0.0, -1.0, 0.0])
        @test isapprox(point(surfaceintersection(infinitecirc, r)), [2.0, 0.0, 3.0], rtol = RTOLERANCE, atol = ATOLERANCE)
        @test surfaceintersection(finitecirc, r) isa EmptyInterval

        annulus = Annulus(0.5, 1.0, SVector(0.0, 1.0, 0.0), SVector(0.0, 0.0, 1.0))
        # through hole
        r = Ray([0.0, 1.0, 1.0], [0.0, -1.0, 0.0])
        @test surfaceintersection(annulus, r) isa EmptyInterval
        r = Ray([0.0, -1.0, 1.0], [0.0, 1.0, 0.0])
        @test surfaceintersection(annulus, r) isa EmptyInterval
        # on edge
        r = Ray([0.0, 1.0, 0.5], [0.0, -1.0, 0.0])
        res = surfaceintersection(annulus, r)
        @test isapprox(point(lower(res)), [0.0, 0.0, 0.5], rtol = RTOLERANCE, atol = ATOLERANCE) && upper(res) isa Infinity
        # through hole 2
        r = Ray([0.3, 1.0, 1.0], [0.0, -1.0, 0.0])
        @test surfaceintersection(infiniterect, r) isa EmptyInterval
        @test surfaceintersection(annulus, r) isa EmptyInterval
        # on finite edge
        r = Ray([1.0, -1.0, 1.0], [0.0, 1.0, 0.0])
        res = surfaceintersection(annulus, r)
        @test lower(res) isa RayOrigin && isapprox(point(upper(res)), [1.0, 0.0, 1.0], rtol = RTOLERANCE, atol = ATOLERANCE)
        # outside finite
        r = Ray([0.9, 1.0, 1.9], [0.0, -1.0, 0.0])
        @test surfaceintersection(annulus, r) isa EmptyInterval
    end # testset Stops

    @testset "Bezier" begin
        surf = TestData.beziersurface()
        accelsurf = AcceleratedParametricSurface(surf)
        numsamples = 100

        samples = samplepoints(numsamples, 0.05, 0.95, 0.05, 0.95)
        missedintersections = 0

        for uv in samples
            pt = point(surf, uv[1], uv[2])
            origin = [0.5, 0.5, 5.0]
            r = Ray(origin, pt .- origin)
            allintersections = surfaceintersection(accelsurf, r)
            if allintersections isa EmptyInterval
                # sampled triangle acceleration structure doesn't guarantee all intersecting rays will be detected.
                # Have to use patch subdivision and convex hull polyhedron to do this.
                missedintersections += 1
            else
                if isa(allintersections, Interval)
                    allintersections = (allintersections,)
                end

                for intsct in allintersections
                    @test isapprox(point(halfspaceintersection(intsct)), pt, rtol = RTOLERANCE, atol = ATOLERANCE)
                    @test upper(intsct) isa Infinity
                end
            end
        end

        if missedintersections > 0
            @warn "$missedintersections out of total $(length(samples)) bezier suface intersections were missed"
        end

        # miss from outside
        r = Ray([0.5, 0.5, 5.0], [0.0, 0.0, 1.0])
        @test surfaceintersection(accelsurf, r) isa EmptyInterval
        r = Ray([5.0, 0.5, -5.0], [0.0, 0.0, -1.0])
        @test surfaceintersection(accelsurf, r) isa EmptyInterval

        # miss from inside
        r = Ray([0.5, 0.5, -5.0], [0.0, 0.0, -1.0])
        res = surfaceintersection(accelsurf, r)
        # TODO!! Fix bezier surface to create halfspace
        # @test lower(res) isa RayOrigin && upper(res) isa Infinity

        # hit from inside
        r = Ray([0.5, 0.5, 0.2], [0.0, 1.0, 0.0])
        res = surfaceintersection(accelsurf, r)
        @test lower(res) isa RayOrigin && isapprox(point(upper(res)), [0.5, 0.8996900152986361, 0.2], rtol = RTOLERANCE, atol = ATOLERANCE)

        # hit from outside
        r = Ray([0.0, 0.5, 0.2], [1.0, 0.0, -1.0])
        res = surfaceintersection(accelsurf, r)
        @test isapprox(point(lower(res)), [0.06398711204047353, 0.5, 0.13601288795952649], rtol = RTOLERANCE, atol = ATOLERANCE) && upper(res) isa Infinity

        # two hits from outside
        r = Ray([0.0, 0.5, 0.25], [1.0, 0.0, 0.0])
        res = surfaceintersection(accelsurf, r)
        @test isapprox(point(lower(res)), [0.12607777053030267, 0.5, 0.25], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(res)), [0.8705887077060419, 0.5, 0.25], rtol = RTOLERANCE, atol = ATOLERANCE)

        surf = TestData.wavybeziersurface()
        accelsurf = AcceleratedParametricSurface(surf)

        # 3 hit starting outside
        r = Ray([0.5, 0.0, 0.0], [0.0, 1.0, 0.0])
        res = surfaceintersection(accelsurf, r)
        @test isa(res, DisjointUnion) && length(res) == 2 && isapprox(point(lower(res[1])), [0.5, 0.10013694786182059, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(res[1])), [0.5, 0.49625, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(lower(res[2])), [0.5, 0.8971357794109067, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE) && (upper(res[2]) isa Infinity)

        # 3 hit starting inside
        r = Ray([0.5, 1.0, 0.0], [0.0, -1.0, 0.0])
        res = surfaceintersection(accelsurf, r)
        @test isa(res, DisjointUnion) && length(res) == 2 && (lower(res[1]) isa RayOrigin) && isapprox(point(upper(res[1])), [0.5, 0.8971357794109067, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(lower(res[2])), [0.5, 0.49625, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(res[2])), [0.5, 0.10013694786182059, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE)

        surf = TestData.verywavybeziersurface()
        accelsurf = AcceleratedParametricSurface(surf, 20)

        # five hits starting outside
        r = Ray([0.9, 0.0, -0.3], [0.0, 1.0, 0.7])
        res = surfaceintersection(accelsurf, r)
        a = isapprox(point(lower(res[1])), [0.9, 0.03172286522032046, -0.2777939943457758], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(res[1])), [0.9, 0.1733979947040411, -0.17862140370717122], rtol = RTOLERANCE, atol = ATOLERANCE)
        b = isapprox(point(lower(res[2])), [0.9, 0.5335760974594397, 0.07350326822160776], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(res[2])), [0.9, 0.7767707607392784, 0.24373953251749486], rtol = RTOLERANCE, atol = ATOLERANCE)
        c = isapprox(point(lower(res[3])), [0.9, 0.9830891958374246, 0.3881624370861975], rtol = RTOLERANCE, atol = ATOLERANCE) && (upper(res[3]) isa Infinity)
        @test isa(res, DisjointUnion) && length(res) == 3 && a && b && c

        # five hits starting inside
        r = Ray([0.9, 1.0, 0.4], [0.0, -1.0, -0.7])
        res = surfaceintersection(accelsurf, r)
        a = (lower(res[1]) isa RayOrigin) && isapprox(point(upper(res[1])), [0.9, 0.9830891958374246, 0.3881624370861975], rtol = RTOLERANCE, atol = ATOLERANCE)
        b = isapprox(point(lower(res[2])), [0.9, 0.7767707607392784, 0.24373953251749486], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(res[2])), [0.9, 0.5335760974594397, 0.07350326822160776], rtol = RTOLERANCE, atol = ATOLERANCE)
        c = isapprox(point(lower(res[3])), [0.9, 0.17339799470404108, -0.1786214037071712], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(res[3])), [0.9, 0.03172286522032046, -0.27779399434577573], rtol = RTOLERANCE, atol = ATOLERANCE)
        @test isa(res, DisjointUnion) && length(res) == 3 && a && b && c

        # 4 hits starting inside
        r = Ray([0.1, 0.0, -0.3], [0.0, 1.0, 0.7])
        res = surfaceintersection(accelsurf, r)
        a = (lower(res[1]) isa RayOrigin) && isapprox(point(upper(res[1])), [0.1, 0.2851860296285551, -0.10036977926001144], rtol = RTOLERANCE, atol = ATOLERANCE)
        b = isapprox(point(lower(res[2])), [0.1, 0.5166793625025807, 0.06167555375180668], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(res[2])), [0.1, 0.7770862508789854, 0.24396037561528983], rtol = RTOLERANCE, atol = ATOLERANCE)
        c = isapprox(point(lower(res[3])), [0.1, 0.98308919558696, 0.3881624369108719], rtol = RTOLERANCE, atol = ATOLERANCE) && (upper(res[3]) isa Infinity)
        @test isa(res, DisjointUnion) && length(res) == 3 && a && b && c

        # 4 hits starting outside
        r = Ray([0.9, 0.9, 0.4], [0.0, -0.9, -0.7])
        res = surfaceintersection(accelsurf, r)
        a = isapprox(point(lower(res[1])), [0.9, 0.736072142615238, 0.2725005553674076], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(res[1])), [0.9, 0.567439326091764, 0.141341698071372], rtol = RTOLERANCE, atol = ATOLERANCE)
        b = isapprox(point(lower(res[2])), [0.9, 0.16601081959179267, -0.1708804736508277], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(res[2])), [0.9, 0.032434058775915924, -0.274773509840954], rtol = RTOLERANCE, atol = ATOLERANCE)
        @test isa(res, DisjointUnion) && length(res) == 2 && a && b
    end # testset Bezier

    @testset "Zernike" begin
        @test_nowarn ZernikeSurface(1.5)
        @test_throws AssertionError ZernikeSurface(0.0)
        @test_throws AssertionError ZernikeSurface(1.5, aspherics = [(1, 0.1)])

        z1 = TestData.zernikesurface1()
        az1 = AcceleratedParametricSurface(z1, 20)
        numsamples = 100

        # test that an on axis ray doesn't miss
        @test !(surfaceintersection(AcceleratedParametricSurface(z1), Ray([0.0, 0.0, 10.0], [0.0, 0.0, -1.0])) isa EmptyInterval)

        samples = samplepoints(numsamples, 0.01, 0.99, 0.01, 2π - 0.01)
        missedintersections = 0

        for uv in samples
            pt = point(z1, uv[1], uv[2])
            origin = [0.0, 0.0, 5.0]
            r = Ray(origin, pt .- origin)
            allintersections = surfaceintersection(az1, r)
            if (allintersections isa EmptyInterval)
                # sampled triangle acceleration structure doesn't guarantee all intersecting rays will be detected.
                # Have to use patch subdivision and convex hull polyhedron to do this.
                missedintersections += 1
            else
                if isa(allintersections, Interval)
                    allintersections = (allintersections,)
                end

                for intsct in allintersections
                    if lower(intsct) isa RayOrigin
                        missedintersections += 1
                    else
                        # closest point should be on the surface, furthest should be on bounding prism
                        @test isapprox(point(lower(intsct)), pt, rtol = RTOLERANCE, atol = ATOLERANCE)
                        @test isa(upper(intsct), Intersection{Float64,3}) || (OpticSim.direction(r)[3] == -1 && isa(upper(intsct), Infinity{Float64}))
                    end
                end
            end
        end

        if missedintersections > 0
            @warn "$missedintersections out of total $(length(samples)) zernike suface intersections were missed"
        end

        # hit from inside
        r = Ray([0.0, 0.0, -0.5], [0.1, 0.2, 0.5])
        intvl = surfaceintersection(az1, r)
        @test (lower(intvl) isa RayOrigin) && isapprox(point(upper(intvl)), [0.12363711269619936, 0.24727422539239863, 0.11818556348099664], rtol = RTOLERANCE, atol = ATOLERANCE)

        # hit from outside
        r = Ray([0.0, 0.0, 0.5], [0.1, 0.2, -0.5])
        intvl = surfaceintersection(az1, r)
        @test isapprox(point(lower(intvl)), [0.07118311042701463, 0.1423662208540292, 0.144084447864927], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(intvl)), [0.6708203932499369, 1.3416407864998738, -2.8541019662496843], rtol = RTOLERANCE, atol = ATOLERANCE)

        # miss from inside
        r = Ray([0.2, 0.2, -0.5], [0.0, 0.0, -1.0])
        intvl = surfaceintersection(az1, r)
        @test (lower(intvl) isa RayOrigin) && (upper(intvl) isa Infinity)

        # miss from outside
        r = Ray([0.2, 2.0, -0.5], [0.0, 0.0, 1.0])
        @test surfaceintersection(az1, r) isa EmptyInterval
        r = Ray([0.2, 2.0, 0.5], [0.0, 0.0, -1.0])
        @test surfaceintersection(az1, r) isa EmptyInterval
        r = Ray([0.2, 0.2, 0.5], [0.0, 0.0, 1.0])
        @test surfaceintersection(az1, r) isa EmptyInterval
        r = Ray([0.2, 2.0, 2.0], [0.0, -1.0, 0.0])
        @test surfaceintersection(az1, r) isa EmptyInterval

        # two hits from outside
        r = Ray([2.0, 0.0, 0.25], [-1.0, 0.0, 0.0])
        hit = surfaceintersection(az1, r)
        a = isapprox(point(lower(hit[1])), [1.5, 0.0, 0.25], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(hit[1])), [1.3571758210851095, 0.0, 0.25], rtol = RTOLERANCE, atol = ATOLERANCE)
        b = isapprox(point(lower(hit[2])), [-1.2665828947521165, 0.0, 0.25], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(hit[2])), [-1.5, 0.0, 0.25], rtol = RTOLERANCE, atol = ATOLERANCE)
        @test isa(hit, DisjointUnion) && length(hit) == 2 && a && b

        # hit cyl from outside
        r = Ray([0.0, 2.0, -0.5], [0.0, -1.0, 0.0])
        intvl = surfaceintersection(az1, r)
        @test isapprox(point(lower(intvl)), [0.0, 1.5, -0.5], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(intvl)), [0.0, -1.5, -0.5], rtol = RTOLERANCE, atol = ATOLERANCE)

        # hit cyl from inside
        r = Ray([0.0, 0.0, -0.5], [0.0, -1.0, 0.0])
        intvl = surfaceintersection(az1, r)
        @test (lower(intvl) isa RayOrigin) && isapprox(point(upper(intvl)), [0.0, -1.5, -0.5], rtol = RTOLERANCE, atol = ATOLERANCE)

        z2 = TestData.zernikesurface2()
        az2 = AcceleratedParametricSurface(z2, 20)

        # two hits
        r = Ray([2.0, -1.0, 0.2], [-1.0, 0.0, 0.0])
        intvl = surfaceintersection(az2, r)
        @test isapprox(point(lower(intvl)), [0.238496383738369, -1.0, 0.2], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(intvl)), [-0.8790518227874484, -1.0, 0.2], rtol = RTOLERANCE, atol = ATOLERANCE)

        # three hits
        r = Ray([-0.7, 1.0, 0.2], [0.0, -1.0, 0.0])
        hit = surfaceintersection(az2, r)
        a = (lower(hit[1]) isa RayOrigin) && isapprox(point(upper(hit[1])), [-0.7, 0.8732489598020176, 0.2], rtol = RTOLERANCE, atol = ATOLERANCE)
        b = isapprox(point(lower(hit[2])), [-0.7, -0.34532502965048606, 0.2], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(hit[2])), [-0.7, -1.1615928003236047, 0.2], rtol = RTOLERANCE, atol = ATOLERANCE)
        @test isa(hit, DisjointUnion) && length(hit) == 2 && a && b

        # four hits
        r = Ray([1.5, 0.1, 0.35], [-1.0, -0.5, 0.0])
        hit = surfaceintersection(az2, r)
        a = isapprox(point(lower(hit[1])), [1.4371467631969683, 0.06857338159848421, 0.35], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(hit[1])), [1.1428338578611368, -0.07858307106943167, 0.35], rtol = RTOLERANCE, atol = ATOLERANCE)
        b = isapprox(point(lower(hit[2])), [0.12916355683491865, -0.5854182215825409, 0.35], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(hit[2])), [-0.7290713614371003, -1.0145356807185504, 0.35], rtol = RTOLERANCE, atol = ATOLERANCE)
        @test isa(hit, DisjointUnion) && length(hit) == 2 && a && b

        # failure case, had a bug where rounding error would cause o2 to fail
        o1 = leaf(AcceleratedParametricSurface(ZernikeSurface(1.4 * 1.15)), translation(0.0, 0.0, 3.0))()
        o2 = leaf(AcceleratedParametricSurface(ZernikeSurface(1.4 * 1.15)), translation(2.0, 0.0, 3.0))()
        r = Ray([-5.0, 0.0, 0.0], [1.0, 0.0, 0.0])
        @test !(surfaceintersection(o1, r) isa EmptyInterval)
        @test !(surfaceintersection(o2, r) isa EmptyInterval)
    end # testset Zernike

    @testset "QType" begin
        @test_nowarn QTypeSurface(1.5)
        @test_throws AssertionError QTypeSurface(0.0)

        q1 = TestData.qtypesurface1()
        aq1 = AcceleratedParametricSurface(q1, 20)
        numsamples = 100

        # test that an on axis ray doesn't miss
        @test !(surfaceintersection(AcceleratedParametricSurface(q1), Ray([0.0, 0.0, 10.0], [0.0, 0.0, -1.0])) isa EmptyInterval)

        samples = samplepoints(numsamples, 0.01, 0.99, 0.01, 2π - 0.01)
        missedintersections = 0

        for uv in samples
            pt = point(q1, uv[1], uv[2])
            origin = [0.0, 0.0, 5.0]
            r = Ray(origin, pt .- origin)
            allintersections = surfaceintersection(aq1, r)
            if (allintersections isa EmptyInterval)
                # sampled triangle acceleration structure doesn't guarantee all intersecting rays will be detected.
                # Have to use patch subdivision and convex hull polyhedron to do this.
                missedintersections += 1
            else
                if isa(allintersections, Interval)
                    allintersections = (allintersections,)
                end

                for intsct in allintersections
                    if lower(intsct) isa RayOrigin
                        missedintersections += 1
                    else
                        # closest point should be on the surface, furthest should be on bounding prism
                        @test isapprox(point(lower(intsct)), pt, rtol = RTOLERANCE, atol = ATOLERANCE)
                        @test isa(upper(intsct), Intersection{Float64,3}) || (OpticSim.direction(r)[3] == -1 && isa(upper(intsct), Infinity{Float64}))
                    end
                end
            end
        end

        if missedintersections > 0
            @warn "$missedintersections out of total $(length(samples)) qtype suface intersections were missed"
        end

        # hit from inside
        r = Ray([0.0, 0.0, -0.5], [0.1, 0.2, 0.5])
        intvl = surfaceintersection(aq1, r)
        @test (lower(intvl) isa RayOrigin) && isapprox(point(upper(intvl)), [0.09758258750208074, 0.19516517500416142, -0.01208706248959643], rtol = RTOLERANCE, atol = ATOLERANCE)

        # hit from outside
        r = Ray([0.0, 0.0, 0.5], [0.1, 0.2, -0.5])
        intvl = surfaceintersection(aq1, r)
        @test isapprox(point(lower(intvl)), [0.10265857811957124, 0.20531715623914257, -0.013292890597856405], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(intvl)), [0.6708203932499369, 1.3416407864998738, -2.8541019662496843], rtol = RTOLERANCE, atol = ATOLERANCE)

        # miss from inside
        r = Ray([0.2, 0.2, -0.5], [0.0, 0.0, -1.0])
        intvl = surfaceintersection(aq1, r)
        @test (lower(intvl) isa RayOrigin) && (upper(intvl) isa Infinity)

        # miss from outside
        r = Ray([0.2, 2.0, -0.5], [0.0, 0.0, 1.0])
        @test surfaceintersection(aq1, r) isa EmptyInterval
        r = Ray([0.2, 2.0, 0.5], [0.0, 0.0, -1.0])
        @test surfaceintersection(aq1, r) isa EmptyInterval
        r = Ray([0.2, 0.2, 0.5], [0.0, 0.0, 1.0])
        @test surfaceintersection(aq1, r) isa EmptyInterval
        r = Ray([0.2, 2.0, 2.0], [0.0, -1.0, 0.0])
        @test surfaceintersection(aq1, r) isa EmptyInterval
    end # testset QType

    @testset "BoundingBox" begin
        function boundingval(a::BoundingBox{T}, axis::Int, plane::Bool) where {T<:Real}
            if axis === 1
                return plane ? a.xmax : a.xmin
            elseif axis === 2
                return plane ? a.ymax : a.ymin
            elseif axis == 3
                return plane ? a.zmax : a.zmin
            else
                throw(ErrorException("Invalid axis: $axis"))
            end
        end

        Random.seed!(SEED)

        axis(x) = (x + 1) ÷ 2
        face(x) = mod(x, 2)
        facevalue(a::BoundingBox, x) = boundingval(a, axis(x), Bool(face(x)))
        boundingindices = ((2, 3), (1, 3), (1, 2))

        function facepoint(a::BoundingBox, facenumber)
            axisindex = axis(facenumber)
            indices = boundingindices[axisindex]
            b1min, b1max = boundingval(a, indices[1], false), boundingval(a, indices[1], true)
            b2min, b2max = boundingval(a, indices[2], false), boundingval(a, indices[2], true)

            step = 0.00001
            pt1val = rand((b1min + step):step:(b1max - step))
            pt2val = rand((b2min + step):step:(b2max - step))

            result = Array{Float64,1}(undef, 3)
            result[indices[1]] = pt1val
            result[indices[2]] = pt2val
            result[axisindex] = facevalue(a, facenumber)
            return result
        end

        bbox = BoundingBox(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0)

        for i in 1:10000
            # randomly pick two planes
            face1 = rand(1:6)
            face2 = rand(1:6)
            while face2 == face1
                face2 = rand(1:6)
            end
            # pick a point on each face
            pt1 = facepoint(bbox, face1)
            pt2 = facepoint(bbox, face2)

            direction = normalize(pt1 - pt2)
            origin = pt2 - 3 * direction

            r = Ray(origin, direction)

            @test doesintersect(bbox, r)
        end

        # should miss
        @test !doesintersect(bbox, Ray([-2.0, -2.0, 0.0], [1.0, 0.0, 0.0]))
        @test !doesintersect(bbox, Ray([-2.0, -2.0, 0.0], [0.0, 1.0, 0.0]))
        @test !doesintersect(bbox, Ray([-2.0, -2.0, 0.0], [0.0, 0.0, 1.0]))
        @test !doesintersect(bbox, Ray([-2.0, 0.0, 0.0], [-1.0, 0.0, 0.0]))
        @test !doesintersect(bbox, Ray([0.0, -2.0, 0.0], [0.0, -1.0, 0.0]))
        @test !doesintersect(bbox, Ray([0.0, 0.0, -2.0], [0.0, 0.0, -1.0]))

        bbox = BoundingBox(-0.5, 0.5, -0.75, 0.75, typemin(Float64), typemax(Float64))

        # starting outside
        r = Ray([-1.0, -1.1, 0.0], [1.0, 1.0, 0.0])
        intscts = surfaceintersection(bbox, r)
        @test isapprox(point(lower(intscts)), [-0.5, -0.6, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(intscts)), [0.5, 0.4, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE)

        # starting inside
        r = Ray([0.0, 0.0, 0.0], [1.0, 1.0, 0.0])
        intscts = surfaceintersection(bbox, r)
        @test lower(intscts) isa RayOrigin && isapprox(point(upper(intscts)), [0.5, 0.5, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE)

        # inside infinite
        r = Ray([0.0, 0.0, 0.0], [0.0, 0.0, 1.0])
        intscts = surfaceintersection(bbox, r)
        @test lower(intscts) isa RayOrigin && upper(intscts) isa Infinity

        # on surface
        r = Ray([0.5, 0.0, 0.0], [0.0, 0.0, 1.0])
        intscts = surfaceintersection(bbox, r)
        @test lower(intscts) isa RayOrigin && upper(intscts) isa Infinity

        # missing
        r = Ray([5.0, 5.0, 5.0], [0.0, 0.0, -1.0])
        intscts = surfaceintersection(bbox, r)
        @test intscts isa EmptyInterval

        r = Ray([5.0, 5.0, 5.0], [0.0, -1.0, 0.0])
        intscts = surfaceintersection(bbox, r)
        @test intscts isa EmptyInterval
    end # testset BoundingBox

    @testset "CSG" begin
        pln1 = Plane([0.0, 0.0, -1.0], [0.0, 0.0, -1.0])
        pln2 = Plane([0.0, 0.0, 1.0], [0.0, 0.0, 1.0])
        r = Ray([0.0, 0, 2.0], [0.0, 0.0, -1.0])
        gen1 = csgintersection(pln1, pln2)
        gen2 = csgunion(pln1, pln2)

        csg = gen1(identitytransform())
        intsct = evalcsg(csg, r)

        intvlpt1 = lower(intsct)
        intvlpt2 = upper(intsct)

        @test isapprox(point(intvlpt1), SVector{3}(0.0, 0.0, 1.0)) && isapprox(point(intvlpt2), SVector{3}(0.0, 0.0, -1.0))

        csg = gen2(identitytransform())
        intsct = evalcsg(csg, r)
        @test (lower(intsct) isa RayOrigin) && (upper(intsct) isa Infinity)

        # test simple rays for intersection, union and difference operations
        # INTERSECTION
        intersection_obj = csgintersection(leaf(Cylinder(0.5, 3.0), OpticSim.rotationd(90.0, 0.0, 0.0)), (leaf(Sphere(1.0))))()
        r = Ray([0.7, 0.0, 0.0], [-1.0, 0.0, 0.0])
        int = OpticSim.evalcsg(intersection_obj, r)
        @test isapprox(point(lower(int)), [0.5, 0.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(int)), [-0.5, 0.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE)

        r = Ray([5.0, 0.0, 0.0], [-1.0, 0.0, 0.0])
        int = OpticSim.evalcsg(intersection_obj, r)
        @test isapprox(point(lower(int)), [0.5, 0.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(int)), [-0.5, 0.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE)

        r = Ray([0.7, 0.0, 0.0], [1.0, 0.0, 0.0])
        int = OpticSim.evalcsg(intersection_obj, r)
        @test (int isa EmptyInterval)

        r = Ray([0.0, 0.0, 0.0], [1.0, 0.0, 0.0])
        int = OpticSim.evalcsg(intersection_obj, r)
        @test (lower(int) isa RayOrigin) && isapprox(point(upper(int)), [0.5, 0.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE)

        r = Ray([0.0, 0.0, 0.0], [0.0, 1.0, 0.0])
        int = OpticSim.evalcsg(intersection_obj, r)
        @test (lower(int) isa RayOrigin) && isapprox(point(upper(int)), [0.0, 1.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE)

        r = Ray([0.0, 2.0, 0.0], [0.0, -1.0, 0.0])
        int = OpticSim.evalcsg(intersection_obj, r)
        @test isapprox(point(lower(int)), [0.0, 1.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(int)), [0.0, -1.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE)

        # UNION
        union_obj = csgunion(leaf(Cylinder(0.5, 3.0)), (leaf(Sphere(1.0))))(OpticSim.translation(1.0, 0.0, 0.0))
        r = Ray([1.0, 0.0, 0.0], [0.0, 0.0, 1.0])
        int = OpticSim.evalcsg(union_obj, r)
        @test (lower(int) isa RayOrigin) && (upper(int) isa Infinity)

        r = Ray([1.0, 0.0, 0.0], [1.0, 0.0, 0.0])
        int = OpticSim.evalcsg(union_obj, r)
        @test (lower(int) isa RayOrigin) && isapprox(point(upper(int)), [2.0, 0.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE)

        r = Ray([1.6, 0.0, 0.0], [1.0, 0.0, 0.0])
        int = OpticSim.evalcsg(union_obj, r)
        @test (lower(int) isa RayOrigin) && isapprox(point(upper(int)), [2.0, 0.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE)

        r = Ray([1.6, 0.0, 0.0], [-1.0, 0.0, 0.0])
        int = OpticSim.evalcsg(union_obj, r)
        @test (lower(int) isa RayOrigin) && isapprox(point(upper(int)), [0.0, 0.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE)

        r = Ray([-1.0, 0.0, 0.0], [1.0, 0.0, 0.0])
        int = OpticSim.evalcsg(union_obj, r)
        @test isapprox(point(lower(int)), [0.0, 0.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(int)), [2.0, 0.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE)

        r = Ray([-1.0, 0.0, 4.0], [1.0, 0.0, 0.0])
        int = OpticSim.evalcsg(union_obj, r)
        @test isapprox(point(lower(int)), [0.5, 0.0, 4.0], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(int)), [1.5, 0.0, 4.0], rtol = RTOLERANCE, atol = ATOLERANCE)

        # DIFFERENCE
        difference_obj = csgdifference(leaf(Cylinder(0.5, 3.0)), leaf(Sphere(1.0), OpticSim.translation(0.75, 0.0, 0.2)))()
        r = Ray([0.25, 0.0, 0.0], [-1.0, 0.0, 0.0])
        int = OpticSim.evalcsg(difference_obj, r)
        @test isapprox(point(lower(int)), [-0.2297958971132712, 0.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(int)), [-0.5, 0.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE)

        r = Ray([0.25, 0.0, 0.0], [1.0, 0.0, 0.0])
        int = OpticSim.evalcsg(difference_obj, r)
        @test int isa EmptyInterval

        r = Ray([0.25, 0.0, 0.0], [0.0, 0.0, 1.0])
        int = OpticSim.evalcsg(difference_obj, r)
        @test isapprox(point(lower(int)), [0.25, 0.0, 1.0660254037844386], rtol = RTOLERANCE, atol = ATOLERANCE) && (upper(int) isa Infinity)

        r = Ray([0.25, 0.0, 1.5], [0.0, 0.0, -1.0])
        int = OpticSim.evalcsg(difference_obj, r)
        @test (lower(int[1]) isa RayOrigin) && isapprox(point(upper(int[1])), [0.25, 0.0, 1.0660254037844386], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(lower(int[2])), [0.25, 0.0, -0.6660254037844386], rtol = RTOLERANCE, atol = ATOLERANCE) && (upper(int[2]) isa Infinity)

        r = Ray([0.25, 0.0, 1.5], [1.0, 0.0, 0.0])
        int = OpticSim.evalcsg(difference_obj, r)
        @test (lower(int) isa RayOrigin) && isapprox(point(upper(int)), [0.5, 0.0, 1.5], rtol = RTOLERANCE, atol = ATOLERANCE)

        r = Ray([0.25, 0.0, 1.5], [0.0, 0.0, 1.0])
        int = OpticSim.evalcsg(difference_obj, r)
        @test (lower(int) isa RayOrigin) && (upper(int) isa Infinity)

        # DisjointUnion result on CSG
        surf = TestData.verywavybeziersurface()
        accelsurf = leaf(AcceleratedParametricSurface(surf, 20))()

        # five hits starting outside
        r = Ray([0.9, 0.0, -0.3], [0.0, 1.0, 0.7])
        res = surfaceintersection(accelsurf, r)
        a = isapprox(point(lower(res[1])), [0.9, 0.03172286522032046, -0.2777939943457758], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(res[1])), [0.9, 0.1733979947040411, -0.17862140370717122], rtol = RTOLERANCE, atol = ATOLERANCE)
        b = isapprox(point(lower(res[2])), [0.9, 0.5335760974594397, 0.07350326822160776], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(res[2])), [0.9, 0.7767707607392784, 0.24373953251749486], rtol = RTOLERANCE, atol = ATOLERANCE)
        c = isapprox(point(lower(res[3])), [0.9, 0.9830891958374246, 0.3881624370861975], rtol = RTOLERANCE, atol = ATOLERANCE) && (upper(res[3]) isa Infinity)
        @test isa(res, DisjointUnion) && length(res) == 3 && a && b && c

        # five hits starting inside
        r = Ray([0.9, 1.0, 0.4], [0.0, -1.0, -0.7])
        res = surfaceintersection(accelsurf, r)
        a = (lower(res[1]) isa RayOrigin) && isapprox(point(upper(res[1])), [0.9, 0.9830891958374246, 0.3881624370861975], rtol = RTOLERANCE, atol = ATOLERANCE)
        b = isapprox(point(lower(res[2])), [0.9, 0.7767707607392784, 0.24373953251749486], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(res[2])), [0.9, 0.5335760974594397, 0.07350326822160776], rtol = RTOLERANCE, atol = ATOLERANCE)
        c = isapprox(point(lower(res[3])), [0.9, 0.17339799470404108, -0.1786214037071712], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(res[3])), [0.9, 0.03172286522032046, -0.27779399434577573], rtol = RTOLERANCE, atol = ATOLERANCE)
        @test isa(res, DisjointUnion) && length(res) == 3 && a && b && c
    end # testset CSG

    @testset "ThinGrating" begin
        angle_from_ray(raydirection) = 90 + atand(raydirection[3], raydirection[2])
        true_diff(order, λ, period, θi) = asind((order * λ / period + sind(θi)))
        period = 3.0
        int = TestData.transmissivethingrating(period, 2)
        for k in 1:50
            for θi in [-5.0, 0.0, 5.0, 10.0]
                for λ in [0.35, 0.55, 1.0]
                    ray = OpticalRay(SVector(0.0, 0.0, 2.0), SVector(0.0, sind(θi), -cosd(θi)), 1.0, λ)
                    raydir, _, _ = OpticSim.processintersection(int, SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 1.0), ray, 20.0, 1.0, true, true)
                    # order is random so just check that it is correct for one of them
                    @test any([isapprox(x, angle_from_ray(raydir), rtol = RTOLERANCE, atol = ATOLERANCE) for x in [true_diff(m, λ, period, θi) for m in -2:2]])
                end
            end
        end
    end # testset ThinGrating

    # TODO Hologram intersection tests

    @testset "Chebyshev" begin
        @test_throws AssertionError ChebyshevSurface(0.0, 1.0, [(1, 2, 1.0)])
        @test_throws AssertionError ChebyshevSurface(1.0, 0.0, [(1, 2, 1.0)])

        z1 = TestData.chebyshevsurface2()
        az1 = AcceleratedParametricSurface(z1, 20)
        numsamples = 100

        # test that an on axis ray doesn't miss
        @test !(surfaceintersection(AcceleratedParametricSurface(z1), Ray([0.0, 0.0, 10.0], [0.0, 0.0, -1.0])) isa EmptyInterval)

        samples = samplepoints(numsamples, -0.95, 0.95, -0.95, 0.95)
        missedintersections = 0

        for uv in samples
            pt = point(z1, uv[1], uv[2])
            origin = [0.0, 0.0, 5.0]
            r = Ray(origin, pt .- origin)
            allintersections = surfaceintersection(az1, r)
            if (allintersections isa EmptyInterval)
                # sampled triangle acceleration structure doesn't guarantee all intersecting rays will be detected.
                # Have to use patch subdivision and convex hull polyhedron to do this.
                missedintersections += 1
            else
                if isa(allintersections, Interval)
                    allintersections = (allintersections,)
                end

                for intsct in allintersections
                    if lower(intsct) isa RayOrigin
                        missedintersections += 1
                    else
                        # closest point should be on the surface, furthest should be on bounding prism
                        @test isapprox(point(lower(intsct)), pt, rtol = RTOLERANCE, atol = ATOLERANCE)
                        @test isa(upper(intsct), Intersection{Float64,3}) || (OpticSim.direction(r)[3] == -1 && isa(upper(intsct), Infinity{Float64}))
                    end
                end
            end
        end

        if missedintersections > 0
            @warn "$missedintersections out of total $(length(samples)) chebychev suface intersections were missed"
        end

        z1 = TestData.chebyshevsurface1()
        az1 = AcceleratedParametricSurface(z1, 20)

        # hit from inside
        r = Ray([0.0, 0.0, -0.5], [0.1, 0.2, 0.5])
        intvl = surfaceintersection(az1, r)
        @test (lower(intvl) isa RayOrigin) && isapprox(point(upper(intvl)), [0.07870079995991423, 0.15740159991982847, -0.10649600020042879], rtol = RTOLERANCE, atol = ATOLERANCE)

        # hit from outside
        r = Ray([0.0, 0.0, 0.5], [0.1, 0.2, -0.5])
        intvl = surfaceintersection(az1, r)
        @test isapprox(point(lower(intvl)), [0.13584702907541094, 0.2716940581508219, -0.1792351453770547], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(intvl)), [1.0, 2.0, -4.5], rtol = RTOLERANCE, atol = ATOLERANCE)

        # miss from inside
        r = Ray([0.2, 0.2, -0.5], [0.0, 0.0, -1.0])
        intvl = surfaceintersection(az1, r)
        @test (lower(intvl) isa RayOrigin) && (upper(intvl) isa Infinity)

        # miss from outside
        r = Ray([0.2, 3.0, -0.5], [0.0, 0.0, 1.0])
        @test surfaceintersection(az1, r) isa EmptyInterval
        r = Ray([0.2, 3.0, 0.5], [0.0, 0.0, -1.0])
        @test surfaceintersection(az1, r) isa EmptyInterval
        r = Ray([0.2, 0.2, 0.5], [0.0, 0.0, 1.0])
        @test surfaceintersection(az1, r) isa EmptyInterval
        r = Ray([0.2, 3.0, 2.0], [0.0, -1.0, 0.0])
        @test surfaceintersection(az1, r) isa EmptyInterval

        # two hits from outside
        r = Ray([0.0, -1.5, 0.25], [-1.0, 0.0, 0.0])
        intvl = surfaceintersection(az1, r)
        @test isapprox(point(lower(intvl)), [-1.030882920068228, -1.5, 0.25], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(intvl)), [-1.786071230727613, -1.5, 0.25], rtol = RTOLERANCE, atol = ATOLERANCE)

        r = Ray([1.5, -0.8, 0.25], [-1.0, 0.0, 0.0])
        hit = surfaceintersection(az1, r)
        a = isapprox(point(lower(hit[1])), [0.21526832427065165, -0.8, 0.25], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(hit[1])), [-0.25523748548950076, -0.8, 0.25], rtol = RTOLERANCE, atol = ATOLERANCE)
        b = isapprox(point(lower(hit[2])), [-1.9273057166986296, -0.8, 0.25], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(hit[2])), [-2.0, -0.8, 0.25], rtol = RTOLERANCE, atol = ATOLERANCE)
        @test isa(hit, DisjointUnion) && length(hit) == 2 && a && b

        # hit cyl from outside
        r = Ray([0.0, 3.0, -0.5], [0.0, -1.0, 0.0])
        intvl = surfaceintersection(az1, r)
        @test isapprox(point(lower(intvl)), [0.0, 2.0, -0.5], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(intvl)), [0.0, -2.0, -0.5], rtol = RTOLERANCE, atol = ATOLERANCE)

        # hit cyl from inside
        r = Ray([0.0, 0.0, -0.5], [0.0, -1.0, 0.0])
        intvl = surfaceintersection(az1, r)
        @test (lower(intvl) isa RayOrigin) && isapprox(point(upper(intvl)), [0.0, -2.0, -0.5], rtol = RTOLERANCE, atol = ATOLERANCE)

    end # testset Chebyshev

    @testset "GridSag" begin
        z1 = TestData.gridsagsurfacebicubic()
        az1 = AcceleratedParametricSurface(z1, 20)

        numsamples = 100

        # test that an on axis ray doesn't miss
        @test !(surfaceintersection(AcceleratedParametricSurface(z1), Ray([0.0, 0.0, 10.0], [0.0, 0.0, -1.0])) isa EmptyInterval)

        samples1 = samplepoints(numsamples, 0.01, 0.99, 0.01, 2π - 0.01)
        missedintersections = 0

        for uv in samples1
            pt = point(z1, uv[1], uv[2])
            origin = [0.0, 0.0, 5.0]
            r = Ray(origin, pt .- origin)
            allintersections = surfaceintersection(az1, r)
            if (allintersections isa EmptyInterval)
                # sampled triangle acceleration structure doesn't guarantee all intersecting rays will be detected.
                # Have to use patch subdivision and convex hull polyhedron to do this.
                missedintersections += 1
            else
                if isa(allintersections, Interval)
                    allintersections = (allintersections,)
                end

                for intsct in allintersections
                    if lower(intsct) isa RayOrigin
                        missedintersections += 1
                    else
                        # closest point should be on the surface, furthest should be on bounding prism
                        @test isapprox(point(lower(intsct)), pt, rtol = RTOLERANCE, atol = ATOLERANCE)
                        @test isa(upper(intsct), Intersection{Float64,3}) || (OpticSim.direction(r)[3] == -1 && isa(upper(intsct), Infinity{Float64}))
                    end
                end
            end
        end

        z2 = TestData.gridsagsurfacebicubiccheby()
        az2 = AcceleratedParametricSurface(z2, 20)

        samples2 = samplepoints(numsamples, -0.95, 0.95, -0.95, 0.95)
        for uv in samples2
            pt = point(z2, uv[1], uv[2])
            origin = [0.0, 0.0, 5.0]
            r = Ray(origin, pt .- origin)
            allintersections = surfaceintersection(az2, r)
            if (allintersections isa EmptyInterval)
                # sampled triangle acceleration structure doesn't guarantee all intersecting rays will be detected.
                # Have to use patch subdivision and convex hull polyhedron to do this.
                missedintersections += 1
            else
                if isa(allintersections, Interval)
                    allintersections = (allintersections,)
                end

                for intsct in allintersections
                    if lower(intsct) isa RayOrigin
                        missedintersections += 1
                    else
                        # closest point should be on the surface, furthest should be on bounding prism
                        @test isapprox(point(lower(intsct)), pt, rtol = RTOLERANCE, atol = ATOLERANCE)
                        @test isa(upper(intsct), Intersection{Float64,3}) || (OpticSim.direction(r)[3] == -1 && isa(upper(intsct), Infinity{Float64}))
                    end
                end
            end
        end

        if missedintersections > 0
            @warn "$missedintersections out of total $(length(samples1) + length(samples2)) gridsag suface intersections were missed"
        end

        # hit from inside
        r = Ray([0.0, 0.0, -0.5], [0.1, 0.2, 0.5])
        intvl = surfaceintersection(az1, r)
        @test (lower(intvl) isa RayOrigin) && isapprox(point(upper(intvl)), [0.13038780645120746, 0.26077561290241486, 0.15193903225603717], rtol = RTOLERANCE, atol = ATOLERANCE)

        # hit from outside
        r = Ray([0.0, 0.0, 0.5], [0.1, 0.2, -0.5])
        intvl = surfaceintersection(az1, r)
        @test isapprox(point(lower(intvl)), [0.046677740288643785, 0.0933554805772876, 0.26661129855678095], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(intvl)), [0.6708203932499369, 1.3416407864998738, -2.8541019662496843], rtol = RTOLERANCE, atol = ATOLERANCE)

        # miss from inside
        r = Ray([0.2, 0.2, -0.5], [0.0, 0.0, -1.0])
        intvl = surfaceintersection(az1, r)
        @test (lower(intvl) isa RayOrigin) && (upper(intvl) isa Infinity)

        # miss from outside
        r = Ray([0.2, 2.0, -0.5], [0.0, 0.0, 1.0])
        @test surfaceintersection(az1, r) isa EmptyInterval
        r = Ray([0.2, 2.0, 0.5], [0.0, 0.0, -1.0])
        @test surfaceintersection(az1, r) isa EmptyInterval
        r = Ray([0.2, 0.2, 0.5], [0.0, 0.0, 1.0])
        @test surfaceintersection(az1, r) isa EmptyInterval
        r = Ray([0.2, 2.0, 2.0], [0.0, -1.0, 0.0])
        @test surfaceintersection(az1, r) isa EmptyInterval

        # two hits from outside
        r = Ray([2.0, 0.0, 0.25], [-1.0, 0.0, 0.0])
        hit = surfaceintersection(az1, r)
        a = isapprox(point(lower(hit[1])), [1.5, 0.0, 0.25], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(hit[1])), [1.394132887629247, 0.0, 0.25], rtol = RTOLERANCE, atol = ATOLERANCE)
        b = isapprox(point(lower(hit[2])), [0.30184444700168467, 0.0, 0.25], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(hit[2])), [-0.6437764440680096, 0.0, 0.25], rtol = RTOLERANCE, atol = ATOLERANCE)
        @test isa(hit, DisjointUnion) && length(hit) == 2 && a && b

        # hit cyl from outside
        r = Ray([0.0, 2.0, -0.5], [0.0, -1.0, 0.0])
        intvl = surfaceintersection(az1, r)
        @test isapprox(point(lower(intvl)), [0.0, 1.5, -0.5], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(intvl)), [0.0, -1.5, -0.5], rtol = RTOLERANCE, atol = ATOLERANCE)

        # hit cyl from inside
        r = Ray([0.0, 0.0, -0.5], [0.0, -1.0, 0.0])
        intvl = surfaceintersection(az1, r)
        @test (lower(intvl) isa RayOrigin) && isapprox(point(upper(intvl)), [0.0, -1.5, -0.5], rtol = RTOLERANCE, atol = ATOLERANCE)
    end # testset GridSag
end # testset intersection
