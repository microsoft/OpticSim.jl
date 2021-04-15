@testset "Visualization" begin
    if get(ENV, "CI", nothing) == "true" && Sys.iswindows()
        # OpenGL is unreliable on headless Windows VMs (Github Action servers). Skip.
        return
    end

    @testcase begin
        # test that this all at least runs
        surf1 = AcceleratedParametricSurface(TestData.beziersurface(), 15)
        surf2 = AcceleratedParametricSurface(TestData.upsidedownbeziersurface(), 15)

        shapes = [
            csgintersection(leaf(surf1, translation(-0.5, -0.5, 0.0)), csgintersection(leaf(Cylinder(0.3, 5.0)), leaf(surf1, Transform{Float64}(0.0, Float64(π), 0.0, 0.5, -0.5, 0.0))))(),
            csgintersection(leaf(surf1, translation(-0.5, -0.5, 0.0)), csgintersection(leaf(Cylinder(0.3, 5.0)), leaf(surf2, translation(-0.5, -0.5, 0.0))))(),
            csgunion(leaf(Cylinder(0.5, 3.0)), (leaf(Sphere(1.0))))(),
            csgintersection(leaf(Cylinder(0.5, 3.0)), (leaf(Sphere(1.0))))(),
            csgdifference(leaf(Cylinder(0.5, 3.0)), leaf(Sphere(1.0)))()
        ]

        for shape in shapes
            @test_nowarn Vis.draw(shape)
        end
    end
end # testset Visualization
