@otestset "Visualization" begin
    using GLMakie
    GLMakie.activate!()
    # test that this all at least runs
    surf1 = AcceleratedParametricSurface(TestData.beziersurface(), 15)
    surf2 = AcceleratedParametricSurface(TestData.upsidedownbeziersurface(), 15)
    m = csgintersection(leaf(surf1, translation(-0.5, -0.5, 0.0)), csgintersection(leaf(Cylinder(0.3, 5.0)), leaf(surf1, RigidBodyTransform{Float64}(0.0, Float64(π), 0.0, 0.5, -0.5, 0.0))))()
    @test_nowarn Vis.draw(m)
    m = csgintersection(leaf(surf1, translation(-0.5, -0.5, 0.0)), csgintersection(leaf(Cylinder(0.3, 5.0)), leaf(surf2, translation(-0.5, -0.5, 0.0))))()
    @test_nowarn Vis.draw(m)
    m = csgunion(leaf(Cylinder(0.5, 3.0)), (leaf(Sphere(1.0))))()
    @test_nowarn Vis.draw(m)
    m = csgintersection(leaf(Cylinder(0.5, 3.0)), (leaf(Sphere(1.0))))()
    @test_nowarn Vis.draw(m)
    m = csgdifference(leaf(Cylinder(0.5, 3.0)), leaf(Sphere(1.0)))()
    @test_nowarn Vis.draw(m)
end # testset Visualization
