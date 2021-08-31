# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

@testset "Visualization" begin
    if get(ENV, "CI", nothing) == "true" && Sys.iswindows() return #OpenGL is unreliable on headless Windows VMs on github. This fails unpredictably and prevents pull requests from being approved. These tests are not essential, so turn them off when running on windows.
    else
        # empty Vis.draw() call to clear Makie @info message:
        # "Info: Makie/Makie is caching fonts, this may take a while. Needed only on first run!"
        Vis.draw()

        # test that this all at least runs
        surf1 = AcceleratedParametricSurface(TestData.beziersurface(), 15)
        surf2 = AcceleratedParametricSurface(TestData.upsidedownbeziersurface(), 15)

        m = (
            leaf(surf1, translation(-0.5, -0.5, 0.0)) ∩
            Cylinder(0.3, 5.0) ∩
            leaf(surf1, Transform{Float64}(0.0, Float64(π), 0.0, 0.5, -0.5, 0.0))
        )()
        Vis.draw(m)
        @test_nowarn Vis.draw(m)

        m = (
            leaf(surf1, translation(-0.5, -0.5, 0.0)) ∩
            Cylinder(0.3, 5.0) ∩
            leaf(surf2, translation(-0.5, -0.5, 0.0))
        )()
        @test_nowarn Vis.draw(m)

        m = (Cylinder(0.5, 3.0) ∪ Sphere(1.0))()
        @test_nowarn Vis.draw(m)

        m = (Cylinder(0.5, 3.0) ∩ Sphere(1.0))()
        @test_nowarn Vis.draw(m)

        m = (Cylinder(0.5, 3.0) - Sphere(1.0))()
        @test_nowarn Vis.draw(m)

        @test_nowarn Vis.draw(Examples.hex3RGB(),[0 1 0;0 0 1])
        @test_nowarn Examples.drawhexrect()
        @test_nowarn Examples.drawhexneighbors()
    end
end # testset Visualization
