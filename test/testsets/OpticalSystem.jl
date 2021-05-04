# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

@testset "OpticalSystem" begin
    @testset "Single threaded trace makes sure function executes properly" begin
        conv = Examples.doubleconvex()
        rays = Emitters.Sources.CompositeSource(translation(-unitZ3()), repeat([Emitters.Sources.Source()], 100))
        trace(conv, rays)
        @test true #just want to verify that the trace function executed properly
    end
end
