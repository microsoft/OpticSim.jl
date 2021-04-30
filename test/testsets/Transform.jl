# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

@testset "Transform" begin
    A = [
        1 2 3 4;
        4 5 6 7;
        8 9 10 11;
        12 13 14 15
    ]
    x,y,z,w = Vec4.([A[:,i] for i in 1:4])
    @test Geometry.Transform(x,y,z,w) == A

    B = [
        1 2 3 4;
        4 5 6 7;
        8 9 10 11;
        0 0 0 1
    ]
    x,y,z,w = Vec3.([B[1:3,i] for i in 1:4])
    @test B == Geometry.Transform(x,y,z,w)

end # testset Allocations
