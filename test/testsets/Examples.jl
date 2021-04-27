# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

# TODO may want to test this but it's very slow on azure
@testset "Examples" begin
    @test_all_no_arg_functions Examples
end # testset Examples
