# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

@testset "Allocations" begin
    # ensure that there are 0 allocations for all the benchmarks
    for b in Benchmarks.all_benchmarks()
        @test Benchmarks.runbenchmark(b, samples = 1, evals = 1).allocs == 0
    end
end # testset Allocations
