@testset "Allocations" begin
    include("../Benchmarks/Benchmarks.jl")
    @testcase "All Benchmarks" begin
        # ensure that there are 0 allocations for all the benchmarks
        for b in Benchmarks.all_benchmarks()
            @test Benchmarks.runbenchmark(b, samples = 1, evals = 1).allocs == 0
        end
    end
end # testset Allocations
