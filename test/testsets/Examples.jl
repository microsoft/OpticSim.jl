# TODO may want to test this but it's very slow on azure
@testset "Examples" begin
    @testcase begin
        @test_all_no_arg_functions Examples
    end
end # testset Examples
