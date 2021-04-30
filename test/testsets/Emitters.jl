# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

using Test
using OpticSim
using OpticSim.Emitters
using OpticSim.Geometry

using Random
using StaticArrays

@testset "Emitters" begin
    @testset "AngularPower" begin
        @testset "Lambertian" begin
            @test typeof(AngularPower.Lambertian()) === AngularPower.Lambertian{Float64}
            @test_throws MethodError AngularPower.Lambertian(String)

            @test apply(AngularPower.Lambertian(), Transform(), 1, Ray(zeros(3), ones(3))) === 1
        end

        @testset "Cosine" begin
            @test AngularPower.Cosine().cosine_exp === one(Float64)
            @test_throws MethodError AngularPower.Cosine(String)

            @test apply(AngularPower.Cosine(), Transform(), 1, Ray(zeros(3), ones(3))) === 0.5773502691896258
        end

        @testset "Gaussian" begin
            @test AngularPower.Gaussian(0, 1).gaussianu === 0
            @test AngularPower.Gaussian(0, 1).gaussianv === 1

            @test apply(AngularPower.Gaussian(0, 1), Transform(), 1, Ray(zeros(3), ones(3))) === 0.7165313105737892
        end
    end

    @testset "Directions" begin
        @testset "Constant" begin
            @test Directions.Constant(0, 0, 0).direction === Vec3(Int)
            @test Directions.Constant(Vec3()).direction === Vec3()
            @test Directions.Constant().direction === unitZ3()

            @test Base.length(Directions.Constant()) === 1
            @test Emitters.generate(Directions.Constant(), 0) === unitZ3()

            @test Base.iterate(Directions.Constant()) === (unitZ3(), 2)
            @test Base.iterate(Directions.Constant(), 2) |> isnothing
            @test Base.getindex(Directions.Constant(), 1) === unitZ3()
            @test Base.getindex(Directions.Constant(), 2) === unitZ3()
            @test Base.firstindex(Directions.Constant()) === 0
            @test Base.lastindex(Directions.Constant()) === 0
            @test Base.copy(Directions.Constant()) === Directions.Constant()
        end

        @testset "RectGrid" begin
            @test Directions.RectGrid(unitX3(), 0., 0., 0, 0).uvec === unitY3()
            @test Directions.RectGrid(unitX3(), 0., 0., 0, 0).vvec === unitZ3()
            @test Directions.RectGrid(0., 0., 0, 0).uvec === unitX3()
            @test Directions.RectGrid(0., 0., 0, 0).vvec === unitY3()

            @test Base.length(Directions.RectGrid(0., 0., 2, 3)) === 6
            @test collect(Directions.RectGrid(ones(Vec3), π/4, π/4, 2, 3)) == [
                [0.6014252485829169, 0.7571812496798511, 1.2608580392339483],
                [1.1448844430815608, 0.4854516524305293, 0.9891284419846265],
                [0.643629619770125, 1.0798265148205617, 1.0798265148205617],
                [1.225225479837374, 0.7890285847869373, 0.7890285847869373],
                [0.6014252485829169, 1.2608580392339483, 0.7571812496798511],
                [1.1448844430815608, 0.9891284419846265, 0.4854516524305293],
            ]
        end

        @testset "UniformCone" begin
            @test Directions.UniformCone(unitX3(), 0., 0).uvec === unitY3()
            @test Directions.UniformCone(unitX3(), 0., 0).vvec === unitZ3()
            @test Directions.UniformCone(0., 0).uvec === unitX3()
            @test Directions.UniformCone(0., 0).vvec === unitY3()

            @test Base.length(Directions.UniformCone(0., 1)) === 1
            Random.seed!(0)
            @test collect(Directions.UniformCone(π/4, 2)) == [
                [0.30348115383395624, -0.6083405920618145, 0.7333627433388552],
                [0.16266571675478964, 0.2733479444418462, 0.9480615833700192],
            ]
        end

        @testset "HexapolarCone" begin
            @test Directions.HexapolarCone(unitX3(), 0., 0).uvec === unitY3()
            @test Directions.HexapolarCone(unitX3(), 0., 0).vvec === unitZ3()
            @test Directions.HexapolarCone(0., 0).uvec === unitX3()
            @test Directions.HexapolarCone(0., 0).vvec === unitY3()

            @test Base.length(Directions.HexapolarCone(0., 1)) === 7
            @test Base.length(Directions.HexapolarCone(0., 2)) === 19
            @test collect(Directions.HexapolarCone(π/4, 1)) == [
                [0.0, 0.0, 1.0],
                [0.7071067811865475, 0.0, 0.7071067811865476],
                [0.3535533905932738, 0.6123724356957945, 0.7071067811865476],
                [-0.3535533905932736, 0.6123724356957946, 0.7071067811865477],
                [-0.7071067811865475, 8.659560562354932e-17, 0.7071067811865476],
                [-0.35355339059327406, -0.6123724356957944, 0.7071067811865476],
                [0.3535533905932738, -0.6123724356957945, 0.7071067811865476],
            ]
        end
    end

    @testset "Origins" begin
        
    end

    @testset "Sources" begin
        
    end

    @testset "Spectrum" begin
        
    end
end # testset Emitters
