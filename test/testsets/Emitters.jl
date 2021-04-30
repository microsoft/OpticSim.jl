# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

using Test
using OpticSim
using OpticSim.Emitters
using OpticSim.Geometry

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
        
    end

    @testset "Origins" begin
        
    end

    @testset "Sources" begin
        
    end

    @testset "Spectrum" begin
        
    end
end # testset Emitters
