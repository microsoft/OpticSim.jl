@testset "Polarization" begin
    function chipman(k̂ᵣ, n̂, nᵢ, nₜ) 
        (sinθᵢ, sinθₜ) = OpticSim.snell(n̂, k̂ᵣ, nᵢ, nₜ)
        rₛ,tₛ,rₚ,tₚ,Tₐ = OpticSim.fresnel(nᵢ, nₜ, sinθᵢ, sinθₜ)
        reflectedray = OpticSim.reflectedray(n̂,k̂ᵣ)
        refractedray = OpticSim.refractedray(nᵢ,nₜ, n̂,k̂ᵣ)
        P = OpticSim.Polarization.Chipman{Float64}(SVector(Complex.((0.0,0.0,1.0))...),SMatrix{3,3,Complex{Float64},9}(I))
        Preflected = OpticSim.Polarization.composepolarization(Complex(rₛ),Complex(rₚ),1.0,n̂,k̂ᵣ,reflectedray,P)
        Prefracted = OpticSim.Polarization.composepolarization(Complex(rₛ),Complex(rₚ),1.0,n̂,k̂ᵣ,refractedray,P)
            
        return Preflected,Prefracted
    end

    function chipmanexample9_4()
        k̂ᵣ = SVector{3,Float64}(0.0, sin(π/6), cos(π/6))
        n̂ = SVector{3,Float64}(0.0, 0.0, 1.0)

        nᵢ = 1.0
        nₜ = 1.5
        return chipman(k̂ᵣ, n̂, nᵢ, nₜ)
    end

    function chipmanexample9_5()
        k̂ = SVector(0.0,sind(10),cosd(10))
        n̂ = SVector(0.0,sind(10),cosd(10))
        jonematrix = SVector{3,3,Complex{Float64},9}(
            
        )
    end

    #test data is from Polarized Light and Optical Systems by Chipman et al pg 333 example 9.4
    @test isapprox(real.(OpticSim.Polarization.pmatrix(chipmanexample9_4()[1])), [
        -.240408 0.0        0.0;
        0       .130825    .501818;
        0       -.501818    -.710275
        ],
        rtol = 1e-5
    )

    #not used yet but will when we write a test that uses aluminum mirrors
    function aluminumfresnel()
        #index of refraction of aluminum at 633nm is 1.374 + 7.620im
        nₜ = 1.374 + 7.620im
        nᵢ = 1.0 #assume air
        s = Vector{Complex{Float64}}(undef,0)
        p = Vector{Complex{Float64}}(undef,0)
    
        for θ in 0.0:.01:π/2
            rₛ,_,rₚ,_,_ = fresnelcomplex(nᵢ,nₜ,sin(θ),0.0)
            push!(s,rₛ)
            push!(p,rₚ)
        end
    
        return s,p
    end

end #testset
