# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

"""
Optimization interface consists of two functions `optimizationvariables` and `updateoptimizationvariables`.
`optimizationvariables` packs variables to be optimized into a vector.
`updateoptimizationvariables` receives a vector of variables and creates a new optical system with the variable values.
"""
module Optimization

using ..OpticSim
using ..OpticSim: detectorsize, temperature, pressure
using ..OpticSim.Emitters

include("Constraints.jl")


"""
    optimizationvariables(a::AxisymmetricOpticalSystem{T}) -> Vector{T}

Pack variables that have been marked to be optimized into a vector in a form suitable for the optimizer.
Variables are marked for optimization by having a true value in the `:OptimizeName` column, where `Name` can be `Radius`, `Thickness` or `Conic`.
"""
function optimizationvariables(a::AxisymmetricOpticalSystem{T}) where {T<:Real}
    prescription = a.prescription

    radiusvariables = in(:OptimizeRadius, propertynames(prescription)) ? prescription[prescription[!, :OptimizeRadius], :Radius] : []
    radiuslower = [typemin(T) for _ in 1:length(radiusvariables)]
    radiusupper = [typemax(T) for _ in 1:length(radiusvariables)]
    thicknessvariables = in(:OptimizeThickness, propertynames(prescription)) ? prescription[prescription[!, :OptimizeThickness], :Thickness] : []
    thicknesslower = [zero(T) for _ in 1:length(thicknessvariables)] # TODO
    thicknessupper = [T(100) for _ in 1:length(thicknessvariables)] # TODO
    conicvariables = in(:OptimizeConic, propertynames(prescription)) ? prescription[prescription[!, :OptimizeConic], :Conic] : []
    coniclower = [typemin(T) for _ in 1:length(radiusvariables)]
    conicupper = [typemax(T) for _ in 1:length(radiusvariables)]

    #convert to Vector{T} because may have missing values and would otherwise get Array{Union{Missing,Float64}} which optimization routines won't like.
    return Vector{T}(vcat(radiusvariables, thicknessvariables, conicvariables)), Vector{T}(vcat(radiuslower, thicknesslower, coniclower)), Vector{T}(vcat(radiusupper, thicknessupper, conicupper))
end

"""
    updateoptimizationvariables(a::AxisymmetricOpticalSystem{T}, optimizationvariables::Vector{S}) -> AxisymmetricOpticalSystem{S}

Creates a new optical system with updated variables corresponding to the optimization variables.
"""
function updateoptimizationvariables(a::AxisymmetricOpticalSystem{T}, optimizationvariables::Vector{S}) where {T<:Real,S<:Real}
    prescription = a.prescription
    #assume these two types of prescription variables will always be present
    prescradii = prescription[!, :Radius]
    prescthicknesses = prescription[!, :Thickness]

    @assert reduce(&, map((x) -> !isnan(x), optimizationvariables)) #make sure no optimization variables are NaN

    conic = in(:Conic, propertynames(prescription)) ? prescription[!, :Conic] : nothing

    radii = Vector{Union{Missing,S}}(undef, 0)
    thicknesses = Vector{Union{Missing,S}}(undef, 0)
    newconic = Vector{Union{Missing,S}}(undef, 0)

    optindex = 1
    if in(:OptimizeRadius, propertynames(prescription))
        radiusflags = prescription[!, :OptimizeRadius]
        for index in eachindex(radiusflags)
            if radiusflags[index]
                push!(radii, optimizationvariables[optindex])
                optindex += 1
            else
                push!(radii, prescradii[index])
            end
        end
    end

    if in(:OptimizeThickness, propertynames(prescription))
        thicknessflags = prescription[!, :OptimizeThickness]
        for index in eachindex(thicknessflags)
            if thicknessflags[index]
                push!(thicknesses, optimizationvariables[optindex])
                optindex += 1
            else
                push!(thicknesses, prescthicknesses[index])
            end
        end
    end

    if conic !== nothing
        if in(:OptimizeConic, propertynames(prescription))
            conicflags = prescription[!, :OptimizeConic]
            for index in eachindex(conicflags)
                if conicflags[index]
                    push!(newconic, optimizationvariables[optindex])
                    optindex += 1
                else
                    push!(newconic, conic[index])
                end
            end
        end
    end

    newprescription = copy(prescription)
    newprescription.Radius = radii
    newprescription.Thickness = thicknesses
    if conic !== nothing
        newprescription.Conic = newconic
    end

    return AxisymmetricOpticalSystem{S}(newprescription, detectorsize(a)..., S, temperature = S(temperature(a)), pressure = S(pressure(a)))
end

end #module
export Optimization
