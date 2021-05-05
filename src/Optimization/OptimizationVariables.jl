# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

"""
    optimizationvariables(a::AxisymmetricOpticalSystem{T}) -> Vector{T}

Pack variables that have been marked to be optimized into a vector in a form suitable for the optimizer. Variables are
marked for optimization by having a true value in the `:OptimizeName` column, where `Name` can be `Radius`, `Thickness`
or `Conic`.
"""
function optimizationvariables(a::AxisymmetricOpticalSystem{T}) where {T<:Real}
    prescription = a.prescription

    radiusvariables = (:OptimizeRadius ∈ propertynames(prescription) ?
        prescription[prescription[!, :OptimizeRadius], :Radius] : [])
    radiuslower = repeat([typemin(T)], length(radiusvariables))
    radiusupper = repeat([typemax(T)], length(radiusvariables))

    thicknessvariables = (:OptimizeThickness ∈ propertynames(prescription) ?
        prescription[prescription[!, :OptimizeThickness], :Thickness] : [])
    thicknesslower = repeat([zero(T)], length(thicknessvariables))
    thicknessupper = repeat([convert(T, 100)], length(thicknessvariables))

    conicvariables = (:OptimizeConic ∈ propertynames(prescription) ?
        prescription[prescription[!, :OptimizeConic], :Conic] : [])
    coniclower = repeat([typemin(T)], length(radiusvariables))
    conicupper = repeat([typemax(T)], length(radiusvariables))

    # convert to Vector{T} because may have missing values and would otherwise get Array{Union{Missing,Float64}} which
    # optimization routines won't like.
    return (
        Vector{T}(vcat(radiusvariables, thicknessvariables, conicvariables)),
        Vector{T}(vcat(radiuslower, thicknesslower, coniclower)),
        Vector{T}(vcat(radiusupper, thicknessupper, conicupper)))
end

"""
    updateoptimizationvariables(a::AxisymmetricOpticalSystem{T}, optimizationvariables::Vector{S}) -> AxisymmetricOpticalSystem{S}

Creates a new optical system with updated variables corresponding to the optimization variables.
"""
function updateoptimizationvariables(a::AxisymmetricOpticalSystem{T}, optimizationvariables::Vector{S}) where {T<:Real,S<:Real}
    prescription = a.prescription
    # assume these two types of prescription variables will always be present
    prescradii = prescription[!, :Radius]
    prescthicknesses = prescription[!, :Thickness]

    @assert reduce(&, map((x) -> !isnan(x), optimizationvariables)) # make sure no optimization variables are NaN

    conic = :Conic ∈ propertynames(prescription) ? prescription[!, :Conic] : nothing

    radii = Vector{Union{Missing,S}}(undef, 0)
    thicknesses = Vector{Union{Missing,S}}(undef, 0)
    newconic = Vector{Union{Missing,S}}(undef, 0)

    optindex = 1
    if :OptimizeRadius ∈ propertynames(prescription)
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

    if :OptimizeThickness ∈ propertynames(prescription)
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

    if :OptimizeConic ∈ propertynames(prescription) && !isnothing(conic)
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

    newprescription = copy(prescription)
    newprescription.Radius = radii
    newprescription.Thickness = thicknesses
    if !isnothing(conic)
        newprescription.Conic = newconic
    end

    return AxisymmetricOpticalSystem{S}(
        newprescription,
        detectorsize(a)...,
        S,
        temperature = convert(S, temperature(a)),
        pressure = convert(S, pressure(a)))
end
