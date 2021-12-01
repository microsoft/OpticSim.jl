# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

"""
    absorption(glass::AbstractGlass, wavelength; temperature=20°C, pressure=1Atm)

Compute the intensity absorption per mm of `glass` at `wavelength`, optionally at specified `temperature` and `pressure`.
Transmission values are linearly interpolated from the adjacent values in the data table of `glass`, if `wavelength` is below the minimum or above the maximum in the table then the nearest value is taken.

Absorption is defined as ``\\frac{-\\log(t)}{\\tau}`` where ``t`` is the transmission value and ``\\tau`` is the thickness, both of which are provided in the data table.

If unitless, arguments are interpretted as μm, °C and Atm respectively.

# Examples
```julia-repl
julia> absorption(GlassCat.SUMITA.LAK7, 700u"nm")
0.0006018072325563021

julia> absorption(GlassCat.SCHOTT.N_BK7, 0.55, temperature = 22.0)
0.00016504471175660636

julia> absorption(GlassCat.SCHOTT.PSK3, 532u"nm", temperature = 25u"°C", pressure = 1.3)
0.00020855284788532435
```
"""
function absorption(glass::Glass, wavelength::Length; temperature::Temperature = TEMP_REF_UNITFUL, pressure::Float64 = PRESSURE_REF)::Float64
    λ = Float64(ustrip(u"μm", wavelength))
    return absorption(glass, λ, temperature = ustrip(Float64, u"°C", temperature), pressure = pressure)
end

function absorption(glass::Glass, λ::T; temperature::T = T(TEMP_REF), pressure::T = T(PRESSURE_REF))::T where {T<:Real}
    # if the glass has no transmission data then assume no absorption
    if glass.transmission === nothing
        return zero(T)
    end

    reference_temp = T(glass.temp)

    # to work out the wavelength at the reference temperature we need the RIs of air at system temp and at reference temp
    n_air_at_sys = absairindex(λ, temperature = temperature, pressure = pressure)
    n_air_at_ref = absairindex(λ, temperature = reference_temp)

    # scale the wavelength to air at the reference temperature/pressure
    λ = λ * (n_air_at_sys / n_air_at_ref)

    tdata = glass.transmission
    N = glass.transmissionN
    if λ < tdata[1][1]
        t = tdata[1][2]
        τ = tdata[1][3]
        return T(-log1p(t - 1.0) / τ)
    elseif λ > tdata[N][1]
        t = tdata[N][2]
        τ = tdata[N][3]
        return T(-log1p(t - 1.0) / τ)
    else
        let λlow = 0.0, tlow = 0.0, τlow = 0.0, λhigh = 0.0, thigh = 0.0, τhigh = 0.0
            for i in 2:N
                if λ <= tdata[i][1]
                    λlow, tlow, τlow = tdata[i - 1]
                    λhigh, thigh, τhigh = tdata[i]
                    break
                end
            end
            λhigh = T(λhigh)
            λlow = T(λlow)
            δλ = λhigh - λlow
            @assert τlow == τhigh
            t = (tlow * (λhigh - λ) / δλ) + (thigh * (λ - λlow) / δλ)
            return -log1p(t - 1.0) / τhigh
        end
    end
end

function absorption(::AirType, ::Length; temperature::Temperature = TEMP_REF_UNITFUL, pressure::Float64 = PRESSURE_REF)::Float64
    return 0.0
end

function absorption(::AirType, ::T; temperature::T = T(TEMP_REF), pressure::T = T(PRESSURE_REF))::T where {T<:Real}
    return zero(T)
end

"""
    index(glass::AbstractGlass, wavelength; temperature=20°C, pressure=1Atm)

Compute the refractive index of `glass` at `wavelength`, optionally at specified `temperature` and `pressure`.
Result is relative to the refractive index of air at given temperature and pressure.

If unitless, arguments are interpretted as μm, °C and Atm respectively.

**This is defined to always equal 1.0 for Air at any temperature and pressure**, use [`absairindex`](@ref) for the absolute refractive index of air at a given temperature and pressure.

# Examples
```julia-repl
julia> index(GlassCat.SUMITA.LAK7, 700u"nm")
1.646494204478318

julia> index(GlassCat.SCHOTT.N_BK7, 0.55, temperature = 22.0)
1.51852824383283

julia> index(GlassCat.HOYA.FF1, 532u"nm", temperature = 25u"°C", pressure = 1.3)
1.5144848290944655
```
"""
function index(glass::Glass, wavelength::Length; temperature::Temperature = TEMP_REF_UNITFUL, pressure::Float64 = PRESSURE_REF)::Float64
    λ = Float64(ustrip(uconvert(u"μm", wavelength)))
    return index(glass, λ, temperature = ustrip(Float64, u"°C", temperature), pressure = pressure)
end

function index(glass::Glass, λ::T; temperature::T = T(TEMP_REF), pressure::T = T(PRESSURE_REF))::T where {T<:Real}
    # all calculations for the material must be done at the refernce temperature
    reference_temp = T(glass.temp)

    # to work out the wavelength at the reference temperature we need the RIs of air at system temp and at reference temp
    n_air_at_sys = absairindex(λ, temperature = temperature, pressure = pressure)
    n_air_at_ref = absairindex(λ, temperature = reference_temp)

    # scale the wavelength to air at the reference temperature/pressure
    λabs = λ * n_air_at_sys
    λ = λabs / n_air_at_ref

    if (λ < glass.λmin) || (λ > glass.λmax)
        error("Cannot calculate an index for the specified wavelength: $λ, valid range: [$(glass.λmin), $(glass.λmax)].\n")
    end

    if glass.dispform == -2
        # Cauchy
        n_rel = T(glass.C1) + (glass.C2 * λ^(-2)) + (glass.C3 * λ^(-4)) + (glass.C4 * λ^(-6)) + (glass.C5 * λ^(-8)) + (glass.C6 * λ^(-10))
    elseif glass.dispform == -1
        # use fitted result from GOptical:
        n_rel = T(glass.Nd + (glass.Nd - one(T)) / glass.Vd * (glass.C1 + glass.C2 / λ + glass.C3 / λ^2 + glass.C4 / λ^3))
    elseif glass.dispform == 1
        # Schott
        formula_rhs = T(glass.C1) + (glass.C2 * λ^2) + (glass.C3 * λ^(-2)) + (glass.C4 * λ^(-4)) + (glass.C5 * λ^(-6)) + (glass.C6 * λ^(-8))
        n_rel = sqrt(formula_rhs)
    elseif glass.dispform == 2
        # Sellmeier1
        formula_rhs = (glass.C1 * λ^2 / (λ^2 - glass.C2)) + (glass.C3 * λ^2 / (λ^2 - glass.C4)) + (glass.C5 * λ^2 / (λ^2 - glass.C6))
        n_rel = sqrt(formula_rhs + one(T))
    elseif glass.dispform == 3
        # Herzberger
        L = one(T) / (λ^2 - T(0.028))
        n_rel = T(glass.C1) + (glass.C2 * L) + (glass.C3 * L^2) + (glass.C4 * λ^2) + (glass.C5 * λ^4) + (glass.C6 * λ^6)
    elseif glass.dispform == 4
        # Sellmeier2
        formula_rhs = T(glass.C1) + (glass.C2 * λ^2 / (λ^2 - (glass.C3)^2)) + (glass.C4 * λ^2 / (λ^2 - (glass.C5)^2))
        n_rel = sqrt(formula_rhs + one(T))
    elseif glass.dispform == 5
        # Conrady
        n_rel = T(glass.C1) + (glass.C2 / λ) + (glass.C3 / λ^3.5)
    elseif glass.dispform == 6
        # Sellmeier3
        formula_rhs = (glass.C1 * λ^2 / (λ^2 - glass.C2)) + (glass.C3 * λ^2 / (λ^2 - glass.C4)) + (glass.C5 * λ^2 / (λ^2 - glass.C6)) + (glass.C7 * λ^2 / (λ^2 - glass.C8))
        n_rel = sqrt(formula_rhs + one(T))
    elseif glass.dispform == 7
        # HandbookOfOptics1
        formula_rhs = T(glass.C1) + (glass.C2 / (λ^2 - glass.C3)) - (glass.C4 * λ^2)
        n_rel = sqrt(formula_rhs)
    elseif glass.dispform == 8
        # HandbookOfOptics2
        formula_rhs = T(glass.C1) + (glass.C2 * λ^2 / (λ^2 - glass.C3)) - (glass.C4 * λ^2)
        n_rel = sqrt(formula_rhs)
    elseif glass.dispform == 9
        # Sellmeier4
        formula_rhs = T(glass.C1) + (glass.C2 * λ^2 / (λ^2 - glass.C3)) + (glass.C4 * λ^2 / (λ^2 - glass.C5))
        n_rel = sqrt(formula_rhs)
    elseif glass.dispform == 10
        # Extended1
        formula_rhs = T(glass.C1) + (glass.C2 * λ^2) + (glass.C3 * λ^(-2)) + (glass.C4 * λ^(-4)) + (glass.C5 * λ^(-6)) + (glass.C6 * λ^(-8)) + (glass.C7 * λ^(-10)) + (glass.C8 * λ^(-12))
        n_rel = sqrt(formula_rhs)
    elseif glass.dispform == 11
        # Sellmeier5
        formula_rhs = (glass.C1 * λ^2 / (λ^2 - glass.C2)) + (glass.C3 * λ^2 / (λ^2 - glass.C4)) + (glass.C5 * λ^2 / (λ^2 - glass.C6)) + (glass.C7 * λ^2 / (λ^2 - glass.C8)) + (glass.C9 * λ^2 / (λ^2 - glass.C10))
        n_rel = sqrt(formula_rhs + one(T))
    elseif glass.dispform == 12
        # Extended2
        formula_rhs = T(glass.C1) + (glass.C2 * λ^2) + (glass.C3 * λ^(-2)) + (glass.C4 * λ^(-4)) + (glass.C5 * λ^(-6)) + (glass.C6 * λ^(-8)) + (glass.C7 * λ^4) + (glass.C8 * λ^6)
        n_rel = sqrt(formula_rhs)
    elseif glass.dispform == 13
        # Extended3
        formula_rhs = T(glass.C1) + (glass.C2 * λ^2) + (glass.C3 * λ^(4)) + (glass.C4 * λ^(-2)) + (glass.C5 * λ^(-4)) + (glass.C6 * λ^(-6)) + (glass.C7 * λ^(-8)) + (glass.C8 * λ^(-10)) + (glass.C9 * λ^(-12))
        n_rel = sqrt(formula_rhs)
    else
        @error "Invalid glass dispersion formula"
    end

    # get the absolute index of the material
    n_abs = n_rel * n_air_at_ref

    # If "TD" is included in the glass data, then include pressure and temperature dependence of the lens
    # environment. From Schott"s technical report "TIE-19: Temperature Coefficient of the Refractive Index".
    # The above "n_rel" data are assumed to be from the reference temperature T_ref. Now we add a small change
    # delta_n to it due to a change in temperature.
    ΔT = temperature - reference_temp
    if !isnan(glass.D₀) && abs(ΔT) > 0.0 && (glass.D₀ != 0 || glass.D₁ != 0 || glass.D₂ != 0 || glass.E₀ != 0 || glass.E₁ != 0)
        Sₜₖ = glass.λₜₖ < 0.0 ? -one(T) : one(T)
        Δn_abs = ((n_rel^2 - one(T)) / (2.0 * n_rel)) * (glass.D₀ * ΔT + glass.D₁ * ΔT^2 + glass.D₂ * ΔT^3 + ((glass.E₀ * ΔT + glass.E₁ * ΔT^2) / (λ^2 - Sₜₖ * glass.λₜₖ^2)))
        n_abs = n_abs + Δn_abs
    end

    # make the index relative to the RI of the air at the system temperature/pressure again
    n_rel = n_abs / n_air_at_sys
    return n_rel
end

function index(::AirType, ::Length; temperature::Temperature = TEMP_REF_UNITFUL, pressure::Float64 = PRESSURE_REF)::Float64
    return 1.0
end

function index(::AirType, ::T; temperature::T = T(TEMP_REF), pressure::T = T(PRESSURE_REF))::T where {T<:Real}
    return one(T)
end

"""
    absairindex(wavelength; temperature=20°C, pressure=1Atm)

Compute the absolute refractive index of air at `wavelength`, optionally at specified `temperature` and `pressure`. If unitless, arguments are interpretted as μm, °C and Atm respectively.

# Examples
```julia-repl
julia> absairindex(700u"nm")
1.000271074905147

julia> absairindex(0.7, temperature=27.0)
1.000264738846504

julia> absairindex(532u"nm", temperature = 25u"°C", pressure = 1.3)
1.0003494991178161
```
"""
function absairindex(wavelength::Length; temperature::Temperature = TEMP_REF_UNITFUL, pressure::Float64 = PRESSURE_REF)::Float64
    # convert to required units
    λ = Float64(ustrip(uconvert(u"μm", wavelength)))
    return absairindex(λ, temperature = ustrip(Float64, u"°C", temperature), pressure = pressure)
end

function absairindex(λ::T; temperature::T = T(TEMP_REF), pressure::T = T(PRESSURE_REF))::T where {T<:Real}
    # convert to required units
    n_ref = one(T) + ((6432.8 + ((2949810.0 * λ^2) / (146.0 * λ^2 - one(T))) + ((25540.0 * λ^2) / (41.0 * λ^2 - one(T)))) * 1e-8)
    n_rel = one(T) + ((n_ref - one(T)) / (one(T) + (temperature - 15.0) * 0.0034785)) * (pressure / PRESSURE_REF)
    return n_rel
end

"""
    polyfit_indices(wavelengths, n_rel; degree=5)

Fit a polynomial to `indices` at `wavelengths`, optionally specifying the `degree` of the polynomial.
Returns tuple of array of fitted indices at wavelengths and the polynomial.
"""
function polyfit_indices(wavelengths::Union{AbstractRange{<:Length},AbstractArray{<:Length,1}}, indices::AbstractArray{<:Number,1}; degree::Int = 5)
    w = ustrip.(uconvert.(u"μm", wavelengths))
    okay = (indices .> 0.0)
    if !any(okay)
        return (ones(Float64, size(w)) .* NaN, nothing)
    end
    xs = range(-1.0, stop = 1.0, length = length(w[okay]))
    poly = fit(xs, indices[okay], degree)
    interp_indices = poly.(xs)
    # ensure output has all entries
    out = ones(Float64, size(w)) .* NaN
    out[okay] = interp_indices
    return (out, poly)
end

"""
    plot_indices(glass::AbstractGlass; polyfit=false, fiterror=false, degree=5, temperature=20°C, pressure=1Atm, nsamples=300, sampling_domain="wavelength")

Plot the refractive index for `glass` for `nsamples` within its valid range of wavelengths, optionally at `temperature` and `pressure`.
`polyfit` will show a polynomial of optionally specified `degree` fitted to the data, `fiterror` will also show the fitting error of the result.
`sampling_domain` specifies whether the samples will be spaced uniformly in "wavelength" or "wavenumber".
"""
function plot_indices(glass::AbstractGlass; polyfit::Bool = false, fiterror::Bool = false, degree::Int = 5, temperature::Temperature = TEMP_REF_UNITFUL, pressure::Float64 = PRESSURE_REF, nsamples::Int = 300, sampling_domain::String = "wavelength")
    if isair(glass)
        wavemin = 380 * u"nm"
        wavemax = 740 * u"nm"
    else
        wavemin = glass.λmin * u"μm"
        wavemax = glass.λmax * u"μm"
    end

    if (sampling_domain == "wavelength")
        waves = range(wavemin, stop = wavemax, length = nsamples)      # wavelength in um
    elseif (sampling_domain == "wavenumber")
        sigma_min = 1.0 / wavemax
        sigma_max = 1.0 / wavemin
        wavenumbers = range(sigma_min, stop = sigma_max, length = nsamples) # wavenumber in um.^-1
        waves = 1.0 ./ wavenumbers
    else
        error("Invalid sampling domain, should be \"wavelength\" or \"wavenumber\"")
    end

    p = plot(xlabel = "wavelength (um)", ylabel = "refractive index")

    f = w -> begin
        try
            return index(glass, w, temperature = temperature, pressure = pressure)
        catch
            return NaN
        end
    end
    indices = [f(w) for w in waves]
    plot!(ustrip.(waves), indices, color = :blue, label = "From Data")

    if polyfit
        (p_indices, _) = polyfit_indices(waves, indices, degree = degree)
        plot!(ustrip.(waves), p_indices, color = :black, markersize = 4, label = "Polyfit")
    end

    if polyfit && fiterror
        err = p_indices - indices
        p2 = plot(xlabel = "wavelength (um)", ylabel = "fit error")
        plot!(ustrip.(waves), err, color = :red, label = "Fit Error")
        p = plot(p, p2, layout = 2)
    end

    plot!(title = "$(glassname(glass)) dispersion")

    gui(p)
end


"""
    drawglassmap(glasscatalog::Module; λ::Length = 550nm, glassfontsize::Integer = 3, showprefixglasses::Bool = false)

Draw a scatter plot of index vs dispersion (the derivative of index with respect to wavelength). Both index and
dispersion are computed at wavelength λ.

Choose glasses to graph using the glassfilterprediate argument. This is a function that receives a Glass object and returns true if the glass should be graphed.

If showprefixglasses is true then glasses with names like `F_BAK7` will be displayed. Otherwise glasses that have a
leading letter prefix followed by an underscore, such as `F_`, will not be displayed.

The index formulas for some glasses may give incorrect results if λ is outside the valid range for that glass. This can
give anomalous results, such as indices less than zero or greater than 6. To filter out these glasses set maximumindex
to a reasonable value such as 3.0.

example: plot only glasses that do not contain the strings "E_" and "J_"

drawglassmap(NIKON,showprefixglasses = true,glassfilterpredicate = (x) -> !occursin("J_",string(x)) && !occursin("E_",string(x)))
"""
function drawglassmap(glasscatalog::Module; λ::Length = 550nm, glassfontsize::Integer = 3, showprefixglasses::Bool = false, minindex = 1.0, maxindex = 3.0, mindispersion = -.3, maxdispersion = 0.0, glassfilterpredicate = (x)->true)
    wavelength = Float64(ustrip(uconvert(μm, λ)))
    indices = Vector{Float64}(undef,0)
    dispersions = Vector{Float64}(undef,0)
    glassnames = Vector{String}(undef,0)

    for name in names(glasscatalog)
        glass = eval(:($glasscatalog.$name))
        glassstring = String(name)
        hasprefix = occursin("_", glassstring)
 
        if typeof(glass) !== Module && (minindex <= index(glass, wavelength) <= maxindex)
            f(x) = index(glass,x)
            g = x -> ForwardDiff.derivative(f, x);
            dispersion = g(wavelength)

            # don't show glasses that have an _ in the name. This prevents cluttering the map with many glasses of
            # similar (index, dispersion).
            if glassfilterpredicate(glass) && (mindispersion <= dispersion <= maxdispersion) && (showprefixglasses || !hasprefix)
                push!(indices, index(glass, wavelength))
                push!(dispersions, dispersion)
                push!(glassnames, String(name))
            end
        end
    end

    font = Plots.font(family = "Sans", pointsize = glassfontsize, color = RGB(0.0,0.0,.4))
    series_annotations = Plots.series_annotations(glassnames, font)
    scatter(
        dispersions,
        indices;
        series_annotations,
        markeralpha = 0.0,
        legends = :none,
        xaxis = "dispersion @$λ",
        yaxis = "index",
        title = "Glass Catalog: $glasscatalog",
        xflip = true) #should use markershape = :none to prevent markers from being drawn but this option doesn't work. Used markeralpha = 0 so the markers are invisible. A hack which works.
end
