# MIT License

# Copyright (c) Microsoft Corporation.

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE

function _modelglass(Nd::Float64, Vd::Float64, ΔPgF::Float64, ID::GlassID)
    # from Schott "TIE-29: Refractive Index and Dispersion"
    a = ΔPgF + 0.6438 - 0.001682 * Vd
    # Using fitting results from https://www.gnu.org/software/goptical/manual/Material_Abbe_class_reference.html
    C1 = a * -6.11873891971188577088 + 1.17752614766485175224
    C2 = a * 18.27315722388047447566 + -8.93204522498095698779
    C3 = a * -14.55275321129051135927 + 7.91015964461522003148
    C4 = a * 3.48385106908642905310 + -1.80321117937358499361
    return Glass(ID, -1, C1, C2, C3, C4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.36, 0.75, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, TEMP_REF, ΔPgF, -1.0, -1.0, 0.0, -1.0, 0, -1.0, nothing, Nd, -1.0, -1.0, 0, Vd, 1, 0.0, 0.0)
end

"""
    glassfromMIL(glasscode::Union{Float64,Int}) -> Glass

Generates a glass object for the given glass code based on U.S. military standard MIL-G-174, see [the MIL specification](http://www.newportglass.com/GeneCd.htm) for further details.

The glass code is a six-digit number specifying the glass according to its refractive index `Nd` at d-light (587.5618nm), and its Abbe number `Vd` also taken at d-light.
The resulting glass code is the value of `Nd - 1` rounded to three digits, followed by `Vd` rounded to three digits, with all decimal points ignored.
For example, N_BK7 has `Nd = 1.5168` and `Vd = 64.17`, giving a six-digit glass code of `517642`.

For `Nd > 1.999` the format `1.123642` can be used representing `Nd = 2.123` and `Vd = 64.2`.

**Accuracy is poor given the low precision of the input parameters**, the mean error to measured data may be significant.
Behavior may differ from other optical simulation tools when using MIL glasses.
The approximate dispersion calculation used these glasses is generally only valid for visible wavelengths, in this case a limit of 360nm to 750nm is imposed.

# Examples
```julia-repl
julia> index(glassfromMIL(517642), 0.5875618)
1.5170003960064509

julia> index(glassfromMIL(1.134642), 0.5875618)
2.1340008686098946
```
"""
function glassfromMIL(glasscode::Int)
    Nd = floor(Int, glasscode / 1000) / 1000 + 1
    Vd = (glasscode - floor(Int, glasscode / 1000) * 1000) / 10
    g = _modelglass(Nd, Vd, 0.0, GlassID(MIL, glasscode))
    MIL_GLASSES[glasscode] = g
    return g
end
function glassfromMIL(glasscode::Float64)
    @assert glasscode > 1.0
    glasscodeid = round(Int, glasscode * 1000000)
    Nd = floor(Int, glasscode * 1000) / 1000 + 1
    Vd = round((glasscode * 1000 - floor(Int, glasscode * 1000)) * 100, digits = 1)
    g = _modelglass(Nd, Vd, 0.0, GlassID(MIL, glasscodeid))
    MIL_GLASSES[glasscodeid] = g
    return g
end

"""
    modelglass(Nd::Float64, Vd::Float64, ΔPgF::Float64) -> Glass

Generates a glass object for the given refractive index at d-light (587.5618nm), `Nd`, the Abbe number also at d-light, `Vd`, and partial dispersion, `ΔPgF`.
The mean error to measured data for these models is typically small - usually < 0.0001.
Behavior may differ from other optical simulation tools when using model glasses.

The approximate dispersion calculation used for these glasses is generally only valid for visible wavelengths, in this case a limit of 360nm to 750nm is imposed.

# Examples
```julia-repl
julia> index(modelglass(1.5168, 64.17, 0.0), 0.5875618)
1.5168003970108495

julia> index(modelglass(1.2344, 61.57, 0.003), 0.678)
1.2329425902693352
```
"""
function modelglass(Nd::Float64, Vd::Float64, ΔPgF::Float64)
    g = _modelglass(Nd, Vd, ΔPgF, GlassID(MODEL, length(MODEL_GLASSES) + 1))
    push!(MODEL_GLASSES, g)
    return g
end

"""
    absorption(glass::AbstractGlass, wavelength; temperature=20°C, pressure=1Atm)

Compute the intensity absorption per mm of `glass` at `wavelength`, optionally at specified `temperature` and `pressure`.
Transmission values are linearly interpolated from the adjacent values in the data table of `glass`, if `wavelength` is below the minimum or above the maximum in the table then the nearest value is taken.

Absorption is defined as ``\\frac{-\\log(t)}{\\tau}`` where ``t`` is the transmission value and ``\\tau`` is the thickness, both of which are provided in the data table.

If unitless, arguments are interpretted as μm, °C and Atm respectively.

# Examples
```julia-repl
julia> absorption(GlassCat.Sumita.LAK7, 700u"nm")
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
julia> index(GlassCat.Sumita.LAK7, 700u"nm")
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
