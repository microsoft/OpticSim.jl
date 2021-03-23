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

using StaticArrays

@enum GlassType MODEL MIL AGF OTHER AIR

"""
Object identifying a glass, containing a type (e.g. `MODEL`, `MIL`, `OTHER` or `AGF`) depending on how the glass is defined, and an integer ID.
Air is `AIR:0`, others are on the form `AGF:N`, for example.
"""
struct GlassID
    type::GlassType
    num::Int
end

Base.show(io::IO, a::GlassID) = print(io, "$(string(a.type)):$(a.num)")
type(a::GlassID) = a.type
num(a::GlassID) = a.num

ismodel(a::GlassID) = a.type === MODEL
isMIL(a::GlassID) = a.type === MIL
isAGF(a::GlassID) = a.type === AGF
isother(a::GlassID) = a.type === OTHER

"""
Abstract type encapsulating all glasses.
"""
abstract type AbstractGlass end

"""
Stores all attributes relating to a glass type specified in an .AGF glass catalog.

Never used directly, instead created using catalog glasses, e.g. `GlassCat.SCHOTT.N_BK7`.

In order to prevent type ambiguities in OpticSim.jl we can't have this type paramaterized.
"""
struct Glass <: AbstractGlass
    ID::GlassID
    dispform::Int
    C1::Float64
    C2::Float64
    C3::Float64
    C4::Float64
    C5::Float64
    C6::Float64
    C7::Float64
    C8::Float64
    C9::Float64
    C10::Float64
    λmin::Float64
    λmax::Float64
    D₀::Float64
    D₁::Float64
    D₂::Float64
    E₀::Float64
    E₁::Float64
    λₜₖ::Float64
    temp::Float64
    ΔPgF::Float64
    PR::Float64
    relcost::Float64
    TCE::Float64
    CR::Float64
    status::Int
    SR::Float64
    transmission::Union{Nothing,SVector{100,SVector{3,Float64}}} # could use Vector on 1.5.X, but on 1.6 this causes allocs so we need to use a SVector...
    transmissionN::Int
    Nd::Float64
    AR::Float64
    FR::Float64
    exclude_sub::Int
    Vd::Float64
    ignore_thermal_exp::Int
    p::Float64
    meltfreq::Int

    function Glass(ID::GlassID, dispform, C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, λmin, λmax, D₀, D₁, D₂, E₀, E₁, λₜₖ, temp, ΔPgF, PR, relcost, TCE, CR, status, SR, transmission, Nd, AR, FR, exclude_sub, Vd, ignore_thermal_exp, p, meltfreq)
        # need a constructor to massage the transmission data
        if transmission === nothing
            transmission_s = nothing
            transmissionN = -1
        else
            fill = [SVector(0.0, 0.0, 0.0) for _ in 1:(100 - length(transmission))]
            transmission_s = SVector{100,SVector{3,Float64}}(transmission..., fill...)
            transmissionN = length(transmission)
        end
        return new(ID, dispform, C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, λmin, λmax, D₀, D₁, D₂, E₀, E₁, λₜₖ, temp, ΔPgF, PR, relcost, TCE, CR, status, SR, transmission_s, transmissionN, Nd, AR, FR, exclude_sub, Vd, ignore_thermal_exp, p, meltfreq)
    end
end

"""
    glassid(g::AbstractGlass) -> GlassID

Get the ID of the glass, see [`GlassID`](@ref).
"""
glassid(g::Glass) = g.ID

"""
    glassname(g::Union{AbstractGlass,GlassID})

Get the name (including catalog) of the glass, or glass with this ID.
"""
glassname(g::Glass) = glassname(g.ID)
function glassname(ID::GlassID)
    if ismodel(ID)
        return "GlassCat.ModelGlass.$(num(ID))"
    elseif isMIL(ID)
        return "GlassCat.GlassFromMIL.$(num(ID))"
    elseif isair(ID)
        return "GlassCat.Air"
    elseif isother(ID)
        return OTHER_GLASS_NAMES[num(ID)]
    else
        return AGF_GLASS_NAMES[num(ID)]
    end
end

function Base.show(io::IO, g::Glass)
    print(io, glassname(g))
end

"""
    glassforid(ID::GlassID)

Get the glass for a given ID.
"""
function glassforid(ID::GlassID)
    if ismodel(ID)
        return MODEL_GLASSES[num(ID)]::Glass
    elseif isMIL(ID)
        return MIL_GLASSES[num(ID)]::Glass
    elseif isair(ID)
        return GlassCat.Air
    elseif isother(ID)
        return OTHER_GLASSES[num(ID)]::Glass
    else
        return AGF_GLASSES[num(ID)]::Glass
    end
end

"""
    info([io::IO], glass::AbstractGlass)

Print out all data associated with `glass` in an easily readable format.

# Examples
```julia-repl
julia> info(GlassCat.RPO.IG4)
ID:                                                AGF:52
Dispersion formula:                                Schott (1)
Dispersion formula coefficients:
     a₀:                                           6.91189161
     a₁:                                           -0.000787956404
     a₂:                                           -4.22296071
     a₃:                                           142.900646
     a₄:                                           -1812.32748
     a₅:                                           7766.33028
Valid wavelengths:                                 3.0μm to 12.0μm
Reference temperature:                              20.0°C
Thermal ΔRI coefficients:
     D₀:                                           3.24e-5
     D₁:                                           0.0
     D₂:                                           0.0
     E₀:                                           0.0
     E₁:                                           0.0
     λₜₖ:                                          0.0
TCE (÷1e-6):                                       20.4
Ignore thermal expansion:                          false
Density (p):                                       4.47g/m³
ΔPgF:                                              0.0
RI at sodium D-Line (587nm):                       1.0
Abbe Number:                                       0.0
Cost relative to N_BK7:                              ?
Status:                                            Standard (0)
Melt frequency:                                    0
Exclude substitution:                              false
```
"""
function info(io::IO, glass::Glass)
    D = glass.dispform

    println(io, "$(rpad("ID:", 50)) $(glass.ID)")

    if (D == -2)
        println(io, "$(rpad("Dispersion formula:", 50)) Cauchy (-2)")
    elseif (D == -1)
        println(io, "$(rpad("Dispersion formula:", 50)) Fitted for model/MIL glass")
    else
        println(io, "$(rpad("Dispersion formula:", 50)) $(DISPFORM_NAMES[D]) ($D)")
    end
    println(io, "Dispersion formula coefficients:")
    if (D == -2) # Cauchy
        println(io, "     $(rpad("A:", 45)) $(glass.C1)")
        println(io, "     $(rpad("B:", 45)) $(glass.C2)")
        println(io, "     $(rpad("C:", 45)) $(glass.C3)")
        println(io, "     $(rpad("D:", 45)) $(glass.C4)")
        println(io, "     $(rpad("E:", 45)) $(glass.C5)")
        println(io, "     $(rpad("F:", 45)) $(glass.C6)")
    elseif (D == -1)
        println(io, "     $(rpad("C₀:", 45)) $(glass.C1)")
        println(io, "     $(rpad("C₁:", 45)) $(glass.C2)")
        println(io, "     $(rpad("C₂:", 45)) $(glass.C3)")
        println(io, "     $(rpad("C₃:", 45)) $(glass.C4)")
    elseif (D == 1) # Schott
        println(io, "     $(rpad("a₀:", 45)) $(glass.C1)")
        println(io, "     $(rpad("a₁:", 45)) $(glass.C2)")
        println(io, "     $(rpad("a₂:", 45)) $(glass.C3)")
        println(io, "     $(rpad("a₃:", 45)) $(glass.C4)")
        println(io, "     $(rpad("a₄:", 45)) $(glass.C5)")
        println(io, "     $(rpad("a₅:", 45)) $(glass.C6)")
    elseif (D == 2)  # Sellmeier1
        println(io, "     $(rpad("K₁:", 45)) $(glass.C1)")
        println(io, "     $(rpad("L₁:", 45)) $(glass.C2)")
        println(io, "     $(rpad("K₂:", 45)) $(glass.C3)")
        println(io, "     $(rpad("L₂:", 45)) $(glass.C4)")
        println(io, "     $(rpad("K₃:", 45)) $(glass.C5)")
        println(io, "     $(rpad("L₃:", 45)) $(glass.C6)")
    elseif (D == 3)  # Herzberger
        println(io, "     $(rpad("A:", 45)) $(glass.C1)")
        println(io, "     $(rpad("B:", 45)) $(glass.C2)")
        println(io, "     $(rpad("C:", 45)) $(glass.C3)")
        println(io, "     $(rpad("D:", 45)) $(glass.C4)")
        println(io, "     $(rpad("E:", 45)) $(glass.C5)")
        println(io, "     $(rpad("F:", 45)) $(glass.C6)")
    elseif (D == 4)  # Sellmeier2
        println(io, "     $(rpad("A:", 45)) $(glass.C1)")
        println(io, "     $(rpad("B₁:", 45)) $(glass.C2)")
        println(io, "     $(rpad("λ₁:", 45)) $(glass.C3)")
        println(io, "     $(rpad("B₂:", 45)) $(glass.C4)")
        println(io, "     $(rpad("λ₂:", 45)) $(glass.C5)")
    elseif (D == 5)  # Conrady
        println(io, "     $(rpad("n₀:", 45)) $(glass.C1)")
        println(io, "     $(rpad("A:", 45)) $(glass.C2)")
        println(io, "     $(rpad("B:", 45)) $(glass.C3)")
    elseif (D == 6)  # Sellmeier3
        println(io, "     $(rpad("K₁:", 45)) $(glass.C1)")
        println(io, "     $(rpad("L₁:", 45)) $(glass.C2)")
        println(io, "     $(rpad("K₂:", 45)) $(glass.C3)")
        println(io, "     $(rpad("L₂:", 45)) $(glass.C4)")
        println(io, "     $(rpad("K₃:", 45)) $(glass.C5)")
        println(io, "     $(rpad("L₃:", 45)) $(glass.C6)")
        println(io, "     $(rpad("K₄:", 45)) $(glass.C7)")
        println(io, "     $(rpad("L₄:", 45)) $(glass.C8)")
    elseif (D == 7) || (D == 8)  # HandbookOfOptics1/2
        println(io, "     $(rpad("A:", 45)) $(glass.C1)")
        println(io, "     $(rpad("B:", 45)) $(glass.C2)")
        println(io, "     $(rpad("C:", 45)) $(glass.C3)")
        println(io, "     $(rpad("D:", 45)) $(glass.C4)")
    elseif (D == 9)  # Sellmeier4
        println(io, "     $(rpad("A:", 45)) $(glass.C1)")
        println(io, "     $(rpad("B:", 45)) $(glass.C2)")
        println(io, "     $(rpad("C:", 45)) $(glass.C3)")
        println(io, "     $(rpad("D:", 45)) $(glass.C4)")
        println(io, "     $(rpad("E:", 45)) $(glass.C5)")
    elseif (D == 10) || (D == 12) # Extended1/2
        println(io, "     $(rpad("a₀:", 45)) $(glass.C1)")
        println(io, "     $(rpad("a₁:", 45)) $(glass.C2)")
        println(io, "     $(rpad("a₂:", 45)) $(glass.C3)")
        println(io, "     $(rpad("a₃:", 45)) $(glass.C4)")
        println(io, "     $(rpad("a₄:", 45)) $(glass.C5)")
        println(io, "     $(rpad("a₅:", 45)) $(glass.C6)")
        println(io, "     $(rpad("a₆:", 45)) $(glass.C7)")
        println(io, "     $(rpad("a₇:", 45)) $(glass.C8)")
    elseif (D == 11)  # Sellmeier5
        println(io, "     $(rpad("K₁:", 45)) $(glass.C1)")
        println(io, "     $(rpad("L₁:", 45)) $(glass.C2)")
        println(io, "     $(rpad("K₂:", 45)) $(glass.C3)")
        println(io, "     $(rpad("L₂:", 45)) $(glass.C4)")
        println(io, "     $(rpad("K₃:", 45)) $(glass.C5)")
        println(io, "     $(rpad("L₃:", 45)) $(glass.C6)")
        println(io, "     $(rpad("K₄:", 45)) $(glass.C7)")
        println(io, "     $(rpad("L₄:", 45)) $(glass.C8)")
        println(io, "     $(rpad("K₅:", 45)) $(glass.C9)")
        println(io, "     $(rpad("L₅:", 45)) $(glass.C10)")
    elseif (D == 13)  # Extended3
        println(io, "     $(rpad("a₀:", 45)) $(glass.C1)")
        println(io, "     $(rpad("a₁:", 45)) $(glass.C2)")
        println(io, "     $(rpad("a₂:", 45)) $(glass.C3)")
        println(io, "     $(rpad("a₃:", 45)) $(glass.C4)")
        println(io, "     $(rpad("a₄:", 45)) $(glass.C5)")
        println(io, "     $(rpad("a₅:", 45)) $(glass.C6)")
        println(io, "     $(rpad("a₆:", 45)) $(glass.C7)")
        println(io, "     $(rpad("a₇:", 45)) $(glass.C8)")
        println(io, "     $(rpad("a₈:", 45)) $(glass.C9)")
    else
        println(io, "     INVALID DISPERSION FORMULA!!")
    end
    println(io, "$(rpad("Valid wavelengths:", 50)) $(glass.λmin)μm to $(glass.λmax)μm")
    println(io, "$(rpad("Reference temperature:", 50)) $(glass.temp)°C")

    if !isnan(glass.D₀) && (glass.D₀ != 0 || glass.D₁ != 0 || glass.D₂ != 0 || glass.E₀ != 0 || glass.E₁ != 0)
        println(io, "Thermal ΔRI coefficients:")
        println(io, "     $(rpad("D₀:", 45)) $(glass.D₀)")
        println(io, "     $(rpad("D₁:", 45)) $(glass.D₁)")
        println(io, "     $(rpad("D₂:", 45)) $(glass.D₂)")
        println(io, "     $(rpad("E₀:", 45)) $(glass.E₀)")
        println(io, "     $(rpad("E₁:", 45)) $(glass.E₁)")
        println(io, "     $(rpad("λₜₖ:", 45)) $(glass.λₜₖ)")
    end

    println(io, "$(rpad("TCE (÷1e-6):", 50)) $(glass.TCE)")
    println(io, "$(rpad("Ignore thermal expansion:", 50)) $(glass.ignore_thermal_exp == 1)")

    println(io, "$(rpad("Density (p):", 50)) $(glass.p)g/m³")
    println(io, "$(rpad("ΔPgF:", 50)) $(glass.ΔPgF)")

    println(io, "$(rpad("RI at sodium D-Line (587nm):", 50)) $(glass.Nd)")
    println(io, "$(rpad("Abbe Number:", 50)) $(glass.Vd)")

    println(io, "$(rpad("Cost relative to N_BK7:", 50)) $(glass.relcost == -1 ? "?" : glass.relcost)")

    if glass.CR != -1 || glass.FR != -1 || glass.SR != -1 || glass.AR != -1 || glass.PR != -1
        println(io, "Environmental resistance:")
        println(io, "     $(rpad("Climate (CR):", 45)) $(glass.CR == -1 ? "?" : glass.CR)")
        println(io, "     $(rpad("Stain (FR):", 45)) $(glass.FR == -1 ? "?" : glass.FR)")
        println(io, "     $(rpad("Acid (SR):", 45)) $(glass.SR == -1 ? "?" : glass.SR)")
        println(io, "     $(rpad("Alkaline (AR):", 45)) $(glass.AR == -1 ? "?" : glass.AR)")
        println(io, "     $(rpad("Phosphate (PR):", 45)) $(glass.PR == -1 ? "?" : glass.PR)")
    end

    println(io, "$(rpad("Status:", 50)) $(STATUS[glass.status + 1]) ($(glass.status))")
    println(io, "$(rpad("Melt frequency:", 50)) $(glass.meltfreq == -1 ? "?" : glass.meltfreq)")
    println(io, "$(rpad("Exclude substitution:", 50)) $(glass.exclude_sub == 1)")

    if glass.transmission !== nothing
        println(io, "Transmission data:")
        println(io, "$(lpad("Wavelength", 15))$(lpad("Transmission", 15))$(lpad("Thickness", 15))")
        for i in 1:(glass.transmissionN)
            λ, t, τ = glass.transmission[i]
            println(io, "$(lpad("$(λ)μm", 15))$(lpad(t, 15))$(lpad("$(τ)mm", 15))")
        end
    end
end

info(g::AbstractGlass) = info(stdout, g)

function docstring(io::IO, glass::Glass)
    D = glass.dispform

    println(io, "**ID**$(glass.ID)\n")

    if (D == -2)
        println(io, "**Dispersion formula**Cauchy (-2)\n")
    elseif (D == -1)
        println(io, "**Dispersion formula**Fitted for model/MIL glass\n")
    else
        println(io, "**Dispersion formula**$(DISPFORM_NAMES[D]) ($D)\n")
    end
    println(io, "**Dispersion formula coefficients**\n")
    println(io, "|Coefficient|Value|")
    println(io, "|:----------|-----|")
    if (D == -2) # Cauchy
        println(io, "|A |$(glass.C1)|")
        println(io, "|B |$(glass.C2)|")
        println(io, "|C |$(glass.C3)|")
        println(io, "|D |$(glass.C4)|")
        println(io, "|E |$(glass.C5)|")
        println(io, "|F |$(glass.C6)|")
    elseif (D == -1)
        println(io, "|C₀|$(glass.C1)|")
        println(io, "|C₁|$(glass.C2)|")
        println(io, "|C₂|$(glass.C3)|")
        println(io, "|C₃|$(glass.C4)|")
    elseif (D == 1) # Schott
        println(io, "|a₀|$(glass.C1)|")
        println(io, "|a₁|$(glass.C2)|")
        println(io, "|a₂|$(glass.C3)|")
        println(io, "|a₃|$(glass.C4)|")
        println(io, "|a₄|$(glass.C5)|")
        println(io, "|a₅|$(glass.C6)|")
    elseif (D == 2)  # Sellmeier1
        println(io, "|K₁|$(glass.C1)|")
        println(io, "|L₁|$(glass.C2)|")
        println(io, "|K₂|$(glass.C3)|")
        println(io, "|L₂|$(glass.C4)|")
        println(io, "|K₃|$(glass.C5)|")
        println(io, "|L₃|$(glass.C6)|")
    elseif (D == 3)  # Herzberger
        println(io, "|A |$(glass.C1)|")
        println(io, "|B |$(glass.C2)|")
        println(io, "|C |$(glass.C3)|")
        println(io, "|D |$(glass.C4)|")
        println(io, "|E |$(glass.C5)|")
        println(io, "|F |$(glass.C6)|")
    elseif (D == 4)  # Sellmeier2
        println(io, "|A |$(glass.C1)|")
        println(io, "|B₁|$(glass.C2)|")
        println(io, "|λ₁|$(glass.C3)|")
        println(io, "|B₂|$(glass.C4)|")
        println(io, "|λ₂|$(glass.C5)|")
    elseif (D == 5)  # Conrady
        println(io, "|n₀|$(glass.C1)|")
        println(io, "|A |$(glass.C2)|")
        println(io, "|B |$(glass.C3)|")
    elseif (D == 6)  # Sellmeier3
        println(io, "|K₁|$(glass.C1)|")
        println(io, "|L₁|$(glass.C2)|")
        println(io, "|K₂|$(glass.C3)|")
        println(io, "|L₂|$(glass.C4)|")
        println(io, "|K₃|$(glass.C5)|")
        println(io, "|L₃|$(glass.C6)|")
        println(io, "|K₄|$(glass.C7)|")
        println(io, "|L₄|$(glass.C8)|")
    elseif (D == 7) || (D == 8)  # HandbookOfOptics1/2
        println(io, "|A |$(glass.C1)|")
        println(io, "|B |$(glass.C2)|")
        println(io, "|C |$(glass.C3)|")
        println(io, "|D |$(glass.C4)|")
    elseif (D == 9)  # Sellmeier4
        println(io, "|A |$(glass.C1)|")
        println(io, "|B |$(glass.C2)|")
        println(io, "|C |$(glass.C3)|")
        println(io, "|D |$(glass.C4)|")
        println(io, "|E |$(glass.C5)|")
    elseif (D == 10) || (D == 12) # Extended1/2
        println(io, "|a₀|$(glass.C1)|")
        println(io, "|a₁|$(glass.C2)|")
        println(io, "|a₂|$(glass.C3)|")
        println(io, "|a₃|$(glass.C4)|")
        println(io, "|a₄|$(glass.C5)|")
        println(io, "|a₅|$(glass.C6)|")
        println(io, "|a₆|$(glass.C7)|")
        println(io, "|a₇|$(glass.C8)|")
    elseif (D == 11)  # Sellmeier5
        println(io, "|K₁|$(glass.C1)|")
        println(io, "|L₁|$(glass.C2)|")
        println(io, "|K₂|$(glass.C3)|")
        println(io, "|L₂|$(glass.C4)|")
        println(io, "|K₃|$(glass.C5)|")
        println(io, "|L₃|$(glass.C6)|")
        println(io, "|K₄|$(glass.C7)|")
        println(io, "|L₄|$(glass.C8)|")
        println(io, "|K₅|$(glass.C9)|")
        println(io, "|L₅|$(glass.C10)|")
    elseif (D == 13)  # Extended3
        println(io, "|a₀|$(glass.C1)|")
        println(io, "|a₁|$(glass.C2)|")
        println(io, "|a₂|$(glass.C3)|")
        println(io, "|a₃|$(glass.C4)|")
        println(io, "|a₄|$(glass.C5)|")
        println(io, "|a₅|$(glass.C6)|")
        println(io, "|a₆|$(glass.C7)|")
        println(io, "|a₇|$(glass.C8)|")
        println(io, "|a₈|$(glass.C9)|")
    else
        println(io, "     INVALID DISPERSION FORMULA!!\n")
    end
    println(io)

    println(io, "**Valid wavelengths**$(glass.λmin)μm to $(glass.λmax)μm\n")
    println(io, "**Reference temperature**$(glass.temp)°C\n")

    if !isnan(glass.D₀) && (glass.D₀ != 0 || glass.D₁ != 0 || glass.D₂ != 0 || glass.E₀ != 0 || glass.E₁ != 0)
        println(io, "**Thermal ΔRI coefficients**\n")
        println(io, "|Coefficient                 |Value      |")
        println(io, "|:---------------------------|-----------|")
        println(io, "|D₀                          |$(glass.D₀)|")
        println(io, "|D₁                          |$(glass.D₁)|")
        println(io, "|D₂                          |$(glass.D₂)|")
        println(io, "|E₀                          |$(glass.E₀)|")
        println(io, "|E₁                          |$(glass.E₁)|")
        println(io, "|λₜₖ                          |$(glass.λₜₖ)|")
        println(io)
    end

    println(io, "**TCE (÷1e-6)**$(glass.TCE)\n")
    println(io, "**Ignore thermal expansion**$(glass.ignore_thermal_exp == 1)\n")

    println(io, "**Density (p)**$(glass.p)g/m³\n")
    println(io, "**ΔPgF**$(glass.ΔPgF)\n")

    println(io, "**RI at sodium D-Line (587nm)**$(glass.Nd)\n")
    println(io, "**Abbe Number**$(glass.Vd)\n")

    println(io, "**Cost relative to N_BK7**$(glass.relcost == -1 ? "?" : glass.relcost)\n")

    if glass.CR != -1 || glass.FR != -1 || glass.SR != -1 || glass.AR != -1 || glass.PR != -1
        println(io, "**Environmental resistance**\n")
        println(io, "|Environment   |Resistance                        |")
        println(io, "|:-------------|----------------------------------|")
        println(io, "|Climate (CR)  |$(glass.CR == -1 ? "?" : glass.CR)|")
        println(io, "|Stain (FR)    |$(glass.FR == -1 ? "?" : glass.FR)|")
        println(io, "|Acid (SR)     |$(glass.SR == -1 ? "?" : glass.SR)|")
        println(io, "|Alkaline (AR) |$(glass.AR == -1 ? "?" : glass.AR)|")
        println(io, "|Phosphate (PR)|$(glass.PR == -1 ? "?" : glass.PR)|")
        println(io)
    end

    println(io, "**Status**$(STATUS[glass.status + 1]) ($(glass.status))\n")
    println(io, "**Melt frequency**$(glass.meltfreq == -1 ? "?" : glass.meltfreq)\n")
    println(io, "**Exclude substitution**$(glass.exclude_sub == 1)\n")

    if glass.transmission !== nothing
        println(io, "**Transmission data**\n")
        println(io, "|Wavelength|Transmission|Thickness|")
        println(io, "|----------|------------|---------|")
        for i in 1:(glass.transmissionN)
            λ, t, τ = glass.transmission[i]
            println(io, "|$(λ)μm|$t|$(τ)mm|")
        end
    end
end
