# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

# runtime glass cats
const MIL_GLASSES = Dict{Int,Glass}()
const MODEL_GLASSES = Vector{Glass}(undef, 0)

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
