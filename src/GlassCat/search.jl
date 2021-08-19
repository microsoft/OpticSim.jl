# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

"""
    glasscatalogs()

Returns the complete list of glass catalogs available from GlassCat.

## Example
```julia-repl
julia> glasscatalogs()
41-element Array{Any,1}:
 OpticSim.GlassCat.AMTIR
 OpticSim.GlassCat.ANGSTROMLINK
 OpticSim.GlassCat.APEL
 OpticSim.GlassCat.ARCHER
 OpticSim.GlassCat.ARTON
 OpticSim.GlassCat.AUER_LIGHTING
 OpticSim.GlassCat.BIREFRINGENT
 ⋮
```
"""
glasscatalogs() = _child_modules(GlassCat)

"""
    glassnames(catalog::Module)

Returns the glass names available from a given catalog.

# Example
```julia-repl
julia> glassnames(GlassCat.CARGILLE)
3-element Array{Any,1}:
 "OG0607"
 "OG0608"
 "OG081160"
```
"""
function glassnames(catalog::Module)
    glass_names = names(catalog, all = true, imported = false)
    glasses = []
    for glass_name in glass_names
        glass_name_str = string(glass_name)
        if !occursin("#", glass_name_str) && glass_name_str != "eval" && glass_name_str != "include" && glass_name != nameof(catalog)
            push!(glasses, glass_name)
        end
    end
    return glasses
end

"""
    glassnames()

Returns the glass names available from all catalogs.

# Example
```julia-repl
julia> glassnames()
6-element Array{Pair{Module,Array{Any,1}},1}:
 OpticSim.GlassCat.CARGILLE => ["OG0607", "OG0608", "OG081160"]
     OpticSim.GlassCat.HOYA => ["BAC4", "BACD11"  …  "TAFD65"]
    OpticSim.GlassCat.NIKON => ["BAF10", "BAF11"  …  "_7054"]
    OpticSim.GlassCat.OHARA => ["L_BAL35", "L_BAL35P"  …  "S_TIM8"]
   OpticSim.GlassCat.SCHOTT => ["AF32ECO", "BAFN6"  …  "SFL6"]
   OpticSim.GlassCat.SUMITA => ["BAF1", "BAF10"  …  "ZNSF8"]
```
"""
glassnames() = [m => glassnames(m) for m in _child_modules(GlassCat)]

"""
    findglass(condition::Function) -> Vector{Glass}

Returns the list of glasses which satisfy `condition` where `condition::(Glass -> Bool)`.

# Example
```julia-repl
julia> findglass(x -> (x.Nd > 2.3 && x.λmin < 0.5 && x.λmax > 0.9))
8-element Array{GlassCat.Glass,1}:
 BIREFRINGENT.TEO2_E
 BIREFRINGENT.PBMOO4
 BIREFRINGENT.LINBO3
 INFRARED.CLEARTRAN_OLD
 INFRARED.CLEARTRAN
 INFRARED.SRTIO3
 INFRARED.ZNS_BROAD
 INFRARED.ZNS_VIS
```
"""
function findglass(condition::Function)
    # TODO - make the condition easier to specify (accessor functions for fields?)
    out = Vector{Glass}(undef, 0)
    for g in AGF_GLASSES
        if condition(g)
            push!(out, g)
        end
    end
    for g in OTHER_GLASSES
        if condition(g)
            push!(out, g)
        end
    end
    return out
end

function _child_modules(m::Module)
    ns = names(m, imported = false, all = true)
    ms = []
    for n in ns
        if n != nameof(m)
            try
                x = Core.eval(m, n)
                if x isa Module
                    push!(ms, x)
                end
            catch
            end
        end
    end
    return ms
end
