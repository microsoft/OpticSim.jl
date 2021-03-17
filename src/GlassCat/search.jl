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

"""
    glasscatalogs()

Prints the complete list of glass catalogs available from GlassCat.

## Example
```julia-repl
julia> glasscatalogs()
41 glass catalogs:
GlassCat.AMTIR
GlassCat.ANGSTROMLINK
GlassCat.APEL
GlassCat.ARCHER
GlassCat.ARTON
GlassCat.AUER_LIGHTING
GlassCat.BIREFRINGENT
...
```
"""
function glasscatalogs()
    println("$(length(_child_modules(GlassCat))) glass catalogs:")
    for c in _child_modules(GlassCat)
        println("$c")
    end
end
export glasscatalogs

"""
    glasses(catalog::Union{Module,Nothing} = nothing)

Prints the glasses available from a given catalog, or if none if specified from all catalogs.

# Example
```julia-repl
julia> glasses(GlassCat.CARGILLE)
GlassCat.CARGILLE - 3 glasses:
	OG0607
	OG0608
	OG081160
```
"""
function glasses(catalog::Union{Module,Nothing} = nothing)
    function printforcat(cat::Module)
        glass_names = names(cat, all = true, imported = false)
        glasses = []
        for glass_name in glass_names
            glass_name_str = string(glass_name)
            if !occursin("#", glass_name_str) && glass_name_str != "eval" && glass_name_str != "include" && glass_name != nameof(cat)
                push!(glasses, "$glass_name_str")
            end
        end
        println("$cat - $(length(glasses)) glasses:")
        for g in glasses
            println("\t$g")
        end
    end
    if catalog === nothing
        for c in _child_modules(GlassCat)
            printforcat(c)
        end
    else
        printforcat(catalog)
    end
end
export glasses

"""
    findglass(condition::Function) -> Vector{Glass}

Returns the list of glasses which satisfy `condition` where `condition::(Glass -> Bool)`.

TODO - make the condition easier to specify (accessor functions for fields?)

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
