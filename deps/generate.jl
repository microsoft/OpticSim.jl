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

using StringEncodings
using StaticArrays
using Unitful
import Unitful: Length, Temperature, Quantity, Units

"""
Generates .jl files: a `mainfile` and several catalog files.

Each catalog file is a module representing a distinct glass catalog (e.g. NIKON, SCHOTT), generated from corresponding
AGF files in `sourcedir`. These are then included and exported in `mainfile`.
"""
function generate_jls(
    sourcenames::Vector{<:AbstractString}, mainfile::AbstractString, jldir::AbstractString, sourcedir::AbstractString
)
    id = 1
    catalogfiles = []
    glassnames = []

    # generate several catalog files (.jl)
    for catalogname in sourcenames
        # parse the sourcefile (.agf) into a catalog (native Julia dictionary)
        sourcefile = joinpath(sourcedir, "$(catalogname).agf")
        catalog = sourcefile_to_catalog(sourcefile)

        # parse the catalog into a module string and write it to a catalog file (.jl)
        id, modstring = catalog_to_modstring(id, catalogname, catalog)
        push!(catalogfiles, "$(catalogname).jl")
        open(joinpath(jldir, catalogfiles[end]), "w") do io
            write(io, modstring)
        end

        # track glass names for lookup purposes
        append!(glassnames, ["$(catalogname).$(glassname)" for glassname in keys(catalog)])
    end

    # generate the parent main file (.jl)
    agfstrings = [
        "export $(join(sourcenames, ", "))",
        "",
        ["include(\"$(catalogfile)\")" for catalogfile in catalogfiles]...,
        "",
        "const AGF_GLASS_NAMES = [$(join(repr.(glassnames), ", "))]",
        "const AGF_GLASSES = [$(join(glassnames, ", "))]",
        ""
    ]
    open(joinpath(jldir, mainfile), "w") do io
        write(io, join(agfstrings, "\n"))
    end
end

"""
Parse a `sourcefile` (.agf) into a native dictionary, where each `kvp = (glassname, glassinfo)` is a glass.
"""
function sourcefile_to_catalog(source_file::AbstractString)
    catalog_dict = Dict{String,Dict{String}}()

    # check whether the file is UTF8 or UTF16 encoded
    is_utf8 = isvalid(readuntil(source_file, " "))
    fo = is_utf8 ? open(source_file, "r") : open(source_file, enc"UTF-16LE", "r")

    # store glass_name and transmission_data as persistent variables in between loops
    glass_name = ""
    transmission_data = nothing

    # read the file
    for line in readlines(fo)
        if strip(line) == "" || length(strip(line)) == 0 || startswith(line, "CC ") || startswith(line, "GC ")
            continue
        end
        if startswith(line, "NM ")
            nm = split(line)

            glass_name = make_valid_name(nm[2])
            transmission_data = Vector{SVector{3,Float64}}(undef, 0)

            catalog_dict[glass_name] = Dict{String,Any}()
            catalog_dict[glass_name]["raw_name"] = nm[2]
            catalog_dict[glass_name]["dispform"] = Int(parse(Float64, nm[3]))
            catalog_dict[glass_name]["Nd"] = parse(Float64, nm[5])
            catalog_dict[glass_name]["Vd"] = parse(Float64, nm[6])
            catalog_dict[glass_name]["exclude_sub"] = length(nm) < 7 ? 0 : Int(parse(Float64, nm[7]))
            catalog_dict[glass_name]["status"] = length(nm) < 8 ? 0 : Int(parse(Float64, nm[8]))
            catalog_dict[glass_name]["meltfreq"] = length(nm) < 9 || "-" ∈ nm ? 0 : Int(parse(Float64, nm[9]))
        elseif startswith(line, "ED ")
            ed = split(line)
            catalog_dict[glass_name]["TCE"] = parse(Float64, ed[2])
            catalog_dict[glass_name]["p"] = parse(Float64, ed[4])
            catalog_dict[glass_name]["ΔPgF"] = parse(Float64, ed[5])
            catalog_dict[glass_name]["ignore_thermal_exp"] = length(ed) < 6 ? 0 : Int(parse(Float64, ed[6]))
        elseif startswith(line, "CD ")
            cd = parse.(Float64, split(line)[2:end])
            catalog_dict[glass_name]["C1"] = get(cd, 1, NaN)
            catalog_dict[glass_name]["C2"] = get(cd, 2, NaN)
            catalog_dict[glass_name]["C3"] = get(cd, 3, NaN)
            catalog_dict[glass_name]["C4"] = get(cd, 4, NaN)
            catalog_dict[glass_name]["C5"] = get(cd, 5, NaN)
            catalog_dict[glass_name]["C6"] = get(cd, 6, NaN)
            catalog_dict[glass_name]["C7"] = get(cd, 7, NaN)
            catalog_dict[glass_name]["C8"] = get(cd, 8, NaN)
            catalog_dict[glass_name]["C9"] = get(cd, 9, NaN)
            catalog_dict[glass_name]["C10"] = get(cd, 10, NaN)
        elseif startswith(line, "TD ")
            td = parse.(Float64, split(line)[2:end])
            catalog_dict[glass_name]["D₀"] = get(td, 1, 0.0)
            catalog_dict[glass_name]["D₁"] = get(td, 2, 0.0)
            catalog_dict[glass_name]["D₂"] = get(td, 3, 0.0)
            catalog_dict[glass_name]["E₀"] = get(td, 4, 0.0)
            catalog_dict[glass_name]["E₁"] = get(td, 5, 0.0)
            catalog_dict[glass_name]["λₜₖ"] = get(td, 6, 0.0)
            catalog_dict[glass_name]["temp"] = get(td, 7, 20.0)
        elseif startswith(line, "OD ")
            od = stringlist_to_floatlist(split(line)[2:end])
            catalog_dict[glass_name]["relcost"] = get(od, 1, -1)
            catalog_dict[glass_name]["CR"] = get(od, 2, -1)
            catalog_dict[glass_name]["FR"] = get(od, 3, -1)
            catalog_dict[glass_name]["SR"] = get(od, 4, -1)
            catalog_dict[glass_name]["AR"] = get(od, 5, -1)
            catalog_dict[glass_name]["PR"] = get(od, 6, -1)
        elseif startswith(line, "LD ")
            ld = parse.(Float64, split(line)[2:end])
            catalog_dict[glass_name]["λmin"] = ld[1]
            catalog_dict[glass_name]["λmax"] = ld[2]
        elseif startswith(line, "IT ")
            it_row = parse.(Float64, split(line)[2:end])
            if length(it_row) == 3 && it_row[1] != 0.0
                entry = SVector(it_row[1], it_row[2], it_row[3])
                push!(transmission_data, entry)
            end
            catalog_dict[glass_name]["transmission"] = transmission_data
        end
    end
    return catalog_dict
end

"""
Make a valid Julia variable name from an arbitrary string.
"""
function make_valid_name(name::AbstractString)
    # remove invalid characters
    name = replace(name, "*" => "_STAR")
    name = replace(name, r"""[ ,.:;?!()&-]""" => "_")
    # cant have module names which are just numbers so add a _ to the start
    if tryparse(Int, "$(name[1])") !== nothing
        name = "_" * name
    end
    return name
end

function stringlist_to_floatlist(x::Vector{<:AbstractString})
    npts = length(x)
    if (npts == 0) || ((npts == 1) && (strip(x[1]) == "-"))
        return (repeat([-1.0], 10))
    end
    res = []
    for a in x
        if (strip(a) == "-")
            push!(res, -1.0)
        else
            try
                push!(res, parse(Float64, a))
            catch
                push!(res, NaN)
            end
        end
    end
    return (res)
end

"""
Convert a `catalog` dict into a `modstring` which can be written to a Julia source file.
"""
function catalog_to_modstring(start_id::Integer, catalogname::AbstractString, catalog::Dict{<:AbstractString})
    id = start_id
    isCI = haskey(ENV, "CI")

    modstrings = [
        "module $catalogname",
        "using ..GlassCat: Glass, GlassID, AGF",
        "export $(join(keys(catalog), ", "))",
        ""
    ]
    for (glassname, glassinfo) in catalog
        argstring = glassinfo_to_argstring(glassinfo, id)
        push!(modstrings, "const $glassname = Glass($argstring)")
        id += 1
    end
    append!(modstrings, ["end # module", ""]) # last "" is for \n at EOF

    return id, join(modstrings, "\n")
end

"""
Convert a `glassinfo` dict into a `docstring` to be prepended to a `Glass` const.
"""
function glassinfo_to_docstring(
    glassinfo::Dict{<:AbstractString}, id::Integer, catalogname::AbstractString, glassname::AbstractString
)
    raw_name = glassinfo["raw_name"] == glassname ? "" : " ($(glassinfo["raw_name"]))"
    pad(str, padding=25) = rpad(str, padding)
    getinfo(key, default=0.0) = get(glassinfo, key, default)

    return join([
        "\"\"\"    $catalogname.$glassname$raw_name",
        "```",
        "$(pad("ID:"))AGF:$id",
        "$(pad("RI @ 587nm:"))$(getinfo("Nd"))",
        "$(pad("Abbe Number:"))$(getinfo("Vd"))",
        "$(pad("ΔPgF:"))$(getinfo("ΔPgF"))",
        "$(pad("TCE (÷1e-6):"))$(getinfo("TCE"))",
        "$(pad("Density:"))$(getinfo("p"))g/m³",
        "$(pad("Valid wavelengths:"))$(getinfo("λmin"))μm to $(getinfo("λmax"))μm",
        "$(pad("Reference Temp:"))$(getinfo("temp", 20.0))°C",
        "```",
        "\"\"\""
    ], "\n")
end

"""
Convert a `glassinfo` dict into an `argstring` to be passed into a `Glass` constructor.
"""
function glassinfo_to_argstring(glassinfo::Dict{<:AbstractString}, id::Integer)
    argstrings = []
    for fn in string.(fieldnames(Glass))
        if fn == "ID"
            push!(argstrings, "GlassID(AGF, $id)")
        elseif fn in ["D₀", "D₁", "D₂", "E₀", "E₁", "λₜₖ"]
            push!(argstrings, repr(get(glassinfo, fn, 0.0)))
        elseif fn == "temp"
            push!(argstrings, repr(get(glassinfo, fn, 20.0)))
        elseif fn == "transmission"
            v = get(glassinfo, "transmission", nothing)
            if isnothing(v)
                push!(argstrings, repr(nothing))
            else
                str = join(["($(join(a, ", ")))" for a in v], ", ")
                push!(argstrings, "[$str]")
            end
        elseif fn == "transmissionN"
            continue
        else
            push!(argstrings, repr(get(glassinfo, fn, NaN)))
        end
    end
    return join(argstrings, ", ")
end
