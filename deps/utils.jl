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

import HTTP
import SHA
import ZipFile

using StringEncodings
using Unitful
using StaticArrays
import Unitful: Length, Temperature, Quantity, Units

const Maybe{T} = Union{T, Nothing}

"""
Build a verified source directory of AGF files according to a specification given by `sources`.

Each `source ∈ sources` is a collection of strings in the format `name, sha256, url [, POST_body]`, where the last
optional string is used to specify data to be sent in a POST request. This allows us to download a greater range of
sources (e.g. SUMITA).
"""
function buildsourcedir(sources::AbstractVector{<:AbstractVector{<:AbstractString}}, source_dir::AbstractString)
    mkpath(SOURCE_DIR)
    for source in sources
        name, sha256 = source[1:2]
        source_file = joinpath(@__DIR__, source_dir, "$(name).agf")
        if !verifysource(source_file, sha256) && length(source) >= 3
            downloadsource(source_file, source[3:end]...)
            verifysource(source_file, sha256)
        end
    end
end

"""
Verify a source file using SHA256, returning true if successful. Otherwise, remove the file and return false.

Note: this isn't a security measure - it's possible to add a file back in after verification.
"""
function verifysource(source_file::AbstractString, sha256::AbstractString)
    if !isfile(source_file)
        @info "[-] Missing file at $source_file"
    elseif sha256 == SHA.bytes2hex(SHA.sha256(read(source_file)))
        @info "[✓] Verified file at $source_file"
        return true
    else
        @info "[x] Removing unverified file at $source_file"
        rm(source_file)
    end
    return false
end

"""
Download and unzip an AGF glass catalog from a publicly available source. Supports POST requests for interactive downloads.
"""
function downloadsource(source_file::AbstractString, url::AbstractString, POST_data::Maybe{AbstractString} = nothing)
    @info "Downloading source file from $url"
    try
        resp = isnothing(POST_data) ? HTTP.get(url) : HTTP.post(url, ["Content-Type" => "application/x-www-form-urlencoded"], POST_data)
        reader = ZipFile.Reader(IOBuffer(resp.body))
        write(source_file, read(reader.files[1]))
    catch e
        @error e
    end
end

########

function load_glass_db(directory::String)
    glass_catalog = Dict{String,Dict}()
    for path in readdir(directory)
        if uppercase(splitext(path)[2]) == ".AGF"
            parse_glass_file!(glass_catalog, joinpath(directory, path))
        end
    end
    return glass_catalog
end

function generate_cat_jl(cat, jlpath)
    catalogs = []
    io = open(jlpath, "w")
    idnum = 1
    glassnames = []
    for nameandcatalog in cat
        catalog_name, catalog = nameandcatalog
        push!(catalogs, catalog_name)
        eval_string = ["module $catalog_name"]
        push!(eval_string, "using ..GlassCat: Glass, GlassID, AGF")
        push!(eval_string, "using StaticArrays: SVector")
        for x in catalog
            let N = 0
                glass_name, glass_info = x
                push!(glassnames, "$catalog_name.$glass_name")
                vals = []
                for fn in string.(fieldnames(Glass))
                    if fn == "ID"
                        push!(vals, "GlassID(AGF, $idnum)")
                    elseif fn in ["D₀", "D₁", "D₂", "E₀", "E₁", "λₜₖ"]
                        push!(vals, repr(get(glass_info, fn, 0.0)))
                    elseif fn == "temp"
                        push!(vals, repr(get(glass_info, fn, 20.0)))
                    elseif fn == "transmission"
                        v = get(glass_info, "transmission", nothing)
                        if isnothing(v)
                            push!(vals, repr(nothing))
                        else
                            str = join(["($(join(a, ", ")))" for a in v], ", ")
                            push!(vals, "[$str]")
                        end
                    elseif fn == "transmissionN"
                        continue
                    else
                        push!(vals, repr(get(glass_info, fn, NaN)))
                    end
                end
                # skip docstrings for CI builds - this prevents missing docstring warnings in makedocs
                if isnothing(get(ENV, "CI", nothing))
                    raw_name = glass_info["raw_name"] == glass_name ? "" : " ($(glass_info["raw_name"]))"
                    doc_string = "\"\"\"    $catalog_name.$glass_name$raw_name\n"
                    doc_string *= "```\n$(rpad("ID:", 25))AGF:$idnum\n"
                    doc_string *= "$(rpad("RI @ 587nm:", 25))$(get(glass_info, "Nd", 0.0))\n"
                    doc_string *= "$(rpad("Abbe Number:", 25))$(get(glass_info, "Vd", 0.0))\n"
                    doc_string *= "$(rpad("ΔPgF:", 25))$(get(glass_info, "ΔPgF", 0.0))\n"
                    doc_string *= "$(rpad("TCE (÷1e-6):", 25))$(get(glass_info, "TCE", 0.0))\n"
                    doc_string *= "$(rpad("Density:", 25))$(get(glass_info, "p", 0.0))g/m³\n"
                    doc_string *= "$(rpad("Valid wavelengths:", 25))$(get(glass_info, "λmin", 0.0))μm to $(get(glass_info, "λmax", 0.0))μm\n"
                    doc_string *= "$(rpad("Reference Temp:", 25))$(get(glass_info, "temp", 20.0))°C\n"
                    doc_string *= "```\n\"\"\""
                    push!(eval_string, doc_string)
                end
                push!(eval_string, "const $glass_name = Glass($(join(vals, ", "))) \n export $glass_name")
            end
            idnum += 1
        end
        push!(eval_string, "end #module \n export $catalog_name \n") # module
        eval_string = join(eval_string, "\n")
        write(io, eval_string * "\n")
    end
    write(io, "const AGF_GLASS_NAMES = [$(join(repr.(glassnames), ", "))]\n")
    write(io, "const AGF_GLASSES = [$(join(glassnames, ", "))]\n")
    close(io)
end

function string_list_to_float_list(x)
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

function parse_glass_file!(glass_catalog, filename::String)
    if !isfile(filename)
        throw(error("AGF file doesn't exist"))
    end
    catalog_name = splitext(basename(filename))[1]
    # remove invalid characters
    catalog_name = replace(catalog_name, r"""[ ,.:;?!()&-]""" => "_")
    try
        # cant have module names which are just numbers so add a _ to the start
        parse(Int, catalog_name[1])
        catalog_name = "_" * catalog_name
    catch
        ()
    end
    glass_catalog[catalog_name] = Dict{String,Any}()
    # check whether the file is UTF8 or UTF16 encoded
    if !isvalid(readuntil(filename, " "))
        fo = open(filename, enc"UTF-16LE", "r")
    else
        fo = open(filename, "r")
    end
    # read the file
    glass_name = ""
    let transmission_data = nothing
        for line in readlines(fo)
            if strip(line) == "" || length(strip(line)) == 0 || startswith(line, "CC ") || startswith(line, "GC ")
                continue
            end
            if startswith(line, "NM ")
                transmission_data = Vector{SVector{3,Float64}}(undef, 0)
                nm = split(line)
                glass_name = nm[2]
                original_glass_name = glass_name
                # remove invalid characters
                glass_name = replace(glass_name, "*" => "_STAR")
                glass_name = replace(glass_name, r"""[ ,.:;?!()&-]""" => "_")
                try
                    # cant have module names which are just numbers so add a _ to the start
                    parse(Int, glass_name[1])
                    glass_name = "_" * glass_name
                catch
                    ()
                end
                glass_catalog[catalog_name][glass_name] = Dict{String,Any}()
                glass_catalog[catalog_name][glass_name]["raw_name"] = original_glass_name
                glass_catalog[catalog_name][glass_name]["dispform"] = Int(parse(Float64, nm[3]))
                glass_catalog[catalog_name][glass_name]["Nd"] = parse(Float64, nm[5])
                glass_catalog[catalog_name][glass_name]["Vd"] = parse(Float64, nm[6])
                if length(nm) < 7
                    glass_catalog[catalog_name][glass_name]["exclude_sub"] = 0
                else
                    glass_catalog[catalog_name][glass_name]["exclude_sub"] = Int(parse(Float64, nm[7]))
                end
                if length(nm) < 8
                    glass_catalog[catalog_name][glass_name]["status"] = 0
                else
                    glass_catalog[catalog_name][glass_name]["status"] = Int(parse(Float64, nm[8]))
                end
                if length(nm) < 9 || "-" ∈ nm
                    glass_catalog[catalog_name][glass_name]["meltfreq"] = 0
                else
                    glass_catalog[catalog_name][glass_name]["meltfreq"] = Int(parse(Float64, nm[9]))
                end
            elseif startswith(line, "ED ")
                ed = split(line)
                glass_catalog[catalog_name][glass_name]["TCE"] = parse(Float64, ed[2])
                glass_catalog[catalog_name][glass_name]["p"] = parse(Float64, ed[4])
                glass_catalog[catalog_name][glass_name]["ΔPgF"] = parse(Float64, ed[5])
                if (length(ed) < 6)
                    glass_catalog[catalog_name][glass_name]["ignore_thermal_exp"] = 0
                else
                    glass_catalog[catalog_name][glass_name]["ignore_thermal_exp"] = Int(parse(Float64, ed[6]))
                end
            elseif startswith(line, "CD ")
                cd = parse.(Float64, split(line)[2:end])
                glass_catalog[catalog_name][glass_name]["C1"] = get(cd, 1, NaN)
                glass_catalog[catalog_name][glass_name]["C2"] = get(cd, 2, NaN)
                glass_catalog[catalog_name][glass_name]["C3"] = get(cd, 3, NaN)
                glass_catalog[catalog_name][glass_name]["C4"] = get(cd, 4, NaN)
                glass_catalog[catalog_name][glass_name]["C5"] = get(cd, 5, NaN)
                glass_catalog[catalog_name][glass_name]["C6"] = get(cd, 6, NaN)
                glass_catalog[catalog_name][glass_name]["C7"] = get(cd, 7, NaN)
                glass_catalog[catalog_name][glass_name]["C8"] = get(cd, 8, NaN)
                glass_catalog[catalog_name][glass_name]["C9"] = get(cd, 9, NaN)
                glass_catalog[catalog_name][glass_name]["C10"] = get(cd, 10, NaN)
            elseif startswith(line, "TD ")
                td = parse.(Float64, split(line)[2:end])
                glass_catalog[catalog_name][glass_name]["D₀"] = get(td, 1, 0.0)
                glass_catalog[catalog_name][glass_name]["D₁"] = get(td, 2, 0.0)
                glass_catalog[catalog_name][glass_name]["D₂"] = get(td, 3, 0.0)
                glass_catalog[catalog_name][glass_name]["E₀"] = get(td, 4, 0.0)
                glass_catalog[catalog_name][glass_name]["E₁"] = get(td, 5, 0.0)
                glass_catalog[catalog_name][glass_name]["λₜₖ"] = get(td, 6, 0.0)
                glass_catalog[catalog_name][glass_name]["temp"] = get(td, 7, 20.0)
            elseif startswith(line, "OD ")
                od = string_list_to_float_list(split(line)[2:end])
                glass_catalog[catalog_name][glass_name]["relcost"] = get(od, 1, -1)
                glass_catalog[catalog_name][glass_name]["CR"] = get(od, 2, -1)
                glass_catalog[catalog_name][glass_name]["FR"] = get(od, 3, -1)
                glass_catalog[catalog_name][glass_name]["SR"] = get(od, 4, -1)
                glass_catalog[catalog_name][glass_name]["AR"] = get(od, 5, -1)
                glass_catalog[catalog_name][glass_name]["PR"] = get(od, 6, -1)
            elseif startswith(line, "LD ")
                ld = parse.(Float64, split(line)[2:end])
                glass_catalog[catalog_name][glass_name]["λmin"] = ld[1]
                glass_catalog[catalog_name][glass_name]["λmax"] = ld[2]
            elseif startswith(line, "IT ")
                it_row = parse.(Float64, split(line)[2:end])
                if length(it_row) == 3 && it_row[1] != 0.0
                    entry = SVector{3,Float64}(it_row[1], it_row[2], it_row[3])
                    push!(transmission_data, entry)
                end
                glass_catalog[catalog_name][glass_name]["transmission"] = transmission_data
            end
        end
    end
end
