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

Modifies `sources` in-place such that only verified sources remain.
"""
function build_source_dir(sources::AbstractVector{<:AbstractVector{<:AbstractString}}, source_dir::AbstractString)
    mkpath(SOURCE_DIR)

    missing_sources = []
    for (i, source) in enumerate(sources)
        name, sha256 = source[1:2]
        source_file = joinpath(@__DIR__, source_dir, "$(name).agf")
        verified = verify_source(source_file, sha256)
        if !verified && length(source) >= 3
            download_source(source_file, source[3:end]...)
            verify_source(source_file, sha256)
        end
        if !verified
            push!(missing_sources, i)
        end
    end

    deleteat!(sources, missing_sources)
end

"""
Verify a source file using SHA256, returning true if successful. Otherwise, remove the file and return false.

Note: this isn't a security measure - it's possible to add a file back in after verification.
"""
function verify_source(source_file::AbstractString, sha256::AbstractString)
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
function download_source(source_file::AbstractString, url::AbstractString, POST_data::Maybe{AbstractString} = nothing)
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

"""
Returns a nested dict where each kvp represents a module (e.g. SCHOTT, NIKON), loaded from verified sources.
"""
function load_verified_sources(verified_source_names::Vector{<:AbstractString}, source_dir::AbstractString)
    return Dict(
        name => parse_source_file(joinpath(source_dir, "$(name).agf"))
        for name in verified_source_names
    )
end

function generate_cat_jl(cat::Dict{<:AbstractString,<:Dict}, jlpath::AbstractString)
    catalogs = []
    io = open(jlpath, "w")
    idnum = 1
    glassnames = []
    for (catalog_name, catalog) in cat
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
                # skip docstrings for CI builds - this prevents a lot of missing docstring warnings in makedocs
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
                push!(eval_string, "const $glass_name = Glass($(join(vals, ", ")))\nexport $glass_name")
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

"""
Parse an AGF source file into a native Julia dictionary, where each kvp is a glass in the AGF catalog.
"""
function parse_source_file(source_file::AbstractString)
    catalog_dict = Dict{String,Any}()

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
            od = string_list_to_float_list(split(line)[2:end])
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
                entry = SVector{3,Float64}(it_row[1], it_row[2], it_row[3])
                push!(transmission_data, entry)
            end
            catalog_dict[glass_name]["transmission"] = transmission_data
        end
    end
    return catalog_dict
end
