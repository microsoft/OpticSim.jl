# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

import HTTP
import SHA
import ZipFile
using Pkg

"""
    add_agf(sourcefile::AbstractString; name::Union{Nothing,AbstractString} = nothing, rebuild::Bool = true)

Adds an already downloaded AGF file to the sourcelist at data/sources.txt, generating the SHA256 checksum automatically.

Optionally provide a `name` for the corresponding module, and `rebuild` AGFGlassCat.jl by default.
"""
function add_agf(sourcefile::AbstractString; name::Union{Nothing,AbstractString} = nothing, rebuild::Bool = true)
    if !isfile(sourcefile)
        @error "AGF file not found at $sourcefile"
        return
    end

    # infer catalog name from sourcefile basename (alphabetical only)
    if name === nothing
        name = uppercase(match(r"^([a-zA-Z]+)\.agf$"i, basename(sourcefile))[1])
    end

    # avoid duplicate catalog names
    if name ∈ first.(split.(readlines(SOURCES_PATH)))
        @error "Adding the catalog name \"$name\" would create a duplicate entry in sources.txt"
        return
    end

    # copy sourcefile to correct location
    mkpath(AGF_DIR)
    cp(sourcefile, joinpath(AGF_DIR, name * ".agf"), force=true)

    # append a corresponding entry to sources.txt
    sha256sum = SHA.bytes2hex(SHA.sha256(read(sourcefile)))
    open(SOURCES_PATH, "a") do io
        write(io, join([name, sha256sum], ' ') * '\n')
    end

    if rebuild
        @info "Re-building OpticSim.jl"
        Pkg.build("OpticSim"; verbose=true)
    end
end

"""
    verify_sources!(sources::AbstractVector{<:AbstractVector{<:AbstractString}}, sourcedir::AbstractString)

Verify a list of `sources` located in `sourcedir`. If AGF files are missing or invalid, try to download them using the
information provided in `sources`.

Each `source ∈ sources` is a collection of strings in the format `name, sha256sum, url [, POST_data]`, where the last
optional string is used to specify data to be sent in a POST request. This allows us to download a greater range of
sources (e.g. Sumita).

Modifies `sources` in-place such that only verified sources remain.
"""
function verify_sources!(sources::AbstractVector{<:AbstractVector{<:AbstractString}}, sourcedir::AbstractString)
    # track missing sources as we go and delete them afterwards to avoid modifying our iterator
    missing_sources = []

    for (i, source) in enumerate(sources)
        name, sha256sum = source[1:2]
        sourcefile = joinpath(sourcedir, "$(name).agf")
        verified = verify_source(sourcefile, sha256sum)
        if !verified && length(source) >= 3
            # try downloading and re-verifying the source if download information is provided (sources[3:end])
            download_source(sourcefile, source[3:end]...)
            verified = verify_source(sourcefile, sha256sum)
        end
        if !verified
            push!(missing_sources, i)
        end
    end

    deleteat!(sources, missing_sources)
end

"""
    verify_source(sourcefile::AbstractString, sha256sum::AbstractString)

Verify a source file using SHA256, returning true if successful. Otherwise, remove the file and return false.
"""
function verify_source(sourcefile::AbstractString, sha256sum::AbstractString)
    if !isfile(sourcefile)
        @info "[-] Missing file at $sourcefile"
    elseif sha256sum == SHA.bytes2hex(SHA.sha256(read(sourcefile)))
        @info "[✓] Verified file at $sourcefile"
        return true
    else
        @info "[x] Removing unverified file at $sourcefile"
        rm(sourcefile)
    end
    return false
end

"""
    download_source(sourcefile::AbstractString, url::AbstractString, POST_data::Union{Nothing,AbstractString} = nothing)

Download and unzip an AGF glass catalog from a publicly available source. Supports POST requests.
"""
function download_source(sourcefile::AbstractString, url::AbstractString, POST_data::Union{Nothing,AbstractString} = nothing)
    @info "Downloading source file from $url"
    try
        headers = ["Content-Type" => "application/x-www-form-urlencoded"]
        resp = POST_data === nothing ? HTTP.get(url) : HTTP.post(url, headers, POST_data)

        # save agf file, unzipping if necessary
        if endswith(url, ".agf")
            agfdata = resp.body
        else
            reader = ZipFile.Reader(IOBuffer(resp.body))
            agfdata = read(reader.files[findfirst(f -> endswith(lowercase(f.name), ".agf"), reader.files)])
        end
        write(sourcefile, agfdata)
    catch e
        @error e
    end
end
