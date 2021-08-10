# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

import HTTP
import SHA
import ZipFile
using Pkg

function match_agf_basename(s::AbstractString)
    m = match(r"^([a-zA-Z]+)\.(agf|AGF)$", basename(s))
    return m === nothing ? nothing : m[1]
end

is_alphabetical(s::AbstractString) = match(r"^([a-zA-Z]+)$", s) !== nothing

is_duplicate(name::AbstractString, sourcefile::AbstractString) = name ∈ first.(split.(readlines(sourcefile)))

"""
    add_agf(agffile; agfdir = AGF_DIR, sourcefile = SOURCES_PATH, name = nothing, rebuild = true)

Copies a file at `agffile` (this can be either a download link or local path) to `agfdir` and appends a corresponding
entry to the source list at `sourcefile`.

If a `name` is not provided for the catalog, an implicit name is derived from `agffile`.

If `rebuild` is true, Pkg.build is called at the end to install the new catalog.
"""
function add_agf(
    agffile::AbstractString;
    agfdir::AbstractString = AGF_DIR,
    sourcefile::AbstractString = SOURCES_PATH,
    name::Union{Nothing, AbstractString} = nothing,
    rebuild::Bool = true
)
    # check name
    if name === nothing
        name = match_agf_basename(agffile)
        if name === nothing
            @error "invalid implicit catalog name \"$(basename(agffile))\". Should be purely alphabetical with a .agf/.AGF extension."
            return
        end
        name = uppercase(name)
    end

    if !is_alphabetical(name)
        @error "invalid catalog name \"$name\". Should be purely alphabetical."
        return
    end

    if is_duplicate(name, sourcefile)
        @error "adding the catalog name \"$name\" would create a duplicate entry in source file $sourcefile"
        return
    end

    # download agffile from url if necessary
    mkpath(agfdir)
    dest = joinpath(agfdir, name * ".agf")

    if isfile(agffile)
        cp(agffile, dest, force=true)
    else
        download_source(dest, agffile)
        if !isfile(dest)
            @error "failed to download from $agffile"
            return
        end
    end

    # append a corresponding entry to the source list at sourcefile
    sha256sum = SHA.bytes2hex(SHA.sha256(read(dest)))
    open(sourcefile, "a") do io
        source = isfile(agffile) ? [name, sha256sum] : [name, sha256sum, agffile]
        write(io, join(source, ' ') * '\n')
    end

    # optional rebuild
    if rebuild
        @info "Re-building OpticSim.jl"
        Pkg.build("OpticSim"; verbose=true)
    end
end

"""
    verify_sources!(sources::AbstractVector{<:AbstractVector{<:AbstractString}}, agfdir::AbstractString)

Verify a list of `sources` located in `agfdir`. If AGF files are missing or invalid, try to download them using the
information provided in `sources`.

Each `source ∈ sources` is a collection of strings in the format `name, sha256sum, url [, POST_data]`, where the last
optional string is used to specify data to be sent in a POST request. This allows us to download a greater range of
sources (e.g. SUMITA).

Modifies `sources` in-place such that only verified sources remain.
"""
function verify_sources!(sources::AbstractVector{<:AbstractVector{<:AbstractString}}, agfdir::AbstractString)
    # track missing sources as we go and delete them afterwards to avoid modifying our iterator
    missing_sources = []

    for (i, source) in enumerate(sources)
        name, sha256sum = source[1:2]
        agffile = joinpath(agfdir, "$(name).agf")
        verified = verify_source(agffile, sha256sum)
        if !verified && length(source) >= 3
            # try downloading and re-verifying the source if download information is provided (sources[3:end])
            download_source(agffile, source[3:end]...)
            verified = verify_source(agffile, sha256sum)
        end
        if !verified
            push!(missing_sources, i)
        end
    end

    deleteat!(sources, missing_sources)
end

"""
    verify_source(agffile::AbstractString, expected_sha256sum::AbstractString)

Verify a source file using SHA256, returning true if successful. Otherwise, remove the file and return false.
"""
function verify_source(agffile::AbstractString, expected_sha256sum::AbstractString)
    if !isfile(agffile)
        @info "[-] Missing file at $agffile"
    else
        sha256sum = SHA.bytes2hex(SHA.sha256(read(agffile)))
        if expected_sha256sum == sha256sum
            @info "[✓] Verified file at $agffile"
            return true
        else
            @info "[x] Removing unverified file at $agffile (expected $expected_sha256sum, got $sha256sum)"
            rm(agffile)
        end
    end
    return false
end

"""
    download_source(dest::AbstractString, url::AbstractString, POST_data::Union{Nothing,AbstractString} = nothing)

Download and unzip an AGF glass catalog from a publicly available source. Supports POST requests.
"""
function download_source(dest::AbstractString, url::AbstractString, POST_data::Union{Nothing,AbstractString} = nothing)
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
        write(dest, agfdata)
    catch e
        @error e
    end
end
