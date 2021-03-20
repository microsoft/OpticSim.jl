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

const AGF_DIR = joinpath(@__DIR__, "downloads", "glasscat") # contains SCHOTT.agf, Sumita.agf, etc.
const GLASSCAT_DIR = joinpath(@__DIR__, "..", "src", "GlassCat") # contains GlassCat.jl (pre-existing)
const JL_DIR = joinpath(GLASSCAT_DIR, "data") # contains AGFGlasscat.jl, SCHOTT.jl, etc.

mkpath(AGF_DIR)
mkpath(JL_DIR)

const SOURCES_PATH = joinpath(@__DIR__, "sources.txt")
const AGFGLASSCAT_PATH = joinpath(JL_DIR, "AGFGlassCat.jl")

include(joinpath(GLASSCAT_DIR, "GlassTypes.jl"))
include("utils.jl")

# Build a source directory using information from sources.txt
sources = [split(line, " ") for line in readlines(SOURCES_PATH)]
verify_sources!(sources, AGF_DIR)

verified_source_names = [source[1] for source in sources]
@info "$(isfile(AGFGLASSCAT_PATH) ? "Re-g" : "G")enerating $AGFGLASSCAT_PATH"
@info "Using sources: $(join(verified_source_names, ", ", " and "))"
generate_jls(verified_source_names, AGFGLASSCAT_PATH, AGF_DIR, JL_DIR)
