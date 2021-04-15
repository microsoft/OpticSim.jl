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

include(joinpath("..", "src", "GlassCat", "constants.jl"))
include(joinpath(GLASSCAT_DIR, "GlassTypes.jl"))
include(joinpath(GLASSCAT_DIR, "sources.jl"))
include(joinpath(GLASSCAT_DIR, "generate.jl"))

mkpath(AGF_DIR)
mkpath(JL_DIR)

# Build/verify a source directory using information from sources.txt
sources = split.(readlines(SOURCES_PATH))
verify_sources!(sources, AGF_DIR)
verified_source_names = first.(sources)

# Use verified sources to generate required .jl files
@info "Using sources: $(join(verified_source_names, ", ", " and "))"
generate_jls(verified_source_names, AGFGLASSCAT_NAME, JL_DIR, AGF_DIR)
