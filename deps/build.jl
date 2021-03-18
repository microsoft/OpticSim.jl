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

const GLASSCAT_ROOT_DIR = joinpath(@__DIR__, "..", "src", "GlassCat")
const GLASSCAT_DATA_DIR = joinpath(GLASSCAT_ROOT_DIR, "data")
const SOURCE_DIR = joinpath("downloads", "glasscat")

const SOURCES_FILE = "sources.txt"
const GLASS_JL_FILE = "AGFGlassCat.jl"
const GLASS_JL_PATH = joinpath(GLASSCAT_DATA_DIR, GLASS_JL_FILE)

include(joinpath(GLASSCAT_ROOT_DIR, "GlassTypes.jl"))
include("utils.jl")

mkpath(SOURCE_DIR)

# Build a source list from SOURCES_FILE
sources = [split(line, " ") for line in readlines(SOURCES_FILE)]

buildsourcedir(sources, SOURCE_DIR)

@info "$(isfile(GLASS_JL_PATH) ? "Re-g" : "G")enerating $GLASS_JL_PATH\n"
generate_cat_jl(load_glass_db(SOURCE_DIR), GLASS_JL_PATH)
