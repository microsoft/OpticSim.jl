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

const GLASSCAT_ROOT = joinpath(@__DIR__, "..", "src", "GlassCat")
const GLASS_JL_FILE = "AGFGlassCat.jl"
const GLASS_JL_PATH = joinpath(GLASSCAT_ROOT, GLASS_JL_FILE)
const DOWNLOAD_DIR = mktempdir()

include(joinpath(GLASSCAT_ROOT, "GlassTypes.jl"))
include("utils.jl")

glasspaths = [joinpath(GLASSCAT_ROOT, "OtherGlassCat.jl")]

# load all AGF files to dictionary and
# for each catalog create a module with structs containing the default values
# this is done as a string then written to a file which is then included
# download glassfiles that are easily accessed from the web so user will have some glasses to start with
downloadcatalogs(DOWNLOAD_DIR)

@info "Generating $GLASS_JL_PATH\n"
generate_cat_jl(load_glass_db(DOWNLOAD_DIR), GLASS_JL_PATH)
push!(glasspaths, GLASS_JL_PATH)

# write a deps/deps.jl file that points to the generated glass cat source files
open("deps.jl", "w") do io
    for glasspath in glasspaths
        line = "include(raw\"$glasspath\")"
        @info "$line > deps.jl"
        println(io, line)
    end
end
