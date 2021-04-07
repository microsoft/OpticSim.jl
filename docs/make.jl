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

using Documenter
using OpticSim, OpticSim.Geometry

makedocs(
    sitename = "OpticSim.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        assets = [asset("assets/logo.svg", class = :ico, islocal = true)],
    ),
    modules = [OpticSim],
    pages = [
        "Home" => "index.md",
        "Examples" => "examples.md",
        "Geometry" => [
            "Basic Types" => "basic_types.md",
            "Primitives" => "primitives.md",
            "CSG" => "csg.md"
        ],
        "Optical" => [
            "Systems" => "systems.md",
            "Emitters" => "emitters.md",
            "Interfaces" => "interfaces.md",
            "Lenses" => "lenses.md"
        ],
        "Visualization" => "vis.md",
        "Glass Functions" => "glasscat.md",
        "Optimization" => "optimization.md",
        "Cloud Execution" => "cloud.md",
        "Reference" => "ref.md",
        "Roadmap" => "roadmap.md"
    ],
    expandfirst = ["glasscat.md", "systems.md", "vis.md"])

deploydocs(
    repo = "github.com/microsoft/OpticSim.jl.git",
    devbranch = "main",
)

# function children(m::Module)
#     ns = names(m, imported = false, all = true)
#     ms = []
#     for n in ns
#         try
#             x = Core.eval(m, n)
#             if x isa Module
#                 if(x != OpticSim.GlassCat)
#                 println("x $x")
#                 push!(ms, x)
#             end
#         end
#         catch
#         end
#     end
#     return ms
# end

#WARNING: I think this code creates a type for each glass name, which overwrites the definition in the src files where each glass name corresponds to an integer indexing into a glass table. Obviously nothing works after this. Needs major surgery to generate glass documentation at the same time as the documentation for everything else.

# # write a code file for the catalog with docstrings
# io = open(joinpath(@__DIR__, "../src/AGFGlassCatDocs.jl"), "w")
# catalogs = children(OpticSim.GlassCat)
# cat_pages = []
# for catname in catalogs
#     println(catalogs)
#     eval_string = ["module $(nameof(catname))"]
#     escape_catalog_name = replace(string(nameof(catname)), "_" => "\\_")
#     push!(cat_pages, escape_catalog_name => "$(nameof(catname)).md")
#     iomd = open(joinpath(@__DIR__, "src/$(nameof(catname)).md"), "w")
#     write(iomd, "# $escape_catalog_name\n\n")
#     write(iomd, "```@raw html\n")
#     write(iomd, "<style>article p {display:flex;justify-content:space-between;}</style>\n")
#     write(iomd, "```\n")
#     write(iomd, "```@docs\n")
#     catalog_module = eval(catname)
#     glass_names = names(catalog_module, all = true, imported = false)
#     for glass_name in glass_names
#         glass_name_str = string(glass_name)
#         if !occursin("#", glass_name_str) && glass_name_str != "eval" && glass_name_str != "include" && glass_name != nameof(catname) && !in(glass_name,(:MODEL, :MIL, :AGF, :OTHER, :AIR))
#             glass = Core.eval(catalog_module, glass_name)

#             temp = typeof(glass)
# println("type of glass $temp")

#             let io = IOBuffer()
#                 OpticSim.GlassCat.docstring(io, glass)
#                 infostr = String(take!(io))
#                 push!(eval_string, "\"\"\"\n$infostr\n\"\"\"")
#             end
#             push!(eval_string, "function $glass_name_str()\nend")
#             write(iomd, "$catname.$glass_name_str\n")
#         end
#     end
#     write(iomd, "```\n")
#     push!(eval_string, "end") # module
#     eval_string = join(eval_string, "\n")
#     write(io, eval_string * "\n")
#     close(iomd)
# end
# close(io)

# # include source with docstrings
# OpticSim.GlassCat.include(joinpath(@__DIR__, "../src/AGFGlassCatDocs.jl"))

# clean up
# rm(joinpath(@__DIR__, "../src/AGFGlassCatDocs.jl"))
# for p in cat_pages
#     rm(joinpath(@__DIR__, "src/" * p[2]))
# end
