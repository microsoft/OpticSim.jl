# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

using Documenter
using OpticSim

makedocs(
    sitename = "OpticSim.jl",
    format = Documenter.HTML(
        # prettyurls = get(ENV, "CI", nothing) == "true",
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
            "Emitters (NEW)" => "emitters_new.md",
            "Interfaces" => "interfaces.md",
            "Lenses" => "lenses.md"
        ],
        "Visualization" => "vis.md",
        "Glass Functions" => "glasscat.md",
        "Optimization" => "optimization.md",
        "Cloud Execution" => "cloud.md",
        "Notebook utilities" => "notebooksutils.md",
        "Reference" => "ref.md",
        "Roadmap" => "roadmap.md"
    ],
    expandfirst = ["glasscat.md", "systems.md", "vis.md"])

deploydocs(
    repo = "github.com/microsoft/OpticSim.jl.git",
    devbranch = "main",
    push_preview = true,
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
