using Documenter
using Optics

makedocs(
    sitename = "Optics.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        assets = [asset("assets/logo.svg", class = :ico, islocal = true)],
    ),
    modules = [Optics],
    pages = [
        "Home" => "index.md",
        "Glass Functions" => "glasscat.md",
        # "Glasses" => cat_pages,
        "Geometry" => [
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
        "Examples" => "examples.md",
        "Optimization" => "optimization.md",
        "Reference" => "ref.md",
        "Roadmap" => "roadmap.md"
    ],
    expandfirst = ["systems.md", "vis.md"])

deploydocs(
    repo = "github.com/microsoft/OpticsSimulation.git",
    devbranch = "main",
)

# function children(m::Module)
#     ns = names(m, imported = false, all = true)
#     ms = []
 
#     for n in ns
#         try
#             x = Core.eval(m, n)
#             if x isa Module
#                 if(x != Optics.GlassCat)
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
# catalogs = children(Optics.GlassCat)
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
#                 Optics.GlassCat.docstring(io, glass)
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
# Optics.GlassCat.include(joinpath(@__DIR__, "../src/AGFGlassCatDocs.jl"))

# clean up
# rm(joinpath(@__DIR__, "../src/AGFGlassCatDocs.jl"))
# for p in cat_pages
#     rm(joinpath(@__DIR__, "src/" * p[2]))
# end
