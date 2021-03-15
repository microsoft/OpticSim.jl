const GLASSCAT_ROOT = joinpath(@__DIR__, "..", "src", "GlassCat")
const GLASS_JL_FILE = "AGFGlassCat.jl"
const GLASS_JL_PATH = GLASS_JL_FILE
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
        println(io, "include(raw\"$(joinpath(pwd(), glasspath))\")")
    end
end
