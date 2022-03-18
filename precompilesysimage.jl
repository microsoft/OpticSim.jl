#julia --project=. --trace-compile=file.jl 
using PackageCompiler
import Pkg

if Sys.iswindows()
    extension = ".dll"
elseif Sys.islinux()
    extension = ".so"
else
    throw(Error("Operating system was not one of: windows linux"))
end

info = Pkg.TOML.parsefile("Project.toml")

create_sysimage(
    collect(keys(info["deps"])),
    sysimage_path = "JuliaSysimage" * extension, 
    precompile_execution_file = "src/precompilation_functions.jl" #the functions in the file will be used to create precompilation statements. All of these functions should then execute very quickly the first time.
    )
