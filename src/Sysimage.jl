module Sysimage
import Pkg, Libdl, PackageCompiler

function compile(sysimage_path = "JuliaSysimage.$(Libdl.dlext)")
    env_to_precompile = joinpath(@__DIR__, "..")
    precompile_execution_file = joinpath(env_to_precompile, "Precompile.jl")
    project_filename = joinpath(env_to_precompile, "Project.toml")
    project = Pkg.API.read_project(project_filename)
    used_packages = Symbol.(collect(keys(project.deps)))
    # don't need these ever after building this so no need to have them in the sysimage
    used_packages = filter(x -> x ∉ [:Libdl, :PackageCompiler, :Pkg], used_packages)
    if Libdl.dlext == "dll"
        @warn "Ignoring packages which use gl dlls on Windows as these cause build errors"
        # see https://github.com/JuliaLang/PackageCompiler.jl/issues/365
        # recently the build leaves a corrupt dll if we include these rather than the error in the issue above
        used_packages = filter(x -> x ∉ [:Makie, :ImageView], used_packages)
    end
    @info "Building a custom sysimage for OpticSim.jl."
    PackageCompiler.create_sysimage(used_packages, sysimage_path = sysimage_path, project = env_to_precompile, precompile_execution_file = precompile_execution_file)
end

end # module Sysimage
