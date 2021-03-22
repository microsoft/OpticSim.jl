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

import Pkg, Libdl, PackageCompiler

function compile(sysimage_path = "JuliaSysimage.$(Libdl.dlext)")
    env_to_precompile = joinpath(@__DIR__, "..")
    precompile_execution_file = joinpath(@__DIR__, "precompile.jl")
    project_filename = joinpath(env_to_precompile, "Project.toml")
    project = Pkg.API.read_project(project_filename)
    used_packages = Symbol.(collect(keys(project.deps)))
    # don't need these ever after building this so no need to have them in the sysimage
    filter!(x -> x ∉ [:Libdl, :PackageCompiler, :Pkg], used_packages)
    if Libdl.dlext == "dll"
        @warn "Ignoring packages which use gl dlls on Windows as these cause build errors"
        # see https://github.com/JuliaLang/PackageCompiler.jl/issues/365
        # recently the build leaves a corrupt dll if we include these rather than the error in the issue above
        used_packages = filter(x -> x ∉ [:Makie, :ImageView], used_packages)
    end
    @info "Building a custom sysimage for OpticSim.jl."
    PackageCompiler.create_sysimage(used_packages, sysimage_path = sysimage_path, project = env_to_precompile, precompile_execution_file = precompile_execution_file)
end
