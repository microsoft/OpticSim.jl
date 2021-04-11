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

module Cloud

using PyCall
using Conda
using Pkg
using Random

"""
    cache_run_config(subscription_id::String, resource_group::String, workspace_name::String, compute_name::String[, path::String])

Writes the AML config information to a file at `path`. If `path` isn't set then the config will be used globally for that OpticSim install.
"""
function cache_run_config(subscription_id::String, resource_group::String, workspace_name::String, compute_name::String, path::String = joinpath(@__DIR__, "amlconf"))
    open(path, "w") do io
        write(io, subscription_id * "\n")
        write(io, resource_group * "\n")
        write(io, workspace_name * "\n")
        write(io, compute_name)
   end
   nothing
end

"""
    get_cached_run_config([path::String])

Reads the AML config information from a file at `path`. If not specified then the global config will be read.
"""
function get_cached_run_config(path::String = joinpath(@__DIR__, "amlconf"))
    open(path, "r") do io
        subscription_id = readline(io)
        resource_group = readline(io)
        workspace_name = readline(io)
        compute_name = readline(io)
        return subscription_id, resource_group, workspace_name, compute_name
    end
end

"""
    submit_run_to_AML(run_name::String, path_to_script::String, script_args::Vector{String} = nothing, sampled_args:Dict{String,Vector{String}} = nothing, config_path::String; hyperdrive_concurrent_runs::Int = 10)
    submit_run_to_AML(run_name::String, path_to_script::String, subscription_id::String, resource_group::String, workspace_name::String, compute_name::String, script_args::Vector{String} = nothing, sampled_args::Dict{String, Vector{String}} = nothing; hyperdrive_concurrent_runs::Int = 10)

Submit a run to AML, `path_to_script` is relative to your local package root (i.e. location of `Project.toml`).
`script_args` are a series of arguments to your script as strings.
`sampled_args` is a dictionary where keys are argument names and values are lists of values (as strings) that that argument will take.
`config_path` is a path to a config file as written by [`cache_run_config`](@ref), if not specified the global config is used. Alternatively this information can be provided directly using the second method above.
`hyperdrive_concurrent_runs` is the maximum number of concurrent runs that will execute on AML (limited by your compute cluster size).
"""
function submit_run_to_AML(run_name::String, path_to_script::String, script_args::Union{Nothing,Vector{String}} = nothing,
                           sampled_args::Union{Nothing,Dict{String,Vector{String}}} = nothing;
                           config_path::String = joinpath(@__DIR__, "amlconf"),
                           hyperdrive_concurrent_runs::Int = 10)
    subscription_id, resource_group, workspace_name, compute_name = get_cached_run_config(config_path)
    submit_run_to_AML(run_name, path_to_script, subscription_id, resource_group, workspace_name, compute_name,
                      script_args, sampled_args, hyperdrive_concurrent_runs=hyperdrive_concurrent_runs)
end

function submit_run_to_AML(run_name::String, path_to_script::String,
                           subscription_id::String, resource_group::String, workspace_name::String, compute_name::String,
                           script_args::Union{Nothing,Vector{String}} = nothing,
                           sampled_args::Union{Nothing,Dict{String,Vector{String}}} = nothing;
                           hyperdrive_concurrent_runs::Int = 10)

    dockerfile = open(joinpath(@__DIR__, "dockerfile")) do file
        read(file, String)
    end

    # add deps to dockerfile
    project_dict = Pkg.TOML.parsefile(Base.active_project())
    packages = collect(keys(project_dict["deps"]))
    sort!(packages)
    pkg_install_cmd = "RUN julia -e \"using Pkg; "
    for package_name in packages
        if "compat" in keys(project_dict) && package_name in keys(project_dict["compat"])
            pkg_install_cmd = pkg_install_cmd * "Pkg.add(name=\\\"" * package_name * "\\\", version=\\\"" * project_dict["compat"][package_name] * "\\\");"
        else
            pkg_install_cmd = pkg_install_cmd * "Pkg.add(\\\"" * package_name * "\\\");"
        end
    end
    if "OpticSim" in packages
        # build OpticSim if it is there
        pkg_install_cmd = pkg_install_cmd * "Pkg.build(\\\"OpticSim\\\");"
    end
    pkg_install_cmd = pkg_install_cmd * "\""

    dockerfile = dockerfile * pkg_install_cmd

    # TODO maybe compile sysimage in docker - would be horribly slow but should speed up import a lot?

    source_directory = joinpath(dirname(Base.active_project()))

    if isfile(joinpath(source_directory, "Manifest.toml")) && !isfile(joinpath(source_directory, ".amlignore"))
        println("No .amlignore file found, creating one")
        open(joinpath(source_directory, ".amlignore"), "w") do io
            write(io, "Manifest.toml\n")
        end
    end

    # set up env for python stuff
    try
        pyimport("azureml.core")
    catch
        # FIXME maybe won't work, might need to restart Julia after this?
        Conda.add("python=3.7")
        Conda.add("pip=20.1.1")
        Conda.pip_interop(true)
        Conda.pip("install", "azureml-sdk")
        Pkg.build("PyCall")
    end

    # copy entry_script from here to source_directory
    entry_script_path = "entry_script_" * randstring() * ".py"
    cp(joinpath(@__DIR__, "entry_script.py"), joinpath(source_directory, entry_script_path))

    py"""
    import os
    import webbrowser
    from azureml.core import Environment, Experiment, Run, Workspace, ScriptRunConfig
    import azureml.train.hyperdrive as hyperdrive
    from azureml.train.hyperdrive.parameter_expressions import choice

    def get_hyperparam_dict(param_dict):
        hyper_param_dict = {}
        num_params = 1
        for key, value in param_dict.items():
            hyper_param_dict[key] = choice(value)
            num_params = num_params * len(value)

        return hyper_param_dict, num_params

    def submit_run(subscription_id, resource_group, workspace_name, compute_name, source_directory, julia_script,
                   script_args, sampled_args, run_name, dockerfile, entry_script_path, hyperdrive_concurrent_runs):
        workspace = Workspace(subscription_id, resource_group, workspace_name)

        compute_target = workspace.compute_targets[compute_name]

        env = Environment("opticsim")
        env.docker.base_image = None
        env.docker.base_dockerfile = dockerfile

        args = [julia_script.replace(os.sep, "/")]
        if script_args is not None:
            args.extend(script_args)

        src = ScriptRunConfig(source_directory=source_directory,
                              script=entry_script_path,
                              arguments=args,
                              compute_target=compute_target,
                              environment=env)
        src.run_config.docker.use_docker = True

        exp_name = os.getlogin() + "-opticsim"
        experiment = Experiment(workspace, exp_name)

        if sampled_args is not None:
            sampling_params, num_params = get_hyperparam_dict(sampled_args)
            param_sampling = hyperdrive.GridParameterSampling(sampling_params)
            hyperdrive_run_config = hyperdrive.HyperDriveConfig(run_config=src,
                                                                hyperparameter_sampling=param_sampling,
                                                                max_concurrent_runs=hyperdrive_concurrent_runs,
                                                                primary_metric_name="",
                                                                primary_metric_goal=hyperdrive.PrimaryMetricGoal.MINIMIZE,
                                                                max_total_runs=num_params)
            run_object = experiment.submit(hyperdrive_run_config, tags={"run_name": run_name})
        else:
            run_object = experiment.submit(src, tags={"run_name": run_name})

        webbrowser.open_new(run_object.get_portal_url())
    """

    py"submit_run"(subscription_id, resource_group, workspace_name, compute_name, source_directory, path_to_script,
                   script_args, sampled_args, run_name, dockerfile, entry_script_path, hyperdrive_concurrent_runs)

    # remove the entry script
    rm(joinpath(source_directory, entry_script_path))
end

export submit_run_to_AML, cache_run_config, get_cached_run_config

end