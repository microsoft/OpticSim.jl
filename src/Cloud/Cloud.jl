module Cloud

using PyCall
using Conda
using Pkg
using Random

function cache_run_config(subscription_id::String, resource_group::String, workspace_name::String, compute_name::String, path::String = joinpath(@__DIR__, "amlconf"))
    open(path, "w") do io
        write(io, subscription_id * "\n")
        write(io, resource_group * "\n")
        write(io, workspace_name * "\n")
        write(io, compute_name)
   end
   nothing
end

function get_cached_run_config(path::String = joinpath(@__DIR__, "amlconf"))
    open(path, "r") do io
        subscription_id = readline(io)
        resource_group = readline(io)
        workspace_name = readline(io)
        compute_name = readline(io)
        return subscription_id, resource_group, workspace_name, compute_name
    end
end

function submit_run_to_AML(run_name::String, path_to_script::String, script_args::Union{Nothing,Vector{String}} = nothing, config_path::String = joinpath(@__DIR__, "amlconf"))
    subscription_id, resource_group, workspace_name, compute_name = get_cached_run_config(config_path)
    submit_run_to_AML(run_name, path_to_script, script_args, subscription_id, resource_group, workspace_name, compute_name)
end

function submit_run_to_AML(run_name::String, path_to_script::String, script_args::Union{Nothing,Vector{String}},
                           subscription_id::String, resource_group::String, workspace_name::String, compute_name::String)

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

    if isfile(joinpath(source_directory, "Manifest.toml"))
        println("WARNING: Manifest.toml must be included in .amlignore file.")
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
    # import azureml.train.hyperdrive as hyperdrive
    # from azureml.train.hyperdrive.parameter_expressions import choice

    def submit_run(subscription_id, resource_group, workspace_name, compute_name, source_directory, julia_script, script_args, run_name, dockerfile, entry_script_path):
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

        # TODO support hyperdrive

        exp_name = os.getlogin() + "-opticsim"
        experiment = Experiment(workspace, exp_name)
        run_object = experiment.submit(src, tags={"run_name": run_name})

        webbrowser.open_new(run_object.get_portal_url())
    """

    py"submit_run"(subscription_id, resource_group, workspace_name, compute_name, source_directory, path_to_script, script_args, run_name, dockerfile, entry_script_path)

    # remove the entry script
    rm(joinpath(source_directory, entry_script_path))
end

export submit_run_to_AML, cache_run_config, get_cached_run_config

end