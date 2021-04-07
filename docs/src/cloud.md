# Cloud Execution

## Azure

A key benefit and design motivation of OpticSim is being able to execute many simulations/optimizations at once.
This is best enabled through the use of cloud computing services such as Azure.

As part of the base package, we provide support for cloud execution using an [Azure Machine Learning](https://azure.microsoft.com/en-gb/free/machine-learning) workspace.

To use this functionally you'll first need to set up an AML workspace with a compute cluster. Then you'll need to provide a few bits of information to OpticSim:

- Subscription ID, of the form `XXXXXXXX-XXX-XXX-XXXX-XXXXXXXXXXXX`
- Resource group name
- Workspace name
- Compute cluster name

This information can be cached either to a specific file, or globally:

```@docs
Optics.Cloud.cache_run_config
Optics.Cloud.get_cached_run_config
```

You should also include an `.amlignore` file in the root of your project.
This is similar to a `.gitignore` file and should include any files which should not be uploaded to AML as part of your source snapshot, for examples `test/`.
**`Manifest.toml` must be listed in your `.amlignore` file.**
If an `.amlignore` doesn't already exist then one will be created on the first submission of a run to AML.

Once everything is configured, you can submit a run:

```@docs
Optics.Cloud.submit_run_to_AML
```

To retrieve outputs from your run simply write files to the `outputs/` directory and the files will automatically appear as part of the AML run.

### Examples

```julia
using Optics.Cloud

cache_run_config([subscription_id], [resource_group_name], [workspace_name], [compute_name], [path_to_config])

submit_run_to_AML("example-run", [path_to_script], ["--arg1", "1", "--arg2", "2"], nothing, [path_to_config])

submit_run_to_AML("example-hyperdrive-run", [path_to_script], ["--arg1", "1"], Dict("--arg2" => ["1", "2", "3"]), [path_to_config])
```

## Other Cloud Services

Currently no other services are supported, though it should be reasonably straightforward to add similar functionality to that for AML.
