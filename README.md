<p align="center">
  <a href="https://microsoft.github.io/OpticSim.jl/dev/">
    <img src=docs/src/assets/logo.svg height=128px style="text-align:center">
  </a>
</p>

# OpticSim.jl

<table>
<thead>
  <tr>
    <th>Documentation</th>
    <th>Build Status</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td>
      <a href="https://microsoft.github.io/OpticSim.jl/stable/">
        <img src="https://img.shields.io/badge/docs-stable-blue.svg" alt="docs stable">
      </a>
      <a href="https://microsoft.github.io/OpticSim.jl/dev/">
        <img src="https://img.shields.io/badge/docs-dev-blue.svg" alt="docs dev">
      </a>
    </td>
    <td>
      <a href="https://github.com/microsoft/OpticSim.jl/actions/workflows/CI.yml">
        <img src="https://github.com/microsoft/OpticSim.jl/workflows/CI/badge.svg" alt="CI action">
      </a>
      <a href="https://codecov.io/gh/microsoft/OpticSim.jl">
        <img src="https://codecov.io/gh/microsoft/OpticSim.jl/branch/main/graph/badge.svg?token=9QxvIHt5F5" alt="codecov">
      </a>
    </td>
  </tr>
</tbody>
</table>

OpticSim.jl is a [Julia](https://julialang.org/) package for geometric optics (ray tracing) simulation and optimization of complex optical systems developed by the Microsoft Research Interactive Media Group and the Microsoft Hardware Architecture Incubation Team (HART).

It is designed to allow optical engineers to create optical systems procedurally and then to simulate and optimize them. Unlike Zemax, Code V, or other interactive optical design systems OpticSim.jl has limited support for interactivity, primarily in the tools for visualizing optical systems.

A large variety of surface types are supported, and these can be composed into complex 3D objects through the use of constructive solid geometry (CSG). A substantial catalog of optical materials is provided through the GlassCat submodule.

This software provides extensive control over the modelling, simulation, visualization and optimization of optical systems. It is especially suited for designs that have a procedural architecture.

# Installation

Before you can use the software you will need to download glass files. See the documentation for detailed information about how to do this.

*Warning*: During installation OpticSim automatically downloads glass catalogs from a variety of public sources. The Schott website keeps moving their catalog on their website so our software can't find it to download. This caused all our examples to fail because they use Schott glasses. We have replaced the glasses in the examples with hard coded glass files so all examples now work (if they don't please file an issue). 

**If you want to use Schott glasses you will have to manually download and install the Schott catalog.**

## Contributing

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)

This project welcomes contributions and suggestions.  Most contributions require you to agree to a
Contributor License Agreement (CLA) declaring that you have the right to, and actually do, grant us
the rights to use your contribution. For details, visit <https://cla.opensource.microsoft.com>.

When you submit a pull request, a CLA bot will automatically determine whether you need to provide
a CLA and decorate the PR appropriately (e.g., status check, comment). Simply follow the instructions
provided by the bot. You will only need to do this once across all repositories using our CLA.

This project has adopted the [Microsoft Open Source Code of Conduct](https://opensource.microsoft.com/codeofconduct/).
For more information see the [Code of Conduct FAQ](https://opensource.microsoft.com/codeofconduct/faq/) or
contact [opencode@microsoft.com](mailto:opencode@microsoft.com) with any additional questions or comments.

## Trademarks

This project may contain trademarks or logos for projects, products, or services. Authorized use of Microsoft 
trademarks or logos is subject to and must follow 
[Microsoft's Trademark & Brand Guidelines](https://www.microsoft.com/en-us/legal/intellectualproperty/trademarks/usage/general).
Use of Microsoft trademarks or logos in modified versions of this project must not cause confusion or imply Microsoft sponsorship.
Any use of third-party trademarks or logos are subject to those third-party's policies.
