# GlassCat

Julia module for importing and using AGF glass specifications.

The entire AGF glass catalog is specified in `AGFGlassCat.jl`. This Julia source file is generated automatically when `] build OpticSim` is called. The build script downloads AGF files to `deps/downloads/glasscat/` and then uses these to generate corresponding Julia source files at `src/GlassCat/data/`. These steps are run automatically on setup when the package is first installed using `] add OpticSim`, creating a sufficient working environment for our example/test code.

Adding new AGF sources is done by editing `deps/sources.txt`. Minimally, you must provide a name (e.g. SCHOTT) and sha256sum for the AGF file, which can then be placed manually into `deps/downloads/glasscat/[NAME].agf`. Instead of manually sourcing the AGF file, you can also provide a download link for the build script. `deps/sources.txt` already contains examples of all possible use cases. After updating the file, execute `] build OpticSim` to rebuild `AGFGlassCat.jl`.

Source names will be used as module names, so follow the standard convections: alphanumeric, no leading numbers, begin with an uppercase letter. Furthermore, optical systems in the examples file expect glass catalogs to have specific names. The default setup includes HOYA, NIKON, OHARA, SCHOTT and Sumita; changing these names could break some examples.

Glass types are accessed like so: `OpticSim.GlassCat.CATALOG_NAME.GLASS_NAME`, e.g.

```julia
OpticSim.GlassCat.Sumita.LAK7
OpticSim.GlassCat.SCHOTT.PK3
```

All glasses and catalogs are exported in their respective modules, so it is possible to invoke `using` calls for convenience, e.g.

```julia
using OpticSim
GlassCat.Sumita.LAK7
using OpticSim.GlassCat
SCHOTT.PK3
using OpticsSim.GlassCat.SCHOTT
N_BK7
```

Autocompletion can be used to see available catalogs and glasses. All catalog glasses are of type [`OpticSim.GlassCat.Glass`](@ref).
Note that special characters in glass/catalog names are replaced with `_`.
There is a special type and constant value for air: [`OpticSim.GlassCat.Air`](@ref).

[Unitful.jl](https://github.com/PainterQubits/Unitful.jl) is used to manage units, meaning any valid unit can be used for all arguments, e.g., wavelength can be passed in as μm or nm (or cm, mm, m, etc.).
Non-unitful options are also available, in which case units are assumed to be μm, °C and Atm for length, temperature and pressure respectively.

`TEMP_REF` and `PRESSURE_REF` are constants:

```julia
const TEMP_REF = 20.0 # °C
const PRESSURE_REF = 1.0 # Atm
```

## Types

```@docs
OpticSim.GlassCat.AbstractGlass
OpticSim.GlassCat.Glass
OpticSim.GlassCat.Air
OpticSim.GlassCat.GlassID
```

## Functions

```@docs
OpticSim.GlassCat.index
OpticSim.GlassCat.absairindex
OpticSim.GlassCat.absorption
```

---

```@docs
OpticSim.GlassCat.glassfromMIL
OpticSim.GlassCat.modelglass
```

---

```@docs
OpticSim.GlassCat.glasscatalogs
OpticSim.GlassCat.glasses
OpticSim.GlassCat.info
OpticSim.GlassCat.findglass
OpticSim.GlassCat.isair
```

---

```@docs
OpticSim.GlassCat.glassname
OpticSim.GlassCat.glassid
OpticSim.GlassCat.glassforid
```

---

```@docs
OpticSim.GlassCat.polyfit_indices
OpticSim.GlassCat.plot_indices
```
