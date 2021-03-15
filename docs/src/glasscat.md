# GlassCat

Julia module to import and use AGF glass specifications.

The AGF glass catalog is stored in a Julia source file `AGFOptics.GlassCat.jl`. If the source file isn't already present then it will be generated automatically when `using Optics` is called.
Any AGF file in the directory specified by environment variable `GLASS_CAT_DIR` will be parsed, or if this isn't specified the directory `~/Dev/AGFGlassCat` is used.
If a new AGF file is added to this directory, manually delete `AGFOptics.GlassCat.jl` and execute `using Optics` in the Julia REPL to trigger the update of `AGFOptics.GlassCat.jl`.

Optical systems in the examples file expect glass catalogs to have specific names. For example, if you download the Schott glass catalog from the Schott website and unzip the file it will have a name like this: schottzemax-20190109.agf. 

By default the AGF parsing code gives the glass catalog in `AGFOptics.GlassCat.jl` the name associated with the file. But the examples expect the Schott glass catalog to be named SCHOTT. If you want these examples to run correctly you must rename this file to SCHOTT.agf Similarly these glass catalogs should be renamed: [NIKON](https://www.nikon.com/products/optical-glass/assets/pdf/nikon_zemax_data.zip),[OHARA](https://www.oharacorp.com/xls/OHARA_201130_CATALOG.zip),[HOYA](https://hoyaoptics.com/wp-content/uploads/2019/10/HOYA20170401.zip)

Glass types are accessed like so: `Optics.GlassCat.CATALOG_NAME.GLASS_NAME`, e.g.

```julia
Optics.GlassCat.Sumita.LAK7
Optics.GlassCat.SCHOTT.PK3
```

Autocompletion can be used to see available catalogs and glasses. All catalog glasses are of type [`Optics.GlassCat.Glass`](@ref).
Note that special characters in glass/catalog names are replaced with `_`.
There is a special type and constant value for air: [`Optics.GlassCat.Air`](@ref).

[Unitful.jl](https://github.com/PainterQubits/Unitful.jl) is used to manage units, meaning any valid unit can be used for all arguments, e.g., wavelength can be passed in as μm or nm (or cm, mm, m, etc.).
Non-unitful options are also available, in which case units are assumed to be μm, °C and Atm for length, temperature and pressure respectively.

`TEMP_REF` and `PRESSURE_REF` are constants:

```julia
const TEMP_REF = 20.0 # °C
const PRESSURE_REF = 1.0 # Atm
```

## Types

```@docs
Optics.GlassCat.AbstractGlass
Optics.GlassCat.Glass
Optics.GlassCat.Air
Optics.GlassCat.GlassID
```

## Functions

```@docs
Optics.GlassCat.index
Optics.GlassCat.absairindex
Optics.GlassCat.absorption
```

---

```@docs
Optics.GlassCat.glassfromMIL
Optics.GlassCat.modelglass
```

---

```@docs
Optics.GlassCat.glasscatalogs
Optics.GlassCat.glasses
Optics.GlassCat.info
Optics.GlassCat.findglass
Optics.GlassCat.isair
```

---

```@docs
Optics.GlassCat.glassname
Optics.GlassCat.glassid
Optics.GlassCat.glassforid
```

---

```@docs
Optics.GlassCat.polyfit_indices
Optics.GlassCat.plot_indices
```
