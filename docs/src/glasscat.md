# GlassCat

This submodule is used to download, parse, install and manage AGF glass specifications for use in OpticSim.

The central configuration file for GlassCat is located at `src/GlassCat/data/sources.txt`, which ships with the
following default entries.

```
NIKON a49714470fa875ad4fd8d11cbc0edbf1adfe877f42926d7612c1649dd9441e75 https://www.nikon.com/products/components/assets/pdf/nikon_zemax_data.zip
OHARA 0c9021bf11b8d4e660012191818685ad3110d4f9490699cabdc89aae1fd26d2e https://www.oharacorp.com/xls/OHARA_201130_CATALOG.zip
HOYA b02c203e5a5b7a8918cc786badf0a9db1fe2572372c1c163dc306b0a8a908854 http://www.hoya-opticalworld.com/common/agf/HOYA20210105.agf
SCHOTT e9aabbb8ebff116ba0c106a71afd86e72f2a397ac9bc447469129e325e795f4e https://www.schott.com/d/advanced_optics/6959f9a4-0e4f-4ef2-a302-2347468a82f5/1.31/schott-optical-glass-overview-zemax-format.zip
SUMITA c1093e42a1d08acbe30698aba730161e3b43f8f0d50533f65de8b6b11100fdc8 https://www.sumita-opt.co.jp/en/wp/wp-content/themes/sumita-opt-en/download-files.php files%5B%5D=new_sumita.agf&category=data
```

Each line corresponds to one AGF source, which is described by 2 to 4 space-delimited columns. The first column provides
the installed module name for the catalog, e.g. `GlassCat.NIKON`. The second column is the expected SHA256 checksum for
the AGF file.

The final two columns are optional, specifying download instructions for acquiring the zipped AGF files
automatically from the web. The fourth column allows us to use POST requests to acquire files from interactive sites.

When `] build OpticSim` is run, the sources are verified and parsed into corresponding Julia files. These are then
included in OpticSim via `AGFGlassCat.jl`. These steps are run automatically when the package is first installed using
`] add OpticSim`, creating a sufficient working environment for our examples and tests.

## Adding glass catalogs
`sources.txt` can be edited freely to add more glass catalogs. However, this is a somewhat tedious process, so we have a
convenience function for adding a locally downloaded AGF file to the source list.

```@docs
OpticSim.GlassCat.add_agf
```

## Using installed glasses

Glass types are accessed like so: `OpticSim.GlassCat.CATALOG_NAME.GLASS_NAME`, e.g.

```julia
OpticSim.GlassCat.SUMITA.LAK7
OpticSim.GlassCat.SCHOTT.PK3
```

All glasses and catalogs are exported in their respective modules, so it is possible to invoke `using` calls for convenience, e.g.

```julia
using OpticSim
GlassCat.SUMITA.LAK7
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
OpticSim.GlassCat.glassnames
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
OpticSim.GlassCat.drawglassmap
```

---

```@docs
OpticSim.GlassCat.verify_sources!
OpticSim.GlassCat.verify_source
OpticSim.GlassCat.download_source
```
