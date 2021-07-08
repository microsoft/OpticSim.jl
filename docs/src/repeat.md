# Repeating Structures

The Repeat module contains functions for creating regular repeated patterns. This could be pixels in a display grid, or mirrors in an active optics telescope.

Patterns are described by a set of lattice vectors. The ([`LatticeBasis`])(@ref sources)) type allows you to create bases in any dimension. The currently defined lattices are all 2D. There are functions for visualizing lattices defined in the visualization package.

```@example highlight
mdparse(@code_string OpticSim.Examples.drawhexregion()) # hide
```
```@example example
using OpticSim, OpticSim.Examples, OpticSim.Repeat; drawhexregion(); Vis.save("assets/repeat_example_hexregion.png") # hide
nothing #hide
```

![Hex regions visualization](assets/repeat_example_hexregion.png)