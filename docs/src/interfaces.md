# Optical Interfaces

Every [`Surface`](@ref) must have an `OpticalInterface` associated with it to defined the behavior of any ray when it intersects that surface.

```@docs
OpticSim.OpticalInterface
OpticSim.NullInterface
FresnelInterface
ParaxialInterface
ThinGratingInterface
HologramInterface
MultiHologramInterface
```

The critical behavior of each interface is defined in the `processintersection` function:

```@docs
OpticSim.processintersection
```
