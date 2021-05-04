# Optical Systems

## Assemblies

All systems are made up of a `LensAssembly` which contains all the optical components in the system, excluding any sources (see [Emitters](@ref)) and the detector (see below).

```@docs
LensAssembly
```

## Images

The detector image is stored within the system as a `HierarchicalImage` for memory efficiency.

```@docs
HierarchicalImage
OpticSim.reset!
OpticSim.sum!
```

## Systems

There are two types of `AbstractOpticalSystem` which can be used depending on the requirements.

```@docs
AbstractOpticalSystem
CSGOpticalSystem
AxisymmetricOpticalSystem
temperature
pressure
detectorimage
resetdetector!
assembly
semidiameter
```

## Tracing

We can trace an individual [`OpticalRay`](@ref) through the system (or directly through a [`LensAssembly`](@ref)), or we can trace using an [`OpticalRayGenerator`](@ref) to create a large number of rays.

```@docs
trace
traceMT
tracehits
tracehitsMT
OpticSim.LensTrace
```
