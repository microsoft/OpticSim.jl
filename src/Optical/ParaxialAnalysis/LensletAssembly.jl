"""Encapsulates the lens assembly and the display for one lenslet."""
struct LensletAssembly{T}
    lens #need to figure out how to type this
    transform::Geometry.Transform
    display #composite emitter array maybe? See if emitter arrays are greedily or lazily constructed.
end

lens(a::LensletAssembly) = a.lens
display(a::LensletAssembly) = a.display
transform(a::LensletAssembly) = a.transform
