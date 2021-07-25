"""Encapsulates the lens assembly and the display for one lenslet."""
struct LensletAssembly{T}
    lens::LensAssembly
    transform::Geometry.Transform
    display #composite emitter array maybe? See if emitter arrays are greedily or lazily constructed.
end

lens(a::LensLetAssembly) = a.lens
display(a::LensLetAssembly) = a.display
transform(a::LensLetAssembly) = a.transform
