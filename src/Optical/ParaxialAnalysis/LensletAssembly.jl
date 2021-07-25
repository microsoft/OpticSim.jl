"""Encapsulates the lens assembly and the display for one lenslet."""

struct Display{T}
    surface::Rectangle{T}
    #pixels::
end


struct LensletAssembly{T}
    lens::ParaxialLens{T}
    transform::Geometry.Transform{T}
    display #composite emitter array maybe? See if emitter arrays are greedily or lazily constructed.
end

lens(a::LensletAssembly) = a.lens
display(a::LensletAssembly) = a.display
transform(a::LensletAssembly) = a.transform
