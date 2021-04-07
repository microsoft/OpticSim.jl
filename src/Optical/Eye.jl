# MIT License

# Copyright (c) Microsoft Corporation.

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE

"""
    ModelEye(assembly::LensAssembly{T}, nsamples::Int = 17; pupil_radius::T = 3.0, detpixels::Int = 1000, transform::Transform{T} = identitytransform(T))

Geometrically accurate model of the human eye focussed at infinity with variable `pupil_radius`.
The eye is added to the provided `assembly` to create a [`CSGOpticalSystem`](@ref) with the retina of the eye as the detector.

The eye can be positioned in the scene using the `transform` argument and the resolution of the detector specified with `detpixels`.
By default the eye is directed along the positive z-axis with the vertex of the cornea at the origin.

`nsamples` determines the resolution at which accelerated surfaces within the eye are triangulated.
"""
function ModelEye(assembly::LensAssembly{T}; nsamples::Int = 17, pupil_radius::T = 3.0, detpixels::Int = 1000, transform::Transform{T} = identitytransform(T), detdatatype::Type{D} = Float32) where {T<:Real,D<:Real}
    cornea_front = AcceleratedParametricSurface(ZernikeSurface(6.9, radius = -7.8, conic = -0.5), nsamples, interface = FresnelInterface{T}(OpticSim.GlassCat.EYE.CORNEA, OpticSim.GlassCat.Air))
    cornea_rear = Plane(SVector{3,T}(0, 0, -1), SVector{3,T}(0, 0, -3.2))
    cornea_csg = csgintersection(cornea_front, cornea_rear)
    cornea = cornea_csg(transform)

    anterior_chamber_front = AcceleratedParametricSurface(ZernikeSurface(6.0, radius = -6.7, conic = -0.3), nsamples, interface = FresnelInterface{T}(OpticSim.GlassCat.EYE.AQUEOUS, OpticSim.GlassCat.EYE.CORNEA))
    anterior_chamber_rear = Plane(SVector{3,T}(0, 0, -1), SVector(0.0, 0.0, -3.2), vishalfsizeu = 6.0, vishalfsizev = 6.0, interface = FresnelInterface{T}(OpticSim.GlassCat.EYE.AQUEOUS, OpticSim.GlassCat.EYE.VITREOUS))
    anterior_chamber_csg = csgintersection(anterior_chamber_front, anterior_chamber_rear)
    anterior_chamber = leaf(anterior_chamber_csg, translation(0.0, 0.0, -0.52))(transform)

    pupil = Annulus(pupil_radius, 5.9, rotate(transform, SVector{3,T}(0, 0, 1)), transform * SVector{3,T}(0.0, 0.0, -3.62))

    lens_front = AcceleratedParametricSurface(ZernikeSurface(5.1, radius = -10.0), nsamples, interface = FresnelInterface{T}(OpticSim.GlassCat.EYE.LENS, OpticSim.GlassCat.EYE.VITREOUS))
    lens_rear = AcceleratedParametricSurface(ZernikeSurface(5.1, radius = 6.0, conic = -3.25), nsamples, interface = FresnelInterface{T}(OpticSim.GlassCat.EYE.LENS, OpticSim.GlassCat.EYE.VITREOUS))
    lens_barrel = Cylinder(5.0, 5.0, interface = FresnelInterface{T}(OpticSim.GlassCat.EYE.LENS, OpticSim.GlassCat.EYE.VITREOUS))
    lens_csg = csgintersection(csgdifference(leaf(lens_front), leaf(lens_rear, translation(0.0, 0.0, -3.7))), leaf(lens_barrel, translation(0.0, 0.0, -2.0)))
    lens = leaf(lens_csg, translation(0.0, 0.0, -3.72))(transform)

    vitreous_chamber_csg = csgintersection(csgintersection(leaf(Sphere(11.0, interface = FresnelInterface{T}(OpticSim.GlassCat.EYE.VITREOUS, OpticSim.GlassCat.Air, reflectance = 0.0, transmission = 0.0)), translation(0.0, 0.0, -13.138998863513297)), leaf(Plane(0.0, 0.0, 1.0, 0.0, 0.0, -3.72, interface = FresnelInterface{T}(OpticSim.GlassCat.EYE.VITREOUS, OpticSim.GlassCat.EYE.AQUEOUS)))), leaf(Plane(0.0, 0.0, -1.0, 0.0, 0.0, -24 + (11 - sqrt(11^2 - 10^2)), interface = NullInterface(T))))
    vitreous_chamber = vitreous_chamber_csg(transform)

    new_assembly = OpticSim.LensAssembly(assembly, cornea, anterior_chamber, pupil, vitreous_chamber, lens)
    retina = SphericalCap(11.0, asin(10 / 11), rotate(transform, SVector{3,T}(0, 0, 1)), transform * SVector{3,T}(0.0, 0.0, -24.0), interface = FresnelInterface{T}(OpticSim.GlassCat.EYE.VITREOUS, OpticSim.GlassCat.EYE.VITREOUS, reflectance = zero(T), transmission = zero(T)))

    return CSGOpticalSystem(new_assembly, retina, detpixels, detpixels, D)
end
export ModelEye

#! format: off
"""
    ArizonaEye(::Type{T} = Float64; accommodation::T = 0.0)

The popular Arizona eye model taken from [this definition](https://photonengr.com/wp-content/uploads/kbasefiles/ArizonaEyeModel.pdf).
The `accommodation` of the eye can be varied in this model.
Returns a `DataFrame` specifying the prescription of the eye model.
"""
function ArizonaEye(::Type{T} = Float64; accommodation::T = 0.0) where {T<:Real}
    # from https://photonengr.com/wp-content/uploads/kbasefiles/ArizonaEyeModel.pdf also in our documentation/papers directory
    return DataFrame(
        Name = ["Air", "Cornea", "Aqueous", "Lens", "Vitreous", "Retina"],
        Surface = [:Object, 1, 2, 3, 4, :Image],
        Radius = [Inf64, 7.8, 6.5, 12.0 - 0.4*accommodation, -5.22 + 0.2*accommodation, -13.4],
        Conic = [missing, -.25, -.25, -7.52 + 1.29*accommodation, -1.35 - 0.43*accommodation, missing],
        Material = [OpticSim.GlassCat.Air, OpticSim.GlassCat.modelglass(1.377, 57.1, 0.0), OpticSim.GlassCat.modelglass(1.337, 61.3, 0.0), OpticSim.GlassCat.modelglass(1.42 + 0.0026 * accommodation - 0.00022 * accommodation^2, 51.9, 0.0), OpticSim.GlassCat.modelglass(1.336, 61.1, 0.0), missing],
        SemiDiameter = [Inf64, 5.3, 5.3, 5.3, 5.3, 5.3],
        Thickness = [Inf64, 0.55, 2.97 - 0.04*accommodation, 3.77 + 0.04accommodation, 16.713, missing],
    )
end
export ArizonaEye
#! format: on
