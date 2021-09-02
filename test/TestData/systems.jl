function cooketriplet(::Type{T} = Float64) where {T<:Real}
    return AxisymmetricOpticalSystem{T}(
        DataFrame(
            SurfaceType = ["Object", "Standard", "Standard", "Standard", "Stop", "Standard", "Standard", "Image"],
            Radius = [Inf, 26.777, 66.604, -35.571, 35.571, 35.571, -26.777, Inf],
            Thickness = [Inf, 4.0, 2.0, 4.0, 2.0, 4.0, 44.748, missing],
            Material = [Air, SCHOTT.N_SK16, Air, SCHOTT.N_SF2, Air, SCHOTT.N_SK16, Air, missing],
            SemiDiameter = [Inf, 8.580, 7.513, 7.054, 6.033, 7.003, 7.506, 15.0]
        )
    )
end

function cooketripletfirstelement(::Type{T} = Float64) where {T<:Real}
    return AxisymmetricOpticalSystem{T}(
        DataFrame(
            SurfaceType = ["Object", "Standard", "Standard", "Image"],
            Radius = [Inf, -35.571, 35.571, Inf],
            Thickness = [Inf, 4.0, 44.748, missing],
            Material = [Air, SCHOTT.N_SK16, Air, missing],
            SemiDiameter = [Inf, 7.054, 6.033, 15.0]
        )
    )
end

function convexplano(::Type{T} = Float64) where {T<:Real}
    return AxisymmetricOpticalSystem{T}(
        DataFrame(
            SurfaceType = ["Object", "Standard", "Standard", "Image"],
            Radius = [Inf, 60.0, Inf, Inf],
            Thickness = [Inf, 10.0, 57.8, missing],
            Material = [Air, SCHOTT.N_BK7, Air, missing],
            SemiDiameter = [Inf, 9.0, 9.0, 15.0]
        )
    )
end

function doubleconvex(
    frontradius::T, rearradius::T;
    temperature::Unitful.Temperature = GlassCat.TEMP_REF_UNITFUL, pressure::T = T(PRESSURE_REF)
) where {T<:Real}
    return AxisymmetricOpticalSystem{T}(
        DataFrame(
            SurfaceType = ["Object", "Standard", "Standard", "Image"],
            Radius = [T(Inf64), frontradius, rearradius, T(Inf64)],
            Thickness = [T(Inf64), T(10.0), T(57.8), missing],
            Material = [Air, SCHOTT.N_BK7, Air, missing],
            SemiDiameter = [T(Inf64), T(9.0), T(9.0), T(15.0)]
        );
        temperature,
        pressure
    )
end

function doubleconvex(
    ::Type{T} = Float64; temperature::Unitful.Temperature = GlassCat.TEMP_REF_UNITFUL, pressure::T = T(PRESSURE_REF)
) where {T<:Real}
    return AxisymmetricOpticalSystem{T}(
        DataFrame(
            SurfaceType = ["Object", "Standard", "Standard", "Image"],
            Radius = [Inf64, 60, -60, Inf64],
            Thickness = [Inf64, 10.0, 57.8, missing],
            Material = [Air, SCHOTT.N_BK7, Air, missing],
            SemiDiameter = [Inf64, 9.0, 9.0, 15.0]
        );
        temperature,
        pressure
    )
end

function doubleconcave(::Type{T} = Float64) where {T<:Real}
    return AxisymmetricOpticalSystem{T}(
        DataFrame(SurfaceType = ["Object", "Standard", "Standard", "Image"],
            Radius = [Inf64, -41.0, 41.0, Inf64],
            Thickness = [Inf64, 10.0, 57.8, missing],
            Material = [Air, SCHOTT.N_BK7, Air, missing],
            SemiDiameter = [Inf64, 9.0, 9.0, 15.0]
        )
    )
end

function planoconcaverefl(::Type{T} = Float64) where {T<:Real}
    return AxisymmetricOpticalSystem{T}(
        DataFrame(SurfaceType = ["Object", "Standard", "Standard", "Image"],
            Radius = [Inf64, Inf64, -41.0, Inf64],
            Thickness = [Inf64, 10.0, -57.8, missing],
            Material = [Air, SCHOTT.N_BK7, Air, missing],
            SemiDiameter = [Inf64, 9.0, 9.0, 25.0],
            Reflectance = [missing, missing, 1.0, missing]
        )
    )
end

function concaveplano(::Type{T} = Float64) where {T<:Real}
    return AxisymmetricOpticalSystem{T}(
        DataFrame(SurfaceType = ["Object", "Standard", "Standard", "Image"],
            Radius = [Inf64, -41.0, Inf64, Inf64],
            Thickness = [Inf64, 10.0, 57.8, missing],
            Material = [Air, SCHOTT.N_BK7, Air, missing],
            SemiDiameter = [Inf64, 9.0, 9.0, 15.0]
        )
    )
end

function planoplano(::Type{T} = Float64) where {T<:Real}
    return AxisymmetricOpticalSystem{T}(
        DataFrame(SurfaceType = ["Object", "Standard", "Standard", "Image"],
            Radius = [Inf64, Inf64, Inf64, Inf64],
            Thickness = [Inf64, 10.0, 57.8, missing],
            Material = [Air, SCHOTT.N_BK7, Air, missing],
            SemiDiameter = [Inf64, 9.0, 9.0, 15.0]
        )
    )
end

"""
    planarshapes()

Returns a CSGOpticalSystem containing a vertical stack (along the z-axis) of alternating planar shapes and a detector
placed at the end of the stack. This is for the purpose of benchmarking allocations when planar shapes are stored in a
shared Vector within LensAssembly.
"""
function planarshapes()
    interface = FresnelInterface{Float64}(SCHOTT.N_BK7, Air; interfacemode=Transmit)
    s = 10.0
    surfacenormal = SVector(0., 0., 1.)

    elements = []
    for i = 1:9
        centrepoint = SVector(0., 0., -i)
        if i % 3 == 1
            push!(elements, Ellipse(s, s, surfacenormal, centrepoint; interface))
        elseif i % 3 == 2
            push!(elements, Hexagon(s, surfacenormal, centrepoint; interface))
        elseif i % 3 == 0
            push!(elements, Rectangle(s, s, surfacenormal, centrepoint; interface))
        end
    end

    lensassembly = LensAssembly(elements...)
    detector = Rectangle(s, s, surfacenormal, SVector(0., 0., -10.); interface = opaqueinterface())

    return CSGOpticalSystem(lensassembly, detector)
end

zernikesystem() = CSGOpticalSystem(LensAssembly(zernikelens()), Rectangle(40.0, 40.0, [0.0, 0.0, 1.0], [0.0, 0.0, -67.8], interface = opaqueinterface()))
asphericsystem() = CSGOpticalSystem(LensAssembly(asphericlens()), Rectangle(25.0, 25.0, [0.0, 0.0, 1.0], [0.0, 0.0, -67.8], interface = opaqueinterface()))
conicsystemZ() = CSGOpticalSystem(LensAssembly(coniclensZ()), Rectangle(25.0, 25.0, [0.0, 0.0, 1.0], [0.0, 0.0, -67.8], interface = opaqueinterface()))
conicsystemQ() = CSGOpticalSystem(LensAssembly(coniclensQ()), Rectangle(25.0, 25.0, [0.0, 0.0, 1.0], [0.0, 0.0, -67.8], interface = opaqueinterface()))
qtypesystem() = CSGOpticalSystem(LensAssembly(qtypelens()), Rectangle(25.0, 25.0, [0.0, 0.0, 1.0], [0.0, 0.0, -67.8], interface = opaqueinterface()))
chebyshevsystem() = CSGOpticalSystem(LensAssembly(chebyshevlens()), Rectangle(40.0, 40.0, [0.0, 0.0, 1.0], [0.0, 0.0, -67.8], interface = opaqueinterface()))
gridsagsystem() = CSGOpticalSystem(LensAssembly(gridsaglens()), Rectangle(40.0, 40.0, [0.0, 0.0, 1.0], [0.0, 0.0, -67.8], interface = opaqueinterface()))
