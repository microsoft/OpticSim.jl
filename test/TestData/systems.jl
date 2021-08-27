cooketriplet(::Type{T} = Float64) where {T<:Real} = AxisymmetricOpticalSystem{T}(
    DataFrame(SurfaceType = ["Object", "Standard", "Standard", "Standard", "Stop", "Standard", "Standard", "Image"],
              Radius = [Inf, 26.777, 66.604, -35.571, 35.571, 35.571, -26.777, Inf],
              Thickness = [Inf, 4.0, 2.0, 4.0, 2.0, 4.0, 44.748, missing],
              Material = [OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_SK16, OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_SF2, OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_SK16, OpticSim.GlassCat.Air, missing],
              SemiDiameter = [Inf, 8.580, 7.513, 7.054, 6.033, 7.003, 7.506, 15.0]))

cooketripletfirstelement(::Type{T} = Float64) where {T<:Real} = AxisymmetricOpticalSystem(
    DataFrame(SurfaceType = ["Object", "Standard", "Standard", "Image"],
              Radius = [Inf, -35.571, 35.571, Inf],
              Thickness = [Inf, 4.0, 44.748, missing],
              Material = [OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_SK16, OpticSim.GlassCat.Air, missing],
              SemiDiameter = [Inf, 7.054, 6.033, 15.0]))

convexplano(::Type{T} = Float64) where {T<:Real} = AxisymmetricOpticalSystem{T}(
    DataFrame(SurfaceType = ["Object", "Standard", "Standard", "Image"],
              Radius = [Inf, 60.0, Inf, Inf],
              Thickness = [Inf, 10.0, 57.8, missing],
              Material = [OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, missing],
              SemiDiameter = [Inf, 9.0, 9.0, 15.0]))

doubleconvex(frontradius::T,rearradius::T) where{T<:Real} =
AxisymmetricOpticalSystem{T}(
    DataFrame(SurfaceType = ["Object", "Standard", "Standard", "Image"],
              Radius = [T(Inf64), frontradius, rearradius, T(Inf64)],
              Thickness = [T(Inf64), T(10.0), T(57.8), missing],
              Material = [OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, missing],
              SemiDiameter = [T(Inf64), T(9.0), T(9.0), T(15.0)]),
              temperature = OpticSim.GlassCat.TEMP_REF_UNITFUL, pressure = T(OpticSim.GlassCat.PRESSURE_REF))

doubleconvex(::Type{T} = Float64; temperature::Unitful.Temperature = OpticSim.GlassCat.TEMP_REF_UNITFUL, pressure::Float64 = OpticSim.GlassCat.PRESSURE_REF) where {T<:Real} =
AxisymmetricOpticalSystem{T}(
    DataFrame(SurfaceType = ["Object", "Standard", "Standard", "Image"],
              Radius = [Inf64, 60, -60, Inf64],
              Thickness = [Inf64, 10.0, 57.8, missing],
              Material = [OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, missing],
              SemiDiameter = [Inf64, 9.0, 9.0, 15.0]),
              temperature = temperature, pressure = T(pressure))

doubleconcave(::Type{T} = Float64) where {T<:Real} = AxisymmetricOpticalSystem{T}(
    DataFrame(SurfaceType = ["Object", "Standard", "Standard", "Image"],
              Radius = [Inf64, -41.0, 41.0, Inf64],
              Thickness = [Inf64, 10.0, 57.8, missing],
              Material = [OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, missing],
              SemiDiameter = [Inf64, 9.0, 9.0, 15.0]))

planoconcaverefl(::Type{T} = Float64) where {T<:Real} = AxisymmetricOpticalSystem{T}(
    DataFrame(SurfaceType = ["Object", "Standard", "Standard", "Image"],
              Radius = [Inf64, Inf64, -41.0, Inf64],
              Thickness = [Inf64, 10.0, -57.8, missing],
              Material = [OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, missing],
              SemiDiameter = [Inf64, 9.0, 9.0, 25.0],
              Reflectance = [missing, missing, 1.0, missing]))

concaveplano(::Type{T} = Float64) where {T<:Real} = AxisymmetricOpticalSystem{T}(
    DataFrame(SurfaceType = ["Object", "Standard", "Standard", "Image"],
              Radius = [Inf64, -41.0, Inf64, Inf64],
              Thickness = [Inf64, 10.0, 57.8, missing],
              Material = [OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, missing],
              SemiDiameter = [Inf64, 9.0, 9.0, 15.0]))

planoplano(::Type{T} = Float64) where {T<:Real} = AxisymmetricOpticalSystem{T}(
    DataFrame(SurfaceType = ["Object", "Standard", "Standard", "Image"],
              Radius = [Inf64, Inf64, Inf64, Inf64],
              Thickness = [Inf64, 10.0, 57.8, missing],
              Material = [OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, missing],
              SemiDiameter = [Inf64, 9.0, 9.0, 15.0]))

zernikesystem() = CSGOpticalSystem(LensAssembly(zernikelens()), Rectangle(40.0, 40.0, [0.0, 0.0, 1.0], [0.0, 0.0, -67.8], interface = opaqueinterface()))
asphericsystem() = CSGOpticalSystem(LensAssembly(asphericlens()), Rectangle(25.0, 25.0, [0.0, 0.0, 1.0], [0.0, 0.0, -67.8], interface = opaqueinterface()))
conicsystemZ() = CSGOpticalSystem(LensAssembly(coniclensZ()), Rectangle(25.0, 25.0, [0.0, 0.0, 1.0], [0.0, 0.0, -67.8], interface = opaqueinterface()))
conicsystemQ() = CSGOpticalSystem(LensAssembly(coniclensQ()), Rectangle(25.0, 25.0, [0.0, 0.0, 1.0], [0.0, 0.0, -67.8], interface = opaqueinterface()))
qtypesystem() = CSGOpticalSystem(LensAssembly(qtypelens()), Rectangle(25.0, 25.0, [0.0, 0.0, 1.0], [0.0, 0.0, -67.8], interface = opaqueinterface()))
chebyshevsystem() = CSGOpticalSystem(LensAssembly(chebyshevlens()), Rectangle(40.0, 40.0, [0.0, 0.0, 1.0], [0.0, 0.0, -67.8], interface = opaqueinterface()))
gridsagsystem() = CSGOpticalSystem(LensAssembly(gridsaglens()), Rectangle(40.0, 40.0, [0.0, 0.0, 1.0], [0.0, 0.0, -67.8], interface = opaqueinterface()))
