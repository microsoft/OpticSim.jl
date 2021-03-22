struct AirType <: AbstractGlass end
Base.show(io::IO, ::AirType) = print(io, "Air")
glassid(::AirType) = GlassID(AIR, 0)
glassname(::AirType) = "GlassCat.Air"

isair(::AirType) = true
isair(::AbstractGlass) = false
isair(a::GlassID) = a.type === AIR

function info(io::IO, ::AirType)
    println(io, "GlassCat.Air")
    println(io, "Material representing air, RI is always 1.0 at system temperature and pressure, absorption is always 0.0.")
end
