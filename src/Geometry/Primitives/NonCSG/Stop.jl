abstract type StopShape end
"""
    CircularStopShape <: StopShape
"""
abstract type CircularStopShape <: StopShape end
"""
    RectangularStopShape <: StopShape
"""
abstract type RectangularStopShape <: StopShape end

"""
    StopSurface{T} <: Surface{T}

Abstract type to encapsulate any surfaces acting as a stop.
"""
abstract type StopSurface{T} <: Surface{T} end

"""
    InfiniteStop{T,P<:StopShape} <: Surface{T}

Stop surface with infinite extent (outside of the aperture).
`P` refers to the shape of the aperture.
"""
struct InfiniteStop{T,P<:StopShape} <: StopSurface{T}
    plane::Plane{T,3}
    uvec::SVector{3,T}
    vvec::SVector{3,T}
    aphalfsizeu::T
    aphalfsizev::T
end
export InfiniteStop

function surfaceintersection(stop::InfiniteStop{T,CircularStopShape}, r::AbstractRay{T,3}) where {T<:Real}
    interval = surfaceintersection(stop.plane, r)
    if interval isa EmptyInterval{T}
        return EmptyInterval(T) # no ray plane intersection
    else
        intersect = closestintersection(interval)
        if intersect === nothing
            # inside plane but no hit
            return EmptyInterval(T)
        else
            d = point(intersect) - centroid(stop)
            if all(iszero.(d)) || norm(d) < stop.aphalfsizeu
                return EmptyInterval(T) # no ray plane intersection
            else
                return interval # forward interval from plane
            end
        end
    end
end

function surfaceintersection(stop::InfiniteStop{T,RectangularStopShape}, r::AbstractRay{T,3}) where {T<:Real}
    interval = surfaceintersection(stop.plane, r)
    if interval isa EmptyInterval{T}
        return EmptyInterval(T) # no ray plane intersection
    else
        intersect = closestintersection(interval)
        if intersect === nothing
            # inside plane but no hit
            return EmptyInterval(T)
        else
            dif = point(intersect) - centroid(stop)
            if abs(dot(dif, stop.uvec)) < stop.aphalfsizeu && abs(dot(dif, stop.vvec)) < stop.aphalfsizev
                return EmptyInterval(T) # no ray plane intersection
            else
                return interval # forward interval from plane
            end
        end
    end
end

interface(::InfiniteStop{T}) where {T<:Real} = opaqueinterface(T)
normal(r::InfiniteStop{T}) where {T<:Real} = normal(r.plane)
centroid(r::InfiniteStop{T}) where {T<:Real} = r.plane.pointonplane

"""
    FiniteStop{T,P<:StopShape,Q<:StopShape} <: Surface{T}

Stop surface with finite extent.
`P` refers to the shape of the aperture and `Q` represents the shape of the bounds of the stop surface.
"""
struct FiniteStop{T,P<:StopShape,Q<:StopShape} <: StopSurface{T}
    plane::Plane{T,3}
    uvec::SVector{3,T}
    vvec::SVector{3,T}
    innerhalfsizeu::T
    innerhalfsizev::T
    outerhalfsizeu::T
    outerhalfsizev::T
end
export FiniteStop

interface(::FiniteStop{T}) where {T<:Real} = opaqueinterface(T)
normal(r::FiniteStop{T}) where {T<:Real} = normal(r.plane)
normal(r::FiniteStop{T}, ::T, ::T) where {T<:Real} = normal(r.plane)
centroid(r::FiniteStop{T}) where {T<:Real} = r.plane.pointonplane

function surfaceintersection(stop::FiniteStop{T,CircularStopShape,CircularStopShape}, r::AbstractRay{T,3}) where {T<:Real}
    interval = surfaceintersection(stop.plane, r)
    if interval isa EmptyInterval{T}
        return EmptyInterval(T) # no ray plane intersection
    else
        intersect = closestintersection(interval)
        if intersect === nothing
            # inside plane but no hit
            return EmptyInterval(T)
        else
            d = point(intersect) - centroid(stop)
            if !all(iszero.(d)) && stop.innerhalfsizeu <= norm(d) <= stop.outerhalfsizeu
                return interval # forward interval from plane
            else
                return EmptyInterval(T) # no ray plane intersection
            end
        end
    end
end

function surfaceintersection(stop::FiniteStop{T,CircularStopShape,RectangularStopShape}, r::AbstractRay{T,3}) where {T<:Real}
    interval = surfaceintersection(stop.plane, r)
    if interval isa EmptyInterval{T}
        return EmptyInterval(T) # no ray plane intersection
    else
        intersect = closestintersection(interval)
        if intersect === nothing
            # inside plane but no hit
            return EmptyInterval(T)
        else
            p = point(intersect)
            d = p - centroid(stop)
            if all(iszero.(d)) || norm(d) < stop.innerhalfsizeu
                return EmptyInterval(T)
            else
                dif = p - centroid(stop)
                du = abs(dot(dif, stop.uvec))
                dv = abs(dot(dif, stop.vvec))
                if du > stop.outerhalfsizeu || dv > stop.outerhalfsizev
                    return EmptyInterval(T)
                else
                    return interval
                end
            end
        end
    end
end

function surfaceintersection(stop::FiniteStop{T,RectangularStopShape,RectangularStopShape}, r::AbstractRay{T,3}) where {T<:Real}
    interval = surfaceintersection(stop.plane, r)
    if interval isa EmptyInterval{T}
        return EmptyInterval(T) # no ray plane intersection
    else
        intersect = closestintersection(interval)
        if intersect === nothing
            # inside plane but no hit
            return EmptyInterval(T)
        else
            p = point(intersect)
            du = abs(dot(p - centroid(stop), stop.uvec))
            dv = abs(dot(p - centroid(stop), stop.vvec))
            if (du < stop.innerhalfsizeu && dv < stop.innerhalfsizev) || du > stop.outerhalfsizeu || dv > stop.outerhalfsizev
                return EmptyInterval(T)
            else
                return interval
            end
        end
    end
end

"""
    RectangularAperture(aphalfsizeu::T, aphalfsizev::T, surfacenormal::SVector{3,T}, centrepoint::SVector{3,T}; rotationvec::SVector{3,T} = [0.0, 1.0, 0.0])

Creates a rectangular aperture in a plane i.e. `InfiniteStop{T,RectangularStopShape}`.
The rotation of the rectangle around its normal is defined by `rotationvec`.
`rotationvec×surfacenormal` is taken as the vector along the u axis.
"""
function RectangularAperture(aphalfsizeu::T, aphalfsizev::T, surfacenormal::SVector{3,T}, centrepoint::SVector{3,T}; rotationvec::SVector{3,T} = SVector{3,T}(0.0, 1.0, 0.0)) where {T<:Real}
    @assert aphalfsizeu > 0 && aphalfsizev > 0
    p = Plane(surfacenormal, centrepoint, interface = opaqueinterface(T))
    n̂ = normal(p)
    if abs(dot(rotationvec, n̂)) == one(T)
        rotationvec = SVector{3,T}(1.0, 0.0, 0.0)
    end
    uvec = normalize(cross(normalize(rotationvec), n̂))
    vvec = normalize(cross(n̂, uvec))
    return InfiniteStop{T,RectangularStopShape}(p, uvec, vvec, aphalfsizeu, aphalfsizev)
end

"""
    RectangularAperture(innerhalfsizeu::T, innerhalfsizev::T, outerhalfsizeu::T, outerhalfsizev::T, surfacenormal::SVector{3,T}, centrepoint::SVector{3,T}; rotationvec::SVector{3,T} = [0.0, 1.0, 0.0])

Creates a rectangular aperture in a rectangle i.e. `FiniteStop{T,RectangularStopShape,RectangularStopShape}`.
The rotation of the rectangle around its normal is defined by `rotationvec`.
`rotationvec×surfacenormal` is taken as the vector along the u axis.
"""
function RectangularAperture(innerhalfsizeu::T, innerhalfsizev::T, outerhalfsizeu::T, outerhalfsizev::T, surfacenormal::SVector{3,T}, centrepoint::SVector{3,T}; rotationvec::SVector{3,T} = SVector{3,T}(0.0, 1.0, 0.0)) where {T<:Real}
    @assert innerhalfsizeu < outerhalfsizeu && innerhalfsizev < outerhalfsizev && innerhalfsizeu > 0 && innerhalfsizev > 0
    p = Plane(surfacenormal, centrepoint, interface = opaqueinterface(T))
    n̂ = normal(p)
    if abs(dot(rotationvec, n̂)) == one(T)
        rotationvec = SVector{3,T}(1.0, 0.0, 0.0)
    end
    uvec = normalize(cross(normalize(rotationvec), n̂))
    vvec = normalize(cross(n̂, uvec))
    return FiniteStop{T,RectangularStopShape,RectangularStopShape}(p, uvec, vvec, innerhalfsizeu, innerhalfsizev, outerhalfsizeu, outerhalfsizev)
end

export RectangularAperture

"""
    CircularAperture(radius::T, surfacenormal::SVector{3,T}, centrepoint::SVector{3,T})

Creates a circular aperture in a plane i.e. `InfiniteStop{T,CircularStopShape}`.
"""
function CircularAperture(radius::T, surfacenormal::SVector{3,T}, centrepoint::SVector{3,T}) where {T<:Real}
    @assert radius > 0
    p = Plane(surfacenormal, centrepoint, interface = opaqueinterface(T))
    return InfiniteStop{T,CircularStopShape}(p, SVector{3,T}(0.0, 0.0, 0.0), SVector{3,T}(0.0, 0.0, 0.0), radius, radius)
end

"""
    CircularAperture(radius::T, outerhalfsizeu::T, outerhalfsizev::T, surfacenormal::SVector{3,T}, centrepoint::SVector{3,T}; rotationvec::SVector{3,T} = [0.0, 1.0, 0.0])

Creates a circular aperture in a rectangle i.e. `FiniteStop{T,CircularStopShape,RectangularStopShape}`.
The rotation of the rectangle around its normal is defined by `rotationvec`.
`rotationvec×surfacenormal` is taken as the vector along the u axis.
"""
function CircularAperture(radius::T, outerhalfsizeu::T, outerhalfsizev::T, surfacenormal::SVector{3,T}, centrepoint::SVector{3,T}; rotationvec::SVector{3,T} = SVector{3,T}(0.0, 1.0, 0.0)) where {T<:Real}
    @assert radius < outerhalfsizeu && radius < outerhalfsizev && radius > 0
    p = Plane(surfacenormal, centrepoint, interface = opaqueinterface(T))
    n̂ = normal(p)
    if abs(dot(rotationvec, n̂)) == one(T)
        rotationvec = SVector{3,T}(1.0, 0.0, 0.0)
    end
    uvec = normalize(cross(normalize(rotationvec), n̂))
    vvec = normalize(cross(n̂, uvec))
    return FiniteStop{T,CircularStopShape,RectangularStopShape}(p, uvec, vvec, radius, radius, outerhalfsizeu, outerhalfsizev)
end
export CircularAperture

"""
    Annulus(innerradius::T, outerradius::T, surfacenormal::SVector{3,T}, centrepoint::SVector{3,T})

Creates a circular aperture in a circle i.e. `FiniteStop{T,CircularStopShape,CircularStopShape}`.
"""
function Annulus(innerradius::T, outerradius::T, surfacenormal::SVector{3,T}, centrepoint::SVector{3,T}) where {T<:Real}
    @assert outerradius > innerradius && innerradius > 0
    p = Plane(surfacenormal, centrepoint, interface = opaqueinterface(T))
    n̂ = normal(p)
    rotationvec = SVector{3,T}(0.0, 1.0, 0.0)
    if abs(dot(rotationvec, n̂)) == one(T)
        rotationvec = SVector{3,T}(1.0, 0.0, 0.0)
    end
    uvec = normalize(cross(normalize(rotationvec), n̂))
    vvec = normalize(cross(n̂, uvec))
    return FiniteStop{T,CircularStopShape,CircularStopShape}(p, uvec, vvec, innerradius, innerradius, outerradius, outerradius)
end
export Annulus

# override makemesh so that we use as few triangles as possible - these can't be CSG objects so we don't need small triangles

function makemesh(s::FiniteStop{T,CircularStopShape,CircularStopShape}, subdivisions::Int = 30) where {T<:Real}
    point(u::T, v::T) = centroid(s) + ((s.innerhalfsizeu + v * (s.outerhalfsizeu - s.innerhalfsizeu)) * (sin(u) * s.uvec + cos(u) * s.vvec))
    dθ = T(2π) / subdivisions
    tris = Vector{Triangle{T}}(undef, subdivisions * 2)
    @inbounds for i in 0:(subdivisions - 1)
        θ1 = i * dθ - π
        θ2 = (i + 1) * dθ - π
        ip1 = point(θ1, zero(T))
        ip2 = point(θ2, zero(T))
        op1 = point(θ1, one(T))
        op2 = point(θ2, one(T))
        tris[i * 2 + 1] = Triangle(ip1, op2, op1)
        tris[i * 2 + 2] = Triangle(ip1, ip2, op2)
    end
    return TriangleMesh(tris)
end

function makemesh(s::FiniteStop{T,CircularStopShape,RectangularStopShape}, subdivisions::Int = 30) where {T<:Real}
    point(u::T) = centroid(s) + (s.innerhalfsizeu * (sin(u) * s.uvec + cos(u) * s.vvec))
    outeru = s.outerhalfsizeu * s.uvec
    outerv = s.outerhalfsizev * s.vvec
    o00 = centroid(s) - outeru - outerv
    o01 = centroid(s) - outeru + outerv
    o11 = centroid(s) + outeru + outerv
    o10 = centroid(s) + outeru - outerv
    dθ = T(2π) / subdivisions
    tris = Vector{Triangle{T}}(undef, subdivisions + 4)
    @inbounds for i in 0:(subdivisions - 1)
        θ1 = i * dθ - π
        θ2 = (i + 1) * dθ - π
        ip1 = point(θ1)
        ip2 = point(θ2)
        closestcorner = θ1 < -π / 2 ? o00 : θ1 < 0.0 ? o01 : θ1 < π / 2 ? o11 : o10
        tris[i + 1] = Triangle(ip1, ip2, closestcorner)
    end
    tris[end - 3] = Triangle(o00, point(ceil(subdivisions / 4) * dθ - π), o01)
    tris[end - 2] = Triangle(o01, point(ceil(subdivisions / 2) * dθ - π), o11)
    tris[end - 1] = Triangle(o11, point(ceil(3 * subdivisions / 4) * dθ - π), o10)
    tris[end] = Triangle(o10, point(-π), o00)
    return TriangleMesh(tris)
end

function makemesh(s::FiniteStop{T,RectangularStopShape,RectangularStopShape}, ::Int = 0) where {T<:Real}
    outeru = s.outerhalfsizeu * s.uvec
    outerv = s.outerhalfsizev * s.vvec
    inneru = s.innerhalfsizeu * s.uvec
    innerv = s.innerhalfsizev * s.vvec
    o00 = centroid(s) - outeru - outerv
    o01 = centroid(s) - outeru + outerv
    o11 = centroid(s) + outeru + outerv
    o10 = centroid(s) + outeru - outerv
    i00 = centroid(s) - inneru - innerv
    i01 = centroid(s) - inneru + innerv
    i11 = centroid(s) + inneru + innerv
    i10 = centroid(s) + inneru - innerv
    u00 = centroid(s) - outeru - innerv
    u01 = centroid(s) - outeru + innerv
    u11 = centroid(s) + outeru + innerv
    u10 = centroid(s) + outeru - innerv
    t1 = Triangle(o00, o10, u10)
    t2 = Triangle(o00, u10, u00)
    t3 = Triangle(u00, i01, u01)
    t4 = Triangle(u00, i00, i01)
    t5 = Triangle(u01, o11, o01)
    t6 = Triangle(u01, u11, o11)
    t7 = Triangle(i11, i10, u11)
    t8 = Triangle(u11, i10, u10)
    return TriangleMesh([t1, t2, t3, t4, t5, t6, t7, t8])
end

makemesh(::InfiniteStop, ::Int = 0) = nothing
