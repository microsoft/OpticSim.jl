module Geometry     

import ...OpticSim

using LinearAlgebra
using StaticArrays

struct Vec3D{T} <: FieldVector{3, T}
    x::T
    y::T
    z::T
end

unitX(::Type{T} = Float64) where {T<:Real} = Vec3D{T}(one(T), zero(T), zero(T))
unitY(::Type{T} = Float64) where {T<:Real} = Vec3D{T}(zero(T), one(T), zero(T))
unitZ(::Type{T} = Float64) where {T<:Real} = Vec3D{T}(zero(T), zero(T), one(T))
OpticSim.origin(::Type{T} = Float64) where {T<:Real} = Vec3D{T}(zero(T), zero(T), zero(T))
export unitX, unitY, unitZ, origin

function get_orthogonal_vectors(direction::Vec3D{T}) where {T<:Number}
    axis1 = normalize(direction)

    dp = dot(unitX(T), axis1)
    if	( dp >= one(T) - eps(T))
        axis1 = unitX(T);
        axis2 = unitY(T);
        axis3 = unitZ(T);
        return (axis2, axis3)
    elseif  ( dp <= -one(T) + eps(T))
        axis1 = -unitX(T);
        axis2 = -unitY(T);
        axis3 = -unitZ(T);
        return (axis2, axis3)
    end

    axis3 = normalize(cross(axis1, unitX(T)))
    axis2 = normalize(cross(axis3, axis1))
    return (axis2, axis3)
end

#---------------------------------------
# 3D Transform / Local Frame
#---------------------------------------
struct Transform{T} 
    m_origin::Vec3D{T}
    right::Vec3D{T}     # x
    up::Vec3D{T}        # y
    forward::Vec3D{T}   # z
    scale::Vec3D{T}

    # function Transform(::Type{T} = Float64) where {T<:Real}
    #     return new{T}(origin(T), unitX(T), unitY(T), unitZ(T), Vec3D(one(T), one(T), one(T)))
    # end

    function Transform(position::Vec3D{T} = OpticSim.origin(), forward::Vec3D{T} = unitZ(), ::Type{T} = Float64) where {T<:Real}
        forward = normalize(forward)
        right, up = get_orthogonal_vectors(forward)
        return new{T}(position, right, up, forward, Vec3D(one(T), one(T), one(T)))
    end

    function Transform(mat::SMatrix{4, 4, T}) where {T<:Real}
        rot, t, s = decomposeRTS(mat)
        o = t
        r = Vec3D(rot[1:3,1])
        u = Vec3D(rot[1:3,2])
        f = Vec3D(rot[1:3,3])
        return new{T}(o, r, u, f, s)
    end

end


OpticSim.origin(t::Transform{T}) where {T<:Real} = t.m_origin
right(t::Transform{T}) where {T<:Real} = t.right
up(t::Transform{T}) where {T<:Real} = t.up
forward(t::Transform{T}) where {T<:Real} = t.forward
scale(t::Transform{T}) where {T<:Real} = t.scale

function rotation(t::Transform{T}) where {T<:Real}
    return SMatrix{3, 3}(t.right..., t.up..., t.forward...)
end
export rotation

function translation(t::Transform{T}) where {T<:Real}
    return OpticSim.origin(t)
end
export translation

function rotationMat(t::Transform{T}) where {T<:Real}
    return SMatrix{4, 4}(t.right..., 0, t.up..., 0, t.forward..., 0, 0, 0, 0, 1)
end
export rotationMat

function translationMat(t::Transform{T}) where {T<:Real}
    return SMatrix{4, 4}(one(T), zero(T), zero(T), zero(T), zero(T), one(T), zero(T), zero(T), zero(T), zero(T), one(T), zero(T), t.m_origin..., 1)
end
export translationMat

function scaleMat(t::Transform{T}) where {T<:Real}
    return SMatrix{4, 4}(t.scale[1], zero(T), zero(T), zero(T), zero(T), t.scale[2], zero(T), zero(T), zero(T), zero(T), t.scale[3], zero(T), zero(T), zero(T), zero(T), one(T))
end
export scaleMat

function local2worldMat(t::Transform{T}) where {T<:Real}
    return scaleMat(t) * translationMat(t) * rotationMat(t)  
end
export local2worldMat

function decomposeRTS(mat::SMatrix{4, 4, T}) where {T<:Real}
    t = Vec3D{T}(mat[1,4], mat[2, 4], mat[3, 4])
    sx = norm(Vec3D{T}(mat[1,1], mat[2, 1], mat[3, 1]))
    sy = norm(Vec3D{T}(mat[1,2], mat[2, 2], mat[3, 2]))
    sz = norm(Vec3D{T}(mat[1,3], mat[2, 3], mat[3, 3]))
    s = Vec3D{T}(sx, sy, sz)
    rot = SMatrix{4, 4, T}(mat[1,1]/sx, mat[2, 1]/sx, mat[3,1]/sx, 0, mat[1,2]/sy, mat[2, 2]/sy, mat[3,2]/sy, 0, mat[1,3]/sz, mat[2, 3]/sz, mat[3,3]/sz, 0, 0, 0, 0, 1)

    return rot, t, s
end
export decomposeRTS


function Base.:*(x::Transform{T}, y::Transform{T}) where {T<:Real}
    mat = local2worldMat(x) * local2worldMat(y)
    return Transform(mat)
end

function Base.:*(tr::Transform{T}, ray::OpticSim.Ray{T, 3}) where {T<:Real}
    o = local2world(tr, OpticSim.origin(ray))    
    d = rotate(tr, OpticSim.direction(ray))    
    return OpticSim.Ray(o, d)
end

function Base.:*(tr::Transform{T}, ray::OpticSim.OpticalRay{T}) where {T<:Real}
    r = tr * ray.ray
    or = OpticSim.OpticalRay(r, ray.power, ray.wavelength; opl=ray.opl, nhits=ray.nhits, sourcenum=ray.sourcenum, sourcepower=ray.sourcepower)
    return or
    # return OpticSim.OpticalRay(OpticSim.origin(r), OpticSim.direction(r), ray.power, ray.wavelength, ray.opl, ray.nhits, ray.sourcenum, ray.sourcepower)
end



function local2world(t::Transform{T},  v::AbstractVector{T}) where {T<:Real}
    return OpticSim.origin(t) + (v[1] * right(t) + v[2] * up(t) + v[3] * forward(t));
end


function local2worldRBT(t::Transform{T}) where {T<:Real}
    return OpticSim.RigidBodyTransform(rotation(T), OpticSim.origin(t))
end

function rotate(t::Transform{T}, v::AbstractVector{T}) where {T<:Real}
    return rotation(t) * v
end

function rotXmat(angle::T, ::Type{T} = Float64) where {T<:Real} 
    # constructor is column major

    c = cos(angle)
    s = sin(angle)

    return SMatrix{3,3,T,9}(
        one(T), zero(T), zero(T),
        zero(T), c, s,
        zero(T), -s, c
    )
    
end

function rotYmat(angle::T, ::Type{T} = Float64) where {T<:Real} 
    # constructor is column major

    c = cos(angle)
    s = sin(angle)

    return SMatrix{3,3,T,9}(
        c, zero(T), s,
        zero(T), one(T), zero(T),
        -s, zero(T), c
    )
    
end

function rotZmat(angle::T, ::Type{T} = Float64) where {T<:Real} 
    # constructor is column major

    c = cos(angle)
    s = sin(angle)

    return SMatrix{3,3,T,9}(
        c, -s, zero(T),
        s, c, zero(T),
        zero(T), zero(T), one(T)
    )
    
end

function rotAxisMat(axis::SVector{3, T}, angle::T, ::Type{T} = Float64) where {T<:Real} 
    c = cos(-angle);
    s = sin(-angle);
    t = one(T) - c;

    axis = normalize(axis)

    return SMatrix{3,3,T,9}(
        t * axis[1] * axis[1] + c, t * axis[1] * axis[2] + s * axis[3], t * axis[1] * axis[3] - s * axis[2],
        t * axis[1] * axis[2] - s * axis[3], t * axis[2] * axis[2] + c, t * axis[2] * axis[3] + s * axis[1],
        t * axis[1] * axis[3] + s * axis[2], t * axis[2] * axis[3] - s * axis[1], t * axis[3] * axis[3] + c
    )
end


#---------------------------------------
# exports
#---------------------------------------
export  Vec3D, unitX, unitY, unitZ, origin, Transform,
        get_orthogonal_vectors,
        Transform, origin, right, up, forward, scale, rotation, translation, local2world, local2worldMat, local2worldRBT, rotate,
        rotXmat, rotYmat, rotZmat, rotAxisMat

end # module Geometry
