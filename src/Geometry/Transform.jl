# -----------------------------------------------------------------------------------------------
# GEOMETRY
# -----------------------------------------------------------------------------------------------
module Geometry

using StaticArrays
using LinearAlgebra

#region Vec3

"""
    representing a 3D vector
"""
struct Vec3{T} <: FieldVector{3, T}
    _x::T
    _y::T
    _z::T
end
export Vec3

# empty constructor - initialized with zeros
function Vec3(::Type{T} = Float64) where {T<:Real}
    return Vec3{T}(zero(T), zero(T), zero(T))
end

unitX3(::Type{T} = Float64) where {T<:Real} = Vec3{T}(one(T), zero(T), zero(T))
unitY3(::Type{T} = Float64) where {T<:Real} = Vec3{T}(zero(T), one(T), zero(T))
unitZ3(::Type{T} = Float64) where {T<:Real} = Vec3{T}(zero(T), zero(T), one(T))
origin3(::Type{T} = Float64) where {T<:Real} = Vec3{T}(zero(T), zero(T), zero(T))
zero3(::Type{T} = Float64) where {T<:Real} = Vec3{T}(zero(T), zero(T), zero(T))
one3(::Type{T} = Float64) where {T<:Real} = Vec3{T}(one(T), one(T), one(T))
export unitX3, unitY3, unitZ3, origin3, zero3, one3

#endregion Vec3

#region Vec4

"""
    representing a 4D vector
"""
struct Vec4{T} <: FieldVector{4, T}
    _x::T
    _y::T
    _z::T
    _w::T
end
export Vec4

# empty constructor - initialized with zeros
function Vec4(::Type{T} = Float64) where {T<:Real}
    return Vec4{T}(zero(T), zero(T), zero(T), zero(T))
end

# convert vec3 to vec4
function Vec4(v::Vec3{T}) where {T<:Real}
    return Vec4{T}(v..., one(T))
end

function Vec4(v::SVector{3, T}) where {T<:Real}
    return Vec4{T}(v..., one(T))
end

unitX4(::Type{T} = Float64) where {T<:Real} = Vec4{T}(one(T), zero(T), zero(T), zero(T))
unitY4(::Type{T} = Float64) where {T<:Real} = Vec4{T}(zero(T), one(T), zero(T), zero(T))
unitZ4(::Type{T} = Float64) where {T<:Real} = Vec4{T}(zero(T), zero(T), one(T), zero(T))
unitW4(::Type{T} = Float64) where {T<:Real} = Vec4{T}(zero(T), zero(T), zero(T), one(T))
origin4(::Type{T} = Float64) where {T<:Real} = Vec4{T}(zero(T), zero(T), zero(T), zero(T))
zero4(::Type{T} = Float64) where {T<:Real} = Vec4{T}(zero(T), zero(T), zero(T), zero(T))
one4(::Type{T} = Float64) where {T<:Real} = Vec4{T}(one(T), one(T), one(T), one(T))
export unitX4, unitY4, unitZ4, unitW4, origin4, zero4, one4
#endregion Vec4

#region Transform

# utility function that given a vector, return 2 orthogonal vectors to that one 
function get_orthogonal_vectors(direction::Vec3{T}) where {T<:Real}
    axis1 = normalize(direction)

    dp = dot(unitX3(T), axis1)
    # check the special case where the input vector is parallel to the X axis
    if	( dp >= one(T) - eps(T))
        # axis1 = unitX(T);
        axis2 = unitY3(T);
        axis3 = unitZ3(T);
        return (axis2, axis3)
    elseif  ( dp <= -one(T) + eps(T))
        # axis1 = -unitX(T);
        axis2 = -unitY3(T);
        axis3 = -unitZ3(T);
        return (axis2, axis3)
    end

    axis3 = normalize(cross(axis1, unitX3(T)))
    axis2 = normalize(cross(axis3, axis1))
    return (axis2, axis3)
end

#---------------------------------------
# 3D Transform / Local Frame
#---------------------------------------
struct Transform{T} <: FieldMatrix{4, 4, T}
    _xx::T
    _yx::T
    _zx::T
    _wx::T
    _xy::T
    _yy::T
    _zy::T
    _wy::T
    _xz::T
    _yz::T
    _zz::T
    _wz::T
    _xw::T
    _yw::T
    _zw::T
    _ww::T
end
export Transform

"""
    identity([S::Type]) -> Transform{S}

Returns the [`Transform`](@ref) of type `S` (default `Float64`) representing the identity transform.
"""
identityT(::Type{T} = Float64) where {T<:Real} = Transform{T}(
    one(T), zero(T), zero(T), zero(T),
    zero(T), one(T), zero(T), zero(T),
    zero(T), zero(T), one(T), zero(T),
    zero(T), zero(T), zero(T), one(T)
)
export identityT

# for compatability ith the "old" RigidBodyTransform
identitytransform(::Type{T} = Float64) where {T<:Real} = identityT(T)
export identitytransform


"""
    Transform([S::Type]) -> Transform{S}

Returns the [`Transform`](@ref) of type `S` (default `Float64`) representing the identity transform.
"""
function Transform(::Type{T} = Float64) where {T<:Real}
    return identityT(T)
end

"""
    Transform(colx::Vec3{T}, coly::Vec3{T},colz::Vec3{T}, colw::Vec3{T}, ::Type{T} = Float64) where {T<:Real}

Costruct a transform frp, the input columns.     
"""
function Transform(colx::Vec3{T}, coly::Vec3{T},colz::Vec3{T}, colw::Vec3{T} = zero3(T), ::Type{T} = Float64) where {T<:Real}
    return Transform{T}(colx..., zero(T), coly..., zero(T), colz..., zero(T), colw..., one(T))
end

"""
    Transform(colx::Vec3{T}, coly::Vec3{T},colz::Vec3{T}, colw::Vec3{T}, ::Type{T} = Float64) where {T<:Real}

Costruct a transform frp, the input columns.     
"""
function Transform(colx::Vec4{T}, coly::Vec4{T},colz::Vec4{T}, colw::Vec4{T}, ::Type{T} = Float64) where {T<:Real}
    return Transform{T}(colx..., coly..., colz..., colw...)
end

"""
    Transform(origin, forward) -> Transform{S}

Returns the [`Transform`](@ref) of type `S` (default `Float64`) representing the local frame with origin and forward direction. the other 2 axes are computed automaticlly.
"""
function Transform(origin::Vec3{T}, forward::Vec3{T} = unitZ3(); type::Type{T} = Float64) where {T<:Real}
    forward = normalize(forward)
    right, up = get_orthogonal_vectors(forward)
    return Transform(right, up, forward, origin)
end

function Transform{S}(θ::T, ϕ::T, ψ::T, x::T, y::T, z::T; type::Type{S} = Float64) where {T<:Number,S<:Real} 
    temp_transform = Transform(rotmat(S, θ, ϕ, ψ), Vec3{S}(x, y, z))
    return Transform{S}(temp_transform)
end

function Transform(rotation::SMatrix{3,3,T}, translation::SVector{3,T}) where {T<:Real} 
    return Transform(
        rotation[1,1], rotation[2,1], rotation[3,1], zero(T), 
        rotation[1,2], rotation[2,2], rotation[3,2], zero(T), 
        rotation[1,3], rotation[2,3], rotation[3,3], zero(T), 
        translation[1], translation[2], translation[3], one(T))
end

function Transform(rotation::AbstractArray{T,2}, translation::AbstractArray{T,1}) where {T<:Real}
    @assert size(rotation)[1] == size(rotation)[2] == length(translation) == 3
    return Transform(
        rotation[1,1], rotation[2,1], rotation[3,1], zero(T), 
        rotation[1,2], rotation[2,2], rotation[3,2], zero(T), 
        rotation[1,3], rotation[2,3], rotation[3,3], zero(T), 
        translation[1], translation[2], translation[3], one(T))
end


"""
    rotationX(angle::T) where {T<:Real}

Builds a rotation matrix for a rotation around the x-axis. 
Parameters:
    The counter-clockwise `angle` in radians.
"""
function rotationX(angle::T) where {T<:Real}
    c = cos(angle);
    s = sin(angle);

    row1 = unitX4()
    row2 = Vec4(zero(T), c, s, zero(T))
    row3 = Vec4(zero(T), -s, c, zero(T))
    row4 = unitW4()

    # transposing because the constructors treat these vectors as columns instead of rows
    return transpose(Transform(row1, row2, row3, row4))
end
export rotationX

"""
    rotationY(angle::T) where {T<:Real}

Builds a rotation matrix for a rotation around the y-axis. 
Parameters:
    The counter-clockwise `angle` in radians.
"""
function rotationY(angle::T) where {T<:Real}
    c = cos(angle);
    s = sin(angle);

    row1 = Vec4(c, zero(T), -s, zero(T))
    row2 = unitY4()
    row3 = Vec4(s, zero(T), c, zero(T))
    row4 = unitW4()
   
    # transposing because the constructors treat these vectors as columns instead of rows
    return transpose(Transform(row1, row2, row3, row4))
end
export rotationY

"""
    rotationZ(angle::T) where {T<:Real}

Builds a rotation matrix for a rotation around the z-axis. 
Parameters:
    The counter-clockwise `angle` in radians.
"""
function rotationZ(angle::T) where {T<:Real}
    c = cos(angle);
    s = sin(angle);

    row1 = Vec4(c, s, zero(T), zero(T))
    row2 = Vec4(-s, c, zero(T), zero(T))
    row3 = unitZ4()
    row4 = unitW4()
    
    # transposing because the constructors treat these vectors as columns instead of rows
    return transpose(Transform(row1, row2, row3, row4))
end
export rotationZ

function rotation(t::Transform{T}) where {T<:Real}
    rot = SMatrix{3, 3, T}(
        t[1, 1], t[2, 1], t[3, 1], 
        t[1, 2], t[2, 2], t[3, 2],
        t[1, 3], t[2, 3], t[3, 3])
    return Transform(rot, zero3(T))
end

rotate(a::Transform{T}, vector::Union{Vec3{T}, SVector{3,T}}) where {T<:Real} = rotation(a) * vector


"""
    rotation([S::Type], θ::T, ϕ::T, ψ::T) -> Transform{S}

Returns the [`Transform`](@ref) of type `S` (default `Float64`) representing the rotation by `θ`, `ϕ` and `ψ` around the *x*, *y* and *z* axes respectively **in radians**.
"""
rotation(θ::T, ϕ::T, ψ::T) where {T<:Number} = rotation(Float64, θ, ϕ, ψ)
rotation(::Type{S}, θ::T, ϕ::T, ψ::T) where {T<:Number,S<:Real} = Transform(rotmat(S, θ, ϕ, ψ), zeros(SVector{3,S}))
"""
    rotationd([S::Type], θ::T, ϕ::T, ψ::T) -> Transform{S}

Returns the [`Transform`](@ref) of type `S` (default `Float64`) representing the rotation by `θ`, `ϕ` and `ψ` around the *x*, *y* and *z* axes respectively **in degrees**.
"""
rotationd(θ::T, ϕ::T, ψ::T) where {T<:Number} = rotationd(Float64, θ, ϕ, ψ)
rotationd(::Type{S}, θ::T, ϕ::T, ψ::T) where {T<:Number,S<:Real} = Transform(rotmatd(S, θ, ϕ, ψ), zeros(SVector{3,S}))
export translation, rotation, rotationd


"""
    translation(x::T, y::T, z::T) where {T<:Real}

Creates a translation transform
"""
translation(::Type{S}, x::T, y::T, z::T) where {T<:Number,S<:Real} = convert(Transform{S},translation(x, y, z))
function translation(x::T, y::T, z::T) where {T<:Real}
    col1 = unitX4(T)
    col2 = unitY4(T)
    col3 = unitZ4(T)
    col4 = Vec4(x, y, z, one(T))
    
    return Transform(col1, col2, col3, col4, T)
end

"""
    translation(x::T, y::T, z::T) where {T<:Real}

Creates a translation transform
"""
function translation(t::Vec3{T}) where {T<:Real}
    return translation(t...)
end
export translation

"""
    scale(x::T, y::T, z::T) where {T<:Real}

Creates a scaling transform
"""
function scale(x::T, y::T, z::T) where {T<:Real}
    col1 = unitX4(T) * x
    col2 = unitY4(T) * y
    col3 = unitZ4(T) * z
    col4 = unitW4(T)

    return Transform(col1, col2, col3, col4)
end

"""
    scale(s::T) where {T<:Real}

Creates a uniform scaling transform
"""
function scale(s::T) where {T<:Real}
    return scale(s, s, s)
end

"""
    scale(t::Vec3{T}) where {T<:Real}

Creates a scaling transform
"""
function scale(t::Vec3{T}) where {T<:Real}
    return scale(t...)
end
export scale

"""
    local2world(t::Transform{T}) where {T<:Real}

return the transform matrix that takes a point in the local coordinate system to the global one
"""
function local2world(t::Transform{T}) where {T<:Real}
    return t
end
export local2world

"""
    world2local(t::Transform{T}) where {T<:Real}

return the transform matrix that takes a point in the global coordinate system to the local one
"""
function world2local(t::Transform{T}) where {T<:Real}
    return inv(t)
end
export world2local

function Base.:*(t::Transform{T}, v::Vec3{T}) where {T<:Real}
    res = t * Vec4(v)
    if (t[4,4] == one(T))
        return Vec3(res[1], res[2], res[3])
    else    
        return Vec3(res[1]/res[4], res[2]/res[4], res[3]/res[4])
    end
end

function Base.:*(t::Transform{T}, v::SVector{3,T}) where {T<:Real}
    res = t * Vec4(v)
    if (t[4,4] == one(T))
        return SVector(res[1], res[2], res[3])
    else    
        return SVector(res[1]/res[4], res[2]/res[4], res[3]/res[4])
    end
end

# function Base.:*(t::Transform{T}, v::SVector{3,T}) where {T<:Real}
#     res = t * Vec4(v)
#     if (t[4,4] == one(T))
#         return Vec3(res[1], res[2], res[3])
#     else    
#         return Vec3(res[1]/res[4], res[2]/res[4], res[3]/res[4])
#     end
# end

# function Base.:*(t::Transform{T}, v::Vec4{T}) where {T<:Real}
#     res = SMatrix(t) * v
# end


"""
    decomposeRTS(tr::Transform{T}) where {T<:Real}

return a touple containing the rotation matrix, the translation vector and the scale vecto represnting the transform.
"""
function decomposeRTS(tr::Transform{T}) where {T<:Real}
    t = Vec3(tr[1,4], tr[2,4], tr[3,4])
    sx = norm(Vec3(tr[1,1], tr[2,1], tr[3,1]))
    sy = norm(Vec3(tr[1,2], tr[2,2], tr[3,2]))
    sz = norm(Vec3(tr[1,3], tr[2,3], tr[3,3]))
    s = Vec3(sx, sy, sz)
    rot = SMatrix{4, 4, T}(tr[1,1]/sx, tr[2, 1]/sx, tr[3,1]/sx, 0, tr[1,2]/sy, tr[2, 2]/sy, tr[3,2]/sy, 0, tr[1,3]/sz, tr[2, 3]/sz, tr[3,3]/sz, 0, 0, 0, 0, 1)

    return rot, t, s
end
export decomposeRTS


"""
    rotmatbetween([S::Type], a::SVector{3,T}, b::SVector{3,T}) -> SMatrix{3,3,S}

Returns the rotation matrix of type `S` (default `Float64`) representing the rotation between vetors `a` and `b`, i.e. rotation(a,b) * a = b.
"""
rotmatbetween(a::Vec3{T}, b::Vec3{T}) where {T<:Real} = rotmatbetween(Float64, a, b)
function rotmatbetween(::Type{S}, a::Vec3{T}, b::Vec3{T}) where {T<:Real,S<:Real}
    # TODO: Brian, is there a hidden assumption that a and b are normalized?
    v = cross(a, b)
    c = dot(a, b)
    V = SMatrix{3,3,S,9}(0, v[3], -v[2], -v[3], 0, v[1], v[2], -v[1], 0)
    R = I + V + V^2 * one(T) / (one(T) + c)
    return SMatrix{3,3,S,9}(R)
end
rotmatbetween(a::SVector{3,T}, b::SVector{3,T}) where {T<:Real} = rotmatbetween(Float64, Vec3(a), Vec3(b))
function rotmatbetween(type::Type{S}, a::SVector{3,T}, b::SVector{3,T}) where {T<:Real,S<:Real}
    return rotmatbetween(type, Vec3(a), Vec3(b))
end
export rotmatbetween


"""
    rotmatd([S::Type], θ::T, ϕ::T, ψ::T) -> SMatrix{3,3,S}

Returns the rotation matrix of type `S` (default `Float64`) representing the rotation by `θ`, `ϕ` and `ψ` around the *x*, *y* and *z* axes respectively **in degrees**.
"""
rotmatd(θ::T, ϕ::T, ψ::T) where {T<:Number} = rotmat(Float64, deg2rad(θ), deg2rad(ϕ), deg2rad(ψ))
rotmatd(::Type{S}, θ::T, ϕ::T, ψ::T) where {T<:Number,S<:Real} = rotmat(S, deg2rad(θ), deg2rad(ϕ), deg2rad(ψ))
export rotmatd
"""
    rotmat([S::Type], θ::T, ϕ::T, ψ::T) -> SMatrix{3,3,S}

Returns the rotation matrix of type `S` (default `Float64`) representing the rotation by `θ`, `ϕ` and `ψ` around the *x*, *y* and *z* axes respectively **in radians**.
"""
rotmat(θ::T, ϕ::T, ψ::T) where {T<:Number} = rotmat(Float64, θ, ϕ, ψ)
function rotmat(::Type{S}, θ::T, ϕ::T, ψ::T) where {T<:Number,S<:Real}
    sinψ = sin(ψ)
    sinϕ = sin(ϕ)
    sinθ = sin(θ)
    cosψ = cos(ψ)
    cosϕ = cos(ϕ)
    cosθ = cos(θ)
    return SMatrix{3,3,S,9}(cosψ * cosϕ, sinψ * cosϕ, -sinϕ, cosψ * sinϕ * sinθ - sinψ * cosθ, sinψ * sinϕ * sinθ + cosψ * cosθ, cosϕ * sinθ, cosψ * sinϕ * cosθ + sinψ * sinθ, sinψ * sinϕ * cosθ - cosψ * sinθ, cosϕ * cosθ)
end
export rotmat



# the following functions should be moved to their appropriate location after we arrange the modules correctly.
# handling of specific usages of Transform needs to be near the specific types that are handled.

# function Base.:*(a::Transform{T}, r::Ray{T,3})::Ray{T,3} where {T}
#     return Ray(a * origin(r), rotate(a, direction(r)))
# end

# function Base.:*(a::RigidBodyTransform{T}, intsct::Intersection{T,3})::Intersection{T,3} where {T<:Real}
#     u, v = uv(intsct)
#     i = interface(intsct)
#     if VERSION < v"1.6.0-DEV"
#         # TODO REMOVE
#         return @unionsplit OpticalInterface T i Intersection(α(intsct), a * point(intsct), rotate(a, normal(intsct)), u, v, i, flippednormal = flippednormal(intsct))
#     else
#         return Intersection(α(intsct), a * point(intsct), rotate(a, normal(intsct)), u, v, interface(intsct), flippednormal = flippednormal(intsct))
#     end
# end

# function Base.:*(transformation::RigidBodyTransform{T}, a::Interval{T}) where {T<:Real}
#     # looks ridiculous but necessary to dissambiguate the elements of the interval
#     u = upper(a)
#     l = lower(a)
#     if l isa RayOrigin{T}
#         if u isa Infinity{T}
#             return Interval(l, u)
#         else
#             u = transformation * u
#             return Interval(l, u)
#         end
#     else
#         l = transformation * l
#         if u isa Infinity{T}
#             return Interval(l, u)
#         else
#             u = transformation * u
#             return Interval(l, u)
#         end
#     end
# end

# function Base.:*(a::RigidBodyTransform{T}, tmesh::TriangleMesh{T})::TriangleMesh{T} where {T<:Real}
#     newT = Vector{Triangle{T}}(undef, length(tmesh.triangles))
#     @inbounds @simd for i in 1:length(tmesh.triangles)
#         newT[i] = a * tmesh.triangles[i]
#     end
#     return TriangleMesh(newT)
# end

# Base.:*(a::RigidBodyTransform{T}, t::Triangle{T}) where {T<:Real} = Triangle(a * vertex(t, 1), a * vertex(t, 2), a * vertex(t, 3))

# function Base.:*(a::RigidBodyTransform{T}, vector::SVector{3,T}) where {T<:Real}
#     if a.rotation === SMatrix{3,3,T,9}(I)
#         return vector + a.translation
#     else
#         # in julia 0 * Inf = NaN, so if the vector has any infinite value in it then the rotation breaks
#         # have to manually test and correct this - this can still fail when summing Infs of different signs (e.g. Inf + -Inf = NaN)
#         # this happens when the rotation isn't axis aligned so the result becomes less meaningful anyway
#         # in this case we fall back to infinite bounding boxes anyway which are usually clipped later in the csg
#         # TODO better way to do this?
#         if any(isinf.(vector))
#             x = zero(T)
#             y = zero(T)
#             z = zero(T)
#             @inbounds begin
#                 x += abs(a.rotation[1, 1]) <= eps(T) ? zero(T) : a.rotation[1, 1] * vector[1]
#                 x += abs(a.rotation[1, 2]) <= eps(T) ? zero(T) : a.rotation[1, 2] * vector[2]
#                 x += abs(a.rotation[1, 3]) <= eps(T) ? zero(T) : a.rotation[1, 3] * vector[3]
#                 y += abs(a.rotation[2, 1]) <= eps(T) ? zero(T) : a.rotation[2, 1] * vector[1]
#                 y += abs(a.rotation[2, 2]) <= eps(T) ? zero(T) : a.rotation[2, 2] * vector[2]
#                 y += abs(a.rotation[2, 3]) <= eps(T) ? zero(T) : a.rotation[2, 3] * vector[3]
#                 z += abs(a.rotation[3, 1]) <= eps(T) ? zero(T) : a.rotation[3, 1] * vector[1]
#                 z += abs(a.rotation[3, 2]) <= eps(T) ? zero(T) : a.rotation[3, 2] * vector[2]
#                 z += abs(a.rotation[3, 3]) <= eps(T) ? zero(T) : a.rotation[3, 3] * vector[3]
#             end
#             return SVector(x, y, z) + a.translation
#         else
#             return a.rotation * vector + a.translation
#         end
#     end
# end
# Base.:*(a::RigidBodyTransform{T}, b::RigidBodyTransform{T}) where {T<:Real} = RigidBodyTransform(a.rotation * b.rotation, a.translation + a.rotation * b.translation)
# rotate(a::RigidBodyTransform{T}, vector::SVector{3,T}) where {T<:Real} = a.rotation * vector
# Base.inv(a::RigidBodyTransform{T}) where {T<:Real} = RigidBodyTransform(a.rotation', -a.rotation' * a.translation)
# Base.print(io::IO, a::RigidBodyTransform) = print(io, hcat(a.rotation, a.translation))
# Base.collect(a::RigidBodyTransform) = hcat(a.rotation, a.translation)



#endregion Transform

end # module Geometry



