# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

# -----------------------------------------------------------------------------------------------
# GEOMETRY
# -----------------------------------------------------------------------------------------------
using ...OpticSim
using StaticArrays
using LinearAlgebra

#region Vec3

"""
`Vec3{T}` provides an immutable vector of fixed length 3 and type `T`.
    
`Vec3` defines a series of convenience constructors, so you can just type e.g. `Vec3(1, 2, 3)` or `Vec3([1.0, 2.0, 3.0])`. 
It also supports comprehensions, and the `zeros()`, `ones()`, `fill()`, `rand()` and `randn()` functions, such as `Vec3(rand(3))`.
"""
Vec3{T} = SVector{3,T}
export Vec3

# empty constructor - initialized with zeros
Vec3(::Type{T} = Float64) where {T<:Real} = zeros(Vec3{T})

"""
returns the unit vector `[1, 0, 0]`
"""
unitX3(::Type{T} = Float64) where {T<:Real} = Vec3{T}(one(T), zero(T), zero(T))
"""
returns the unit vector `[0, 1, 0]`
"""
unitY3(::Type{T} = Float64) where {T<:Real} = Vec3{T}(zero(T), one(T), zero(T))
"""
returns the unit vector `[0, 0, 1]`
"""
unitZ3(::Type{T} = Float64) where {T<:Real} = Vec3{T}(zero(T), zero(T), one(T))

export unitX3, unitY3, unitZ3

#endregion Vec3

#region Vec4

"""
`Vec4{T}` provides an immutable vector of fixed length 4 and type `T`.
    
`Vec4` defines a series of convenience constructors, so you can just type e.g. `Vec3(1, 2, 3, 4)` or `Vec3([1.0, 2.0, 3.0, 4.0])`. 
It also supports comprehensions, and the `zeros()`, `ones()`, `fill()`, `rand()` and `randn()` functions, such as `Vec4(rand(4))`.
"""
Vec4{T} = SVector{4,T}

# empty constructor - initialized with zeros
Vec4(::Type{T} = Float64) where {T<:Real} = zeros(Vec4{T})
export Vec4

"""
    Vec4(v::SVector{3, T}) where {T<:Real} -> Vec4{T}

Accept `SVector` and create a `Vec4` type [v[1], v[2], v[3], 1]
"""
function Vec4(v::SVector{3, T}) where {T<:Real}
    return Vec4{T}(v[1], v[2], v[3], one(T))
end

""" 
    Vec4(m::SMatrix{3,N,T} where{N,T<:Real} -> SMatrix{3,N,T})
Input is matrix of 3d points, each column is one point. Returns matrix of 3d points with 1 appended in the last row.
"""
function Vec4(m::SMatrix{3,N,T}) where {N,T<:Real}
    return vcat(m,ones(SMatrix{1,N,T}))::SMatrix{4,N,T}
end
"""
returns the unit vector `[1, 0, 0, 0]`
"""
unitX4(::Type{T} = Float64) where {T<:Real} = Vec4{T}(one(T), zero(T), zero(T), zero(T))
"""
returns the unit vector `[0, 1, 0, 0]`
"""
unitY4(::Type{T} = Float64) where {T<:Real} = Vec4{T}(zero(T), one(T), zero(T), zero(T))
"""
returns the unit vector `[0, 0, 1, 0]`
"""
unitZ4(::Type{T} = Float64) where {T<:Real} = Vec4{T}(zero(T), zero(T), one(T), zero(T))
"""
returns the unit vector `[0, 0, 0, 1]`
"""
unitW4(::Type{T} = Float64) where {T<:Real} = Vec4{T}(zero(T), zero(T), zero(T), one(T))

export unitX4, unitY4, unitZ4, unitW4
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
export get_orthogonal_vectors

#---------------------------------------
# 3D Transform / Local Frame
#---------------------------------------
"""
    Transform{S<:Real}

Transform encapsulating rotation, translation and scale in 3D space. Translation happens **after** rotation.

```julia
Transform{S}(θ::T, ϕ::T, ψ::T, x::T, y::T, z::T)
Transform(rotation::SMatrix{3,3,S}, translation::SVector{3,S})
Transform(rotation::AbstractArray{S,2}, translation::AbstractArray{S,1})
```
`θ`, `ϕ` and `ψ` in first constructor are in **radians**.
"""
struct Transform{T} 
    matrix::SMatrix{4,4,T,16}

    """ This is a private internal function. In general don't want to allow users to populate Transform matrices with arbitrary elements. Not to be called by code outside of the Transform module. Don't use Transform{Float64}(...) for example. Instead use Transform(..)"""
    function Transform{T}(a11::T,a21::T,a31::T,a41::T,a12::T,a22::T,a32::T,a42::T,a13::T,a23::T,a33::T,a43::T,a14::T,a24::T,a34::T,a44::T) where{T<:Real} 
        return new{T}(SMatrix{4,4,T,16}(a11,a21,a31,a41,a12,a22,a32,a42,a13,a23,a33,a43,a14,a24,a34,a44))
    end

    Transform{T}(mat::SMatrix{4,4,T,16}) where{T<:Real} = new{T}(mat)
end
export Transform

#functions to make Transform compatible with base matrix API
matrix(a::Transform) = a.matrix

# Base.length(a::Transform) = length(matrix(a))
Base.getindex(a::Transform, indices::Vararg{Int,N}) where{N} = getindex(matrix(a),indices...)
# Base.iterate(a::Transform) = iterate(matrix(a))
# Base.iterate(a::Transform{Float64}, b::Tuple{StaticArrays.SOneTo{16}, Int64}) = iterate(matrix(a),b)

Base.collect(a::Transform) = collect(matrix(a))

Base.:*(transa::Transform{T},transb::Transform{T}) where{T<:Real} = Transform{T}(matrix(transa)*matrix(transb))

function Base.:*(transform::Transform{T}, v::SVector{3,S}) where {T<:Real,S<:Number}
    t = matrix(transform)

    res = t * Vec4(v)
    if (t[4,4] == one(T))
        return SVector{3,S}(res[1], res[2], res[3])
    else    
        return SVector{3,S}(res[1]/res[4], res[2]/res[4], res[3]/res[4])
    end
end

#Transform is not necessarily constrained to be a rigid body transformation so use general invers.
Base.inv(a::Transform{T}) where{T<:Real}= Transform{T}(inv(matrix(a)))

""" The t and m matrices are allowed to be of different element type. This allows transforming a Unitful matrix for example:
```
id = identitytransform()
m = fill(1mm,3,4)
id*m #returns a matrix filled with Unitful quantities. If both matrices had to be the same type this would not work
```
"""
function Base.:*(transform::Transform{T}, m::SMatrix{3,N,S}) where{N,T<:Real,S<:Number}
    res = MMatrix{3,N,T}(undef)
    t = matrix(transform)

    for outcol in 1:N
        for row in 1:3
            sum = T(0)
            for incol in 1:3
                sum += t[row,incol]*m[incol,outcol]
            end
            #implicit 1 w coordinate value
            sum += t[row,4]
            res[row,outcol] = sum
        end
        if t[4,4] != 1
            res[:,outcol] /= t[4,4]
        end
    end
    return SMatrix{3,N,S}(res)
end

""" The t and m matrices are allowed to be of different element type. This allows transforming a Unitful matrix for example:
WARNING: this doesn't work. The translation component of the transform matrix has to be in Unitful units but the rotation part has to be in unitless units for this to work. Only works if one assumes that the translation part of the transform implicitly has the same units as the Unitful vectors being transformed. Brittle and likely to cause obscure bugs.
```
id = identitytransform()
m = fill(1mm,3,4)
id*m #returns a matrix filled with Unitful quantities. If both matrices had to be the same type this would not work
```
"""
Base.:*(transform::Transform{T},m::SMatrix{4,N,S}) where{N,T<:Real,S<:Number} = matrix(transform)*m

Base.transpose(a::Transform{T}) where{T<:Real} = Transform{T}(a.matrix')

# END of functions for compatibility with base matrix API


# for compatibility ith the "old" RigidBodyTransform
"""
identitytransform([S::Type]) -> Transform{S}

Returns the [`Transform`](@ref) of type `S` (default `Float64`) representing the identity transform.
"""
identitytransform(::Type{T} = Float64) where {T<:Real} = Transform{T}(
    one(T), zero(T), zero(T), zero(T),
    zero(T), one(T), zero(T), zero(T),
    zero(T), zero(T), one(T), zero(T),
    zero(T), zero(T), zero(T), one(T)
)
export identitytransform

"""
    Transform([S::Type]) -> Transform{S}

Returns the [`Transform`](@ref) of type `S` (default `Float64`) representing the identity transform.
"""
function Transform(::Type{T} = Float64) where {T<:Real}
    return identitytransform(T)
end

"""
    Transform(colx::Vec3{T}, coly::Vec3{T},colz::Vec3{T}, colw::Vec3{T}, ::Type{T} = Float64) where {T<:Real}

Construct a transform from the input columns.     
"""
function Transform(colx::Vec3{T}, coly::Vec3{T}, colz::Vec3{T}, colw::Vec3{T} = zero(Vec3{T})) where {T<:Real}
    return Transform{T}(vcat(hcat(colx,coly,colz,colw),SMatrix{1,4,T}(zero(T),zero(T),zero(T),one(T)) ))
end


"""
    Transform(colx::Vec3{T}, coly::Vec3{T},colz::Vec3{T}, colw::Vec3{T}, ::Type{T} = Float64) where {T<:Real}

Construct a transform from the input columns.     
"""
function Transform(colx::Vec4{T}, coly::Vec4{T}, colz::Vec4{T}, colw::Vec4{T}) where {T<:Real}
    return Transform{T}(hcat(colx,coly,colz,colw))
end

"""
    Transform(origin, forward) -> Transform{S}

Returns the [`Transform`](@ref) of type `S` (default `Float64`) representing the local frame with origin and forward direction. the other 2 axes are computed automatically.
"""
function Transform(origin::Vec3{T}, forward::Vec3{T} = unitZ3()) where {T<:Real}
    forward = normalize(forward)
    right, up = get_orthogonal_vectors(forward)
    return Transform(right, up, forward, origin)
end

function Transform(θ::T, ϕ::T, ψ::T, x::T, y::T, z::T) where {T<:Number} 
    return Transform(rotmat(T, θ, ϕ, ψ), Vec3{T}(x, y, z))
end

"""
    Transform(rotation::SMatrix{3,3,T}, translation::SVector{3,T}) where {T<:Real} -> Transform{S}

Returns the [`Transform`](@ref) of type `S` (default `Float64`) created by a rotation matrix and translation vector.
"""
function Transform(rotation::SMatrix{3,3,T}, translation::SVector{3,T}) where {T<:Real} 
    return Transform{T}(
        rotation[1,1], rotation[2,1], rotation[3,1], zero(T), 
        rotation[1,2], rotation[2,2], rotation[3,2], zero(T), 
        rotation[1,3], rotation[2,3], rotation[3,3], zero(T), 
        translation[1], translation[2], translation[3], one(T))
end

"""
    Transform(rotation::AbstractArray{T,2}, translation::AbstractArray{T,1}) where {T<:Real} -> Transform{S}

Returns the [`Transform`](@ref) of type `S` (default `Float64`) created by a rotation matrix (3x3) and translation vector of length 3.
"""
function Transform(rotation::AbstractArray{T,2}, translation::AbstractArray{T,1}) where {T<:Real}
    @assert size(rotation)[1] == size(rotation)[2] == length(translation) == 3
    return Transform{T}(
        rotation[1,1], rotation[2,1], rotation[3,1], zero(T), 
        rotation[1,2], rotation[2,2], rotation[3,2], zero(T), 
        rotation[1,3], rotation[2,3], rotation[3,3], zero(T), 
        translation[1], translation[2], translation[3], one(T))
end


# define some utility functions 

"""
    right(t::Transform{<:Real}) -> Vec3

Assuming t is a 3D rigid transform representing a local left-handed coordinate system, this function will return the first column, representing the "X" axis.
"""
right(t::Transform{<:Real}) = normalize(Vec3(t[1,1], t[2,1], t[3,1]))

"""
    up(t::Transform{<:Real}) -> Vec3

Assuming t is a 3D rigid transform representing a local left-handed coordinate system, this function will return the second column, representing the "Y" axis.
"""
up(t::Transform{<:Real}) = normalize(Vec3(t[1,2], t[2,2], t[3,2]))

"""
    forward(t::Transform{<:Real}) -> Vec3

Assuming t is a 3D rigid transform representing a local left-handed coordinate system, this function will return the third column, representing the "Z" axis.
"""
forward(t::Transform{<:Real}) = normalize(Vec3(t[1,3], t[2,3], t[3,3]))

"""
    origin(t::Transform{<:Real}) -> Vec3

Assuming t is a 3D rigid transform representing a local left-handed coordinate system, this function will return the fourth column, containing the translation part of the transform in 3D.
"""
OpticSim.origin(t::Transform{<:Real}) = Vec3(t[1,4], t[2,4], t[3,4])

export right, up, forward, origin


"""
    rotationX(angle::T) where {T<:Real} -> Transform

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
    rotationY(angle::T) where {T<:Real} -> Transform

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
    rotationZ(angle::T) where {T<:Real} -> Transform

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

"""
    rotation(t::Transform{T}) where {T<:Real} -> SMatrix{3,3,T}

returns the rotation part of the transform `t` - a 3x3 matrix.
"""
function rotation(t::Transform{T}) where {T<:Real}
    rot = SMatrix{3, 3, T}(
        t[1, 1], t[2, 1], t[3, 1], 
        t[1, 2], t[2, 2], t[3, 2],
        t[1, 3], t[2, 3], t[3, 3])
    return Transform(rot, zero(Vec3{T}))
end

"""
    rotate(a::Transform{T}, vector::Union{Vec3{T}, SVector{3,T}}) where {T<:Real} -> Vec3{T}

apply the rotation part of the transform `a` to the vector `vector` - this operation is usually used to rotate direction vectors.
"""
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
    
    return Transform(col1, col2, col3, col4)
end

"""
    translation(x::T, y::T, z::T) where {T<:Real}

Creates a translation transform
"""
function translation(t::Vec3{T}) where {T<:Real}
    return translation(t[1], t[2], t[3])
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
    return scale(t[1], [2], [3])
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


    
"""
    decomposeRTS(tr::Transform{T}) where {T<:Real}

return a tuple containing the rotation matrix, the translation vector and the scale vecto represnting the transform.
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
rotmatbetween(a::SVector{3,T}, b::SVector{3,T}) where {T<:Real} = rotmatbetween(Float64, a, b)
function rotmatbetween(::Type{S}, a::SVector{3,T}, b::SVector{3,T}) where {T<:Real,S<:Real}
    # TODO: Brian, is there a hidden assumption that a and b are normalized?
    v = cross(a, b)
    c = dot(a, b)
    V = SMatrix{3,3,S,9}(0, v[3], -v[2], -v[3], 0, v[1], v[2], -v[1], 0)
    R = I + V + V^2 * one(T) / (one(T) + c)
    return SMatrix{3,3,S,9}(R)
end

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


#endregion Transform




