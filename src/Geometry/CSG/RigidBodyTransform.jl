"""
    RigidBodyTransform{S<:Real}

Transform encapsulating rotation and translation in 3D space. Translation happens **after** rotation.

```julia
RigidBodyTransform{S}(θ::T, ϕ::T, ψ::T, x::T, y::T, z::T)
RigidBodyTransform(rotation::SMatrix{3,3,S}, translation::SVector{3,S})
RigidBodyTransform(rotation::AbstractArray{S,2}, translation::AbstractArray{S,1})
```
`θ`, `ϕ` and `ψ` in first constructor are in **radians**.
"""
struct RigidBodyTransform{T<:Real}
    rotation::SMatrix{3,3,T,9}
    translation::SVector{3,T}

    RigidBodyTransform{S}(θ::T, ϕ::T, ψ::T, x::T, y::T, z::T) where {T<:Number,S<:Real} = new{S}(rotmat(S, θ, ϕ, ψ), SVector{3,S}(x, y, z))
    RigidBodyTransform(rotation::SMatrix{3,3,T}, translation::SVector{3,T}) where {T<:Real} = new{T}(rotation, translation)
    function RigidBodyTransform(rotation::AbstractArray{T,2}, translation::AbstractArray{T,1}) where {T<:Real}
        @assert size(rotation)[1] == size(rotation)[2] == length(translation) == 3
        return new{T}(SMatrix{3,3,T}(rotation), SVector{3,T}(translation))
    end
end
export RigidBodyTransform

"""
    identitytransform([S::Type]) -> RigidBodyTransform{S}

Returns the [`RigidBodyTransform`](@ref) of type `S` (default `Float64`) representing the identity transform.
"""
identitytransform(::Type{T} = Float64) where {T<:Real} = RigidBodyTransform(SMatrix{3,3,T,9}(I), zeros(SVector{3,T}))
export identitytransform

"""
    rotmatbetween([S::Type], a::SVector{3,T}, b::SVector{3,T}) -> SMatrix{3,3,S}

Returns the rotation matrix of type `S` (default `Float64`) representing the rotation between vetors `a` and `b`.
"""
rotmatbetween(a::SVector{3,T}, b::SVector{3,T}) where {T<:Number} = rotmatbetween(Float64, a, b)
function rotmatbetween(::Type{S}, a::SVector{3,T}, b::SVector{3,T}) where {T<:Number,S<:Real}
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
rotmatd(θ::T, ϕ::T, ψ::T) where {T<:Number} = rotmat(Float64, θ * π / 180, ϕ * π / 180, ψ * π / 180)
rotmatd(::Type{S}, θ::T, ϕ::T, ψ::T) where {T<:Number,S<:Real} = rotmat(S, θ * π / 180, ϕ * π / 180, ψ * π / 180)
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
"""
    translation([S::Type], x::T, y::T, z::T) -> RigidBodyTransform{S}

Returns the [`RigidBodyTransform`](@ref) of type `S` (default `Float64`) representing the translation by `x`, `y` and `z`.
"""
translation(x::T, y::T, z::T) where {T<:Number} = translation(Float64, x, y, z)
translation(::Type{S}, x::T, y::T, z::T) where {T<:Number,S<:Real} = RigidBodyTransform(SMatrix{3,3,S,9}(I), SVector{3,S}(x, y, z))
"""
    rotation([S::Type], θ::T, ϕ::T, ψ::T) -> RigidBodyTransform{S}

Returns the [`RigidBodyTransform`](@ref) of type `S` (default `Float64`) representing the rotation by `θ`, `ϕ` and `ψ` around the *x*, *y* and *z* axes respectively **in radians**.
"""
rotation(θ::T, ϕ::T, ψ::T) where {T<:Number} = rotation(Float64, θ, ϕ, ψ)
rotation(::Type{S}, θ::T, ϕ::T, ψ::T) where {T<:Number,S<:Real} = RigidBodyTransform(rotmat(S, θ, ϕ, ψ), zeros(SVector{3,S}))
"""
    rotationd([S::Type], θ::T, ϕ::T, ψ::T) -> RigidBodyTransform{S}

Returns the [`RigidBodyTransform`](@ref) of type `S` (default `Float64`) representing the rotation by `θ`, `ϕ` and `ψ` around the *x*, *y* and *z* axes respectively **in degrees**.
"""
rotationd(θ::T, ϕ::T, ψ::T) where {T<:Number} = rotationd(Float64, θ, ϕ, ψ)
rotationd(::Type{S}, θ::T, ϕ::T, ψ::T) where {T<:Number,S<:Real} = RigidBodyTransform(rotmatd(S, θ, ϕ, ψ), zeros(SVector{3,S}))
export translation, rotation, rotationd

function Base.:*(a::RigidBodyTransform{T}, r::Ray{T,3})::Ray{T,3} where {T}
    return Ray(a * origin(r), rotate(a, direction(r)))
end

function Base.:*(a::RigidBodyTransform{T}, intsct::Intersection{T,3})::Intersection{T,3} where {T<:Real}
    u, v = uv(intsct)
    i = interface(intsct)
    if VERSION < v"1.6.0-DEV"
        # TODO REMOVE
        return @unionsplit OpticalInterface T i Intersection(α(intsct), a * point(intsct), rotate(a, normal(intsct)), u, v, i, flippednormal = flippednormal(intsct))
    else
        return Intersection(α(intsct), a * point(intsct), rotate(a, normal(intsct)), u, v, interface(intsct), flippednormal = flippednormal(intsct))
    end
end

function Base.:*(transformation::RigidBodyTransform{T}, a::Interval{T}) where {T<:Real}
    # looks ridiculous but necessary to dissambiguate the elements of the interval
    u = upper(a)
    l = lower(a)
    if l isa RayOrigin{T}
        if u isa Infinity{T}
            return Interval(l, u)
        else
            u = transformation * u
            return Interval(l, u)
        end
    else
        l = transformation * l
        if u isa Infinity{T}
            return Interval(l, u)
        else
            u = transformation * u
            return Interval(l, u)
        end
    end
end

function Base.:*(a::RigidBodyTransform{T}, tmesh::TriangleMesh{T})::TriangleMesh{T} where {T<:Real}
    newT = Vector{Triangle{T}}(undef, length(tmesh.triangles))
    @inbounds @simd for i in 1:length(tmesh.triangles)
        newT[i] = a * tmesh.triangles[i]
    end
    return TriangleMesh(newT)
end

Base.:*(a::RigidBodyTransform{T}, t::Triangle{T}) where {T<:Real} = Triangle(a * vertex(t, 1), a * vertex(t, 2), a * vertex(t, 3))

function Base.:*(a::RigidBodyTransform{T}, vector::SVector{3,T}) where {T<:Real}
    if a.rotation === SMatrix{3,3,T,9}(I)
        return vector + a.translation
    else
        # in julia 0 * Inf = NaN, so if the vector has any infinite value in it then the rotation breaks
        # have to manually test and correct this - this can still fail when summing Infs of different signs (e.g. Inf + -Inf = NaN)
        # this happens when the rotation isn't axis aligned so the result becomes less meaningful anyway
        # in this case we fall back to infinite bounding boxes anyway which are usually clipped later in the csg
        # TODO better way to do this?
        if any(isinf.(vector))
            x = zero(T)
            y = zero(T)
            z = zero(T)
            @inbounds begin
                x += abs(a.rotation[1, 1]) <= eps(T) ? zero(T) : a.rotation[1, 1] * vector[1]
                x += abs(a.rotation[1, 2]) <= eps(T) ? zero(T) : a.rotation[1, 2] * vector[2]
                x += abs(a.rotation[1, 3]) <= eps(T) ? zero(T) : a.rotation[1, 3] * vector[3]
                y += abs(a.rotation[2, 1]) <= eps(T) ? zero(T) : a.rotation[2, 1] * vector[1]
                y += abs(a.rotation[2, 2]) <= eps(T) ? zero(T) : a.rotation[2, 2] * vector[2]
                y += abs(a.rotation[2, 3]) <= eps(T) ? zero(T) : a.rotation[2, 3] * vector[3]
                z += abs(a.rotation[3, 1]) <= eps(T) ? zero(T) : a.rotation[3, 1] * vector[1]
                z += abs(a.rotation[3, 2]) <= eps(T) ? zero(T) : a.rotation[3, 2] * vector[2]
                z += abs(a.rotation[3, 3]) <= eps(T) ? zero(T) : a.rotation[3, 3] * vector[3]
            end
            return SVector(x, y, z) + a.translation
        else
            return a.rotation * vector + a.translation
        end
    end
end
Base.:*(a::RigidBodyTransform{T}, b::RigidBodyTransform{T}) where {T<:Real} = RigidBodyTransform(a.rotation * b.rotation, a.translation + a.rotation * b.translation)
rotate(a::RigidBodyTransform{T}, vector::SVector{3,T}) where {T<:Real} = a.rotation * vector
Base.inv(a::RigidBodyTransform{T}) where {T<:Real} = RigidBodyTransform(a.rotation', -a.rotation' * a.translation)
Base.print(io::IO, a::RigidBodyTransform) = print(io, hcat(a.rotation, a.translation))
Base.collect(a::RigidBodyTransform) = hcat(a.rotation, a.translation)
