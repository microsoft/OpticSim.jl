
# the following functions should be moved to their appropriate location after we arrange the modules correctly.
# handling of specific usages of Transform needs to be near the specific types that are handled.

function Base.:*(a::Transform{T}, r::Ray{T,3})::Ray{T,3} where {T}
    return Ray(a * origin(r), Geometry.rotate(a, direction(r)))
end

function Base.:*(a::Transform{T}, intsct::Intersection{T,3})::Intersection{T,3} where {T<:Real}
    u, v = uv(intsct)
    i = interface(intsct)
    if VERSION < v"1.6.0-DEV"
        # TODO REMOVE
        return @unionsplit OpticalInterface T i Intersection(α(intsct), a * point(intsct), Geometry.rotate(a, normal(intsct)), u, v, i, flippednormal = flippednormal(intsct))
    else
        return Intersection(α(intsct), a * point(intsct), Geometry.rotate(a, normal(intsct)), u, v, interface(intsct), flippednormal = flippednormal(intsct))
    end
end

function Base.:*(transformation::Transform{T}, a::Interval{T}) where {T<:Real}
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

function Base.:*(a::Transform{T}, tmesh::TriangleMesh{T})::TriangleMesh{T} where {T<:Real}
    newT = Vector{Triangle{T}}(undef, length(tmesh.triangles))
    @inbounds @simd for i in 1:length(tmesh.triangles)
        newT[i] = a * tmesh.triangles[i]
    end
    return TriangleMesh(newT)
end

Base.:*(a::Transform{T}, t::Triangle{T}) where {T<:Real} = Triangle(a * vertex(t, 1), a * vertex(t, 2), a * vertex(t, 3))

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






