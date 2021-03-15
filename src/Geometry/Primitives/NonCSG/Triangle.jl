"""
    Triangle{T} <: Surface{T}

Triangular surface, not a valid CSG object.
Primarily used as a component part of [`TriangleMesh`](@ref) or to enable intersection of [`AcceleratedParametricSurface`](@ref)s.
Can never be used directly as an optical surface as it doesn't have an [`OpticalInterface`](@ref).

```julia
Triangle(v1::SVector{3,T}, v2::SVector{3,T}, v3::SVector{3,T}, [uv1::SVector{2,T}, uv2::SVector{2,T}, uv3::SVector{2,T}])
```
"""
struct Triangle{T} <: Surface{T}
    A::SVector{3,T}
    B::SVector{3,T}
    C::SVector{3,T}
    uv1::SVector{2,T}
    uv2::SVector{2,T}
    uv3::SVector{2,T}
    BA::SVector{3,T}
    CA::SVector{3,T}
    normal::SVector{3,T}

    function Triangle(A::SVector{3,T}, B::SVector{3,T}, C::SVector{3,T}, uv1::SVector{2,T}, uv2::SVector{2,T}, uv3::SVector{2,T}) where {T<:Real}
        BA = B - A
        CA = C - A
        crs = cross(BA, CA)
        lcrs = norm(crs)
        # make sure the triangle is valid
        @assert !samepoint(lcrs, zero(T))
        n̂ = crs / lcrs
        return new{T}(A, B, C, uv1, uv2, uv3, BA, CA, n̂)
    end

    function Triangle(v1::AbstractArray{T,1}, v2::AbstractArray{T,1}, v3::AbstractArray{T,1}, uv1::AbstractArray{T,1}, uv2::AbstractArray{T,1}, uv3::AbstractArray{T,1}) where {T<:Real}
        @assert length(v1) == length(v2) == length(v3) == 3
        @assert length(uv1) == length(uv2) == length(uv3) == 2
        return Triangle(SVector{3}(v1), SVector{3}(v2), SVector{3}(v3), SVector{2}(uv1), SVector{2}(uv2), SVector{2}(uv3))
    end

end
export Triangle

Base.show(io::IO, a::Triangle{T}) where {T<:Real} = print(io, "Triangle{$T}($(a.A), $(a.B), $(a.C))")

function Triangle(v1::AbstractArray{T,1}, v2::AbstractArray{T,1}, v3::AbstractArray{T,1}) where {T<:Real}
    return Triangle(v1, v2, v3, zeros(T, 2), zeros(T, 2), zeros(T, 2))
end

function Triangle(v1::SVector{3,T}, v2::SVector{3,T}, v3::SVector{3,T}) where {T<:Real}
    return Triangle(v1, v2, v3, zeros(SVector{2,T}), zeros(SVector{2,T}), zeros(SVector{2,T}))
end

centroid(tri::Triangle{T}) where {T<:Real} = (vertex(tri, 1) + vertex(tri, 2) + vertex(tri, 3)) / 3
validtri(v1::SVector{3,T}, v2::SVector{3,T}, v3::SVector{3,T}) where {T<:Real} = !any(isnan.(v1)) && !any(isnan.(v2)) && !any(isnan.(v3)) && !samepoint(norm(cross(v2 - v1, v3 - v1)), zero(T))
vertex(tri::Triangle{T}, i::Int) where {T<:Real} = i == 1 ? tri.A : i == 2 ? tri.B : i == 3 ? tri.C : throw(ErrorException("Invalid index: $i"))
vertices(tri::Triangle{T}) where {T<:Real} = (tri.A, tri.B, tri.C)
uvvertex(tri::Triangle{T}, i::Int) where {T<:Real} = i == 1 ? tri.uv1 : i == 2 ? tri.uv2 : i == 3 ? tri.uv3 : throw(ErrorException("Invalid index: $i"))
point(tri::Triangle{T}, a1::T, a2::T, a3::T) where {T<:Real} = vertex(tri, 1) * a1 + vertex(tri, 2) * a2 + vertex(tri, 3) * a3
uv(tri::Triangle{T}, a1::T, a2::T, a3::T) where {T<:Real} = uvvertex(tri, 1) * a1 + uvvertex(tri, 2) * a2 + uvvertex(tri, 3) * a3
normal(tri::Triangle{T}) where {T<:Real} = tri.normal

# returns an Interval representing a half-space. If the ray and normal are pointing in opposite directions then the ray is entering the surface and the half-space will be Interval(surface intersection, +inf). If the ray and normal are pointing in the same directions then the ray is leaving the surface and the half-space will be Interval(-inf,surface intersection)
function surfaceintersection(tri::Triangle{T}, r::AbstractRay{T,3}) where {T<:Real}
    # Möller–Trumbore ray-triangle intersection algorithm
    D = direction(r)
    BA = tri.BA
    CA = tri.CA
    h = cross(D, CA)
    a = dot(BA, h)
    if samepoint(a, zero(T))
        return EmptyInterval(T) # This ray is parallel to this triangle
    end
    f = one(T) / a
    s = origin(r) - tri.A
    δ = f * dot(s, h)
    if (δ < zero(T) || δ > one(T))
        return EmptyInterval(T) # outside tri
    end
    q = cross(s, BA)
    β = f * dot(D, q)
    if (β < zero(T) || δ + β > one(T))
        return EmptyInterval(T) # outside tri
    end
    # At this stage we can compute t to find out where the intersection point is on the line.
    t = f * dot(CA, q)
    if (t > zero(T))
        u, v = δ * uvvertex(tri, 2) + β * uvvertex(tri, 3) + (1 - δ - β) * uvvertex(tri, 1)
        # these are the uvcoords want to use as starting point for newton iteration. Don't want the u,v which has already computed for the basis vecors BA,CA.
        n̂ = normal(tri)
        int = Intersection(t, point(r, t), n̂, u, v, NullInterface(T))
        if dot(n̂, D) < zero(T)
            return positivehalfspace(int)
        else
            return rayorigininterval(int)
        end
    else
        return EmptyInterval(T) # wrong side of ray origin
    end
end

function makemesh(t::Triangle{T}, ::Int = 0) where {T<:Real}
    return TriangleMesh([t])
end

##########################################################################################################

"""
    TriangleMesh{T} <: Surface{T}

An array of [`Triangle`](@ref)s forming a mesh.
Used for visualization purposes only.

```julia
TriangleMesh(tris::Vector{Triangle{T}})
```
"""
struct TriangleMesh{T} <: Surface{T}
    # TODO - make this more efficient, e.g. store unique vertices and indices
    triangles::Vector{Triangle{T}}

    function TriangleMesh(tris::Vector{Triangle{T}}) where {T<:Real}
        return new{T}(tris)
    end
end
export TriangleMesh
Base.eltype(::TriangleMesh{T}) where {T<:Real} = Triangle{T}

function makiemesh(tmesh::TriangleMesh{T}) where {T<:Real}
    # TODO probably should unify shared verts at this point or already have this stored in the trimesh
    points = Vector{SVector{3,T}}(undef, length(tmesh.triangles) * 3)
    indices = Array{UInt32,2}(undef, length(tmesh.triangles), 3)
    @inbounds @simd for i in 0:(length(tmesh.triangles) - 1)
        t = tmesh.triangles[i + 1]
        points[i * 3 + 1] = vertex(t, 1)
        points[i * 3 + 2] = vertex(t, 2)
        points[i * 3 + 3] = vertex(t, 3)
        indices[i + 1, :] = [i * 3 + 1, i * 3 + 2, i * 3 + 3]
    end
    return (points, indices)
end
