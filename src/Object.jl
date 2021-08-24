using Base: uniontypes
using Random
using UUIDs

export Object, Properties

uuid() = uuid1(Random.GLOBAL_RNG)

struct Object{T}
    object::T
    id::UUID

    Object(object::T, id::UUID = uuid()) where {T} = new{T}(object, id)
end

Base.show(io::IO, object::Object) = print(io, "[$(object.id)]")

unionobject(T::Type) = Union{T, Object{<:T}}
unionobject(T::Union) = Union{uniontypes(T)..., [Object{<:S} for S in uniontypes(T)]...}

abstract type HasID end
struct WithID <: HasID end
struct WithoutID <: HasID end

hasid(::Type) = WithoutID()
hasid(::Object) = WithID()

const Properties = Dict{UUID,Dict{String,Any}}
