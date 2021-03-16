abstract type AbstractRayInterval{T<:Real} end

"""
    EmptyInterval{T} <: AbstractRayInterval{T}

An interval with no [`Intersection`](@ref)s which is also not infinite.

```
EmptyInterval(T = Float64)
EmptyInterval{T}()
```
"""
struct EmptyInterval{T} <: AbstractRayInterval{T}
    EmptyInterval(::Type{T} = Float64) where {T<:Real} = new{T}()
    EmptyInterval{T}() where {T<:Real} = new{T}()
end
export EmptyInterval

"""
    Interval{T} <: AbstractRayInterval{T}

Datatype representing an interval between two [`IntervalPoint`](@ref)s on a ray.

The lower element can either be [`RayOrigin`](@ref) or an [`Intersection`](@ref).
The upper element can either be an [`Intersection`](@ref) or [`Infinity`](@ref).

```julia
positivehalfspace(int::Intersection) -> Interval with lower = int, upper = Infinity
rayorigininterval(int::Intersection) -> Interval with lower = RayOrigin, upper = int
Interval(low, high)
```

Has the following accessor methods:
```julia
lower(a::Interval{T}) -> Union{RayOrigin{T},Intersection{T,3}}
upper(a::Interval{T}) -> Union{Intersection{T,3},Infinity{T}}
```
"""
struct Interval{R} <: AbstractRayInterval{R}
    lower::Union{RayOrigin{R},Intersection{R,3}} #making this a Union of two concrete types allows the compiler to statically enough space to hold the larger of the two. The entire Interval struct can then be stack allocated.
    upper::Union{Intersection{R,3},Infinity{R}} #same as above

    function Interval(low::S, high::T) where {R<:Real,S<:Union{RayOrigin{R},Intersection{R,3}},T<:Union{Intersection{R,3},Infinity{R}}}
        @assert low <= high
        return new{R}(low, high)
    end
end
export Interval

function Base.show(io::IO, a::Interval{R}) where {R<:Real}
    println(io, Interval{R})
    println(io, "\tLower \n\t\t$(lower(a))")
    println(io, "\tUpper \n\t\t$(upper(a))")
end

Base.eltype(::Interval{R}) where {R} = R

lower(a::Interval) = a.lower
upper(a::Interval) = a.upper

########################################################################################################################

"""
To prevent allocations we have a manually managed pool of arrays of [`Interval`](@ref)s which are used to store values during execution.
The memory is kept allocated and reused across runs of functions like [`trace`](@ref).

`threadedintervalpool` is a global threadsafe pool which is accessed through the functions:
```julia
newinintervalpool!(::Type{T} = Float64, tid::Int = Threads.threadid()) -> Vector{Interval{T}}
indexednewinintervalpool!(::Type{T} = Float64, tid::Int = Threads.threadid()) -> Tuple{Int,Vector{Interval{T}}}
emptyintervalpool!(::Type{T} = Float64, tid::Int = Threads.threadid())
getfromintervalpool([::Type{T} = Float64], id::Int, tid::Int = Threads.threadid()) -> Vector{Interval{T}}
```
"""
struct IntervalPool{T<:Real}
    allocated::Vector{Vector{Interval{T}}}
    unallocated::Vector{Vector{Interval{T}}}

    IntervalPool{T}() where {T<:Real} = new{T}(Vector{Vector{Interval{T}}}(), Vector{Vector{Interval{T}}}())
end

function allocate!(a::IntervalPool{T}) where {T<:Real}
    if length(a.unallocated) === 0
        push!(a.unallocated, Vector{Interval{T}}()) #this will extend the array with a new Vector of intervals.
    end

    temp = pop!(a.unallocated)
    Base.empty!(temp) #resets count in array but doesn't reclaim space
    push!(a.allocated, temp)
    return temp
end

function empty!(a::IntervalPool{T}) where {T<:Real}
    while length(a.allocated) > 0
        temp = pop!(a.allocated)
        Base.empty!(temp) #not strictly necessary since allocate will reset this array to empty before pushing it back on allocated stack. But this ensures the interval pool is in a consistent emptied state upon completion of this function.
        push!(a.unallocated, temp)
    end
end

const threadedintervalpool = [Dict{DataType,IntervalPool}([Float64 => IntervalPool{Float64}()]) for _ in 1:Threads.nthreads()]

function newinintervalpool!(::Type{T} = Float64, tid::Int = Threads.threadid())::Vector{Interval{T}} where {T<:Real}
    if T ∉ keys(threadedintervalpool[tid])
        # if the type of the interval pool has changed then we need to refill it with the correct type
        threadedintervalpool[tid][T] = IntervalPool{T}()
    end
    return allocate!(threadedintervalpool[tid][T])
end

function indexednewinintervalpool!(::Type{T} = Float64, tid::Int = Threads.threadid())::Tuple{Int,Vector{Interval{T}}} where {T<:Real}
    a = newinintervalpool!(T, tid)
    index = length(threadedintervalpool[tid][T].allocated)
    return index, a
end

function emptyintervalpool!(::Type{T} = Float64, tid::Int = Threads.threadid()) where {T}
    if T ∈ keys(threadedintervalpool[tid])
        OpticSim.empty!(threadedintervalpool[tid][T])
    end
end

getfromintervalpool(id::Int, tid::Int = Threads.threadid())::Vector{Interval{Float64}} = getfromintervalpool(Float64, id, tid)

function getfromintervalpool(::Type{T}, id::Int, tid::Int = Threads.threadid())::Vector{Interval{T}} where {T<:Real}
    return threadedintervalpool[tid][T].allocated[id]
end

########################################################################################################################

macro inplaceinsertionsort(array, f)
    esc(quote
        for i in 2:length($array)
            value = $array[i]
            j = i - 1
            while j > 0 && $f($array[j]) > $f(value)
                $array[j + 1] = $array[j]
                j = j - 1
            end
            $array[j + 1] = value
        end
    end)
end

"""
Datatype representing an ordered series of disjoint intervals on a ray.
An arbitrary array of `Interval`s can be input to the constructor and they will automatically be processed into a valid `DisjointUnion` (or a single [`Interval`](@ref) if appropriate).

```julia
DisjointUnion(intervals::AbstractVector{Interval{R}})
```
"""
struct DisjointUnion{T<:Real}
    intervalarrayindex::Int

    function DisjointUnion(a::Interval{R}, b::Interval{R})::Union{DisjointUnion{R},Interval{R}} where {R<:Real}
        temp = newinintervalpool!(R)
        push!(temp, a)
        push!(temp, b)
        return DisjointUnion(temp)
    end

    function DisjointUnion(intervals::AbstractVector{Interval{R}})::Union{DisjointUnion{R},Interval{R}} where {R<:Real}
        if length(intervals) === 1
            return intervals[1]
        end

        @inplaceinsertionsort(intervals, upper)

        index, result = indexednewinintervalpool!(R)

        while length(intervals) > 0
            current = pop!(intervals)

            while length(intervals) > 0
                temp = intervalintersection(current, last(intervals))
                if !(temp isa EmptyInterval{R})
                    current = intervalunion(current, last(intervals))
                    pop!(intervals)
                else
                    break
                end
            end

            push!(result, current)
        end

        if length(result) === 1
            return result[1]
        else
            @inplaceinsertionsort(result, lower)
            return new{R}(index)
        end
    end
end

intervals(a::DisjointUnion{R}) where {R<:Real} = getfromintervalpool(R, a.intervalarrayindex)
function Base.show(io::IO, a::DisjointUnion{R}) where {R<:Real}
    println(io, DisjointUnion{R})
    for i in intervals(a)
        show(io, i)
    end
end

export DisjointUnion

Base.getindex(a::DisjointUnion, index::Int) = intervals(a)[index]
Base.firstindex(a::DisjointUnion) = firstindex(intervals(a))
Base.lastindex(a::DisjointUnion) = lastindex(intervals(a))
Base.length(a::DisjointUnion) = length(intervals(a))
Base.:(==)(a::DisjointUnion, b::DisjointUnion) = intervals(a) == intervals(b) # this works so long as the elements in the intervals array are structs and not arrays.

########################################################################################################################

"""
    isemptyinterval(a) -> Bool

Returns true if `a` is an [`EmptyInterval`](@ref). In performance critical contexts use `a isa EmptyInterval{T}`.

"""
isemptyinterval(::EmptyInterval) = true
isemptyinterval(::Interval) = false
isemptyinterval(::DisjointUnion) = false

"""
    ispositivehalfspace(a) -> Bool

Returns true if `upper(a)` is [`Infinity`](@ref). In performance critical contexts check directly i.e. `upper(a) isa Infinity{T}`.

"""
ispositivehalfspace(a::Interval{T}) where {T<:Real} = upper(a) isa Infinity{T}

"""
    israyorigininterval(a) -> Bool

Returns true if `lower(a)` is [`RayOrigin`](@ref). In performance critical contexts check directly i.e. `lower(a) isa RayOrigin{T}`.

"""
israyorigininterval(a::Interval{T}) where {T<:Real} = lower(a) isa RayOrigin{T}

isinfiniteinterval(a::Interval{T}) where {T<:Real} = israyorigininterval(a) && ispositivehalfspace(a)

positivehalfspace(a::RayOrigin{T}) where {T<:Real} = Interval(a, Infinity(T))
positivehalfspace(a::Intersection{T,N}) where {T<:Real,N} = Interval(a, Infinity(T))

rayorigininterval(a::Intersection{T,N}) where {T<:Real,N} = Interval(RayOrigin(T), a) # special kind of interval that has ray origin as lower bound. α(a::RayOrigin) will always return 0 so CSG operations should work.
rayorigininterval(a::Infinity{T}) where {T<:Real} = Interval(RayOrigin(T), a)

"""
    halfspaceintersection(a::Interval{T}) -> Intersection{T,3}

Returns the [`Intersection`](@ref) from a half space [`Interval`](@ref), throws an error if not a half space.
"""
function halfspaceintersection(a::Interval{T})::Intersection{T,3} where {T<:Real}
    la = lower(a)
    ua = upper(a)
    if ua isa Infinity{T} && !(la isa RayOrigin{T})
        return la
    elseif la isa RayOrigin{T} && !(ua isa Infinity{T})
        return ua
    else
        return throw(ErrorException("Not a half-space: $a"))
    end
end

"""
    closestintersection(a::Union{EmptyInterval{T},Interval{T},DisjointUnion{T}}, ignorenull::Bool = true) -> Union{Nothing,Intersection{T,3}}

Returns the closest [`Intersection`](@ref) from an [`Interval`](@ref) or [`DisjointUnion`](@ref).
Ignores intersection with null interfaces if `ignorenull` is true. Will return `nothing` if there is no valid intersection.
"""
closestintersection(::EmptyInterval, ::Bool = true) = nothing
function closestintersection(a::Interval{T}, ignorenull::Bool = true)::Union{Nothing,Intersection{T,3}} where {T<:Real}
    la = lower(a)
    ua = upper(a)
    if la isa RayOrigin{T}
        if !(ua isa Infinity{T}) && !(ignorenull && interface(ua) isa NullInterface{T})
            return ua
        else
            return nothing
        end
    elseif !(ignorenull && interface(la) isa NullInterface{T})
        return la
    else
        return nothing
    end
end
function closestintersection(a::DisjointUnion{T}, ignorenull::Bool = true)::Union{Nothing,Intersection{T,3}} where {T<:Real}
    for i in intervals(a)
        c = closestintersection(i, ignorenull)
        if c !== nothing
            return c
        end
    end
    return nothing
end
export closestintersection

function reversenormal(a::Interval{T})::Interval{T} where {T<:Real}
    la = lower(a)
    ua = upper(a)
    if la isa RayOrigin{T}
        if ua isa Infinity{T}
            return a
        else
            return Interval(la, reversenormal(ua))
        end
    else
        if ua isa Infinity{T}
            return Interval(reversenormal(la), ua)
        else
            return Interval(reversenormal(la), reversenormal(ua))
        end
    end
end
function reversenormal(a::DisjointUnion{R})::DisjointUnion{R} where {R<:Real}
    intvls = newinintervalpool!(R)
    for i in intervals(a)
        push!(intvls, reversenormal(i))
    end
    return DisjointUnion(intvls)
end
########################################################################################################################

macro intervalintersectionhigh(low)
    esc(quote
        if ua isa Infinity{R}
            if ub isa Infinity{R}
                return Interval($low, Infinity(R))
            else
                if ub <= $low
                    return EmptyInterval(R)
                end
                return Interval($low, ub)
            end
        else
            if ub isa Infinity{R}
                if ua <= $low
                    return EmptyInterval(R)
                end
                return Interval($low, ua)
            else
                high = min(ua, ub)
                if high <= $low
                    return EmptyInterval(R)
                end
                return Interval($low, high)
            end
        end
    end)
end

function intervalintersection(a::Interval{R}, b::Interval{R})::Union{EmptyInterval{R},Interval{R}} where {R<:Real}
    # This method is just doing the below, but to avoid type ambiguities things have to be much more complicated
    # if upper(a) <= lower(b) || upper(b) <= lower(a)
    #     return EmptyInterval(R)
    # end
    # low = max(lower(a), lower(b))
    # high = min(upper(a), upper(b))
    # return Interval(low, high)

    la = lower(a)
    lb = lower(b)
    ua = upper(a)
    ub = upper(b)
    if la isa RayOrigin{R}
        if lb isa RayOrigin{R}
            @intervalintersectionhigh(RayOrigin(R))
        else
            @intervalintersectionhigh(lb)
        end
    else
        if lb isa RayOrigin{R}
            @intervalintersectionhigh(la)
        else
            low = max(la, lb)
            @intervalintersectionhigh(low)
        end
    end
end

intervalintersection(::EmptyInterval{T}, ::EmptyInterval{T}) where {T<:Real} = EmptyInterval(T)
intervalintersection(::EmptyInterval{T}, ::Interval{T}) where {T<:Real} = EmptyInterval(T)
intervalintersection(::Interval{T}, ::EmptyInterval{T}) where {T<:Real} = EmptyInterval(T)
intervalintersection(::EmptyInterval{T}, ::DisjointUnion{T}) where {T<:Real} = EmptyInterval(T)
intervalintersection(::DisjointUnion{T}, ::EmptyInterval{T}) where {T<:Real} = EmptyInterval(T)
intervalintersection(a::Interval{T}, b::DisjointUnion{T}) where {T<:Real} = intervalintersection(a, intervals(b))
intervalintersection(a::DisjointUnion{T}, b::Interval{T}) where {T<:Real} = intervalintersection(b, intervals(a))
intervalintersection(a::DisjointUnion{T}, b::DisjointUnion{T}) where {T<:Real} = intervalintersection(intervals(a), intervals(b))

function intervalintersection(a::Interval{T}, b::AbstractVector{Interval{T}})::Union{EmptyInterval{T},Interval{T},DisjointUnion{T}} where {T<:Real}
    temp = nothing
    int1 = nothing
    for bint in b
        intsct = intervalintersection(a, bint)
        if !(intsct isa EmptyInterval{T})
            if int1 === nothing
                int1 = intsct
            elseif temp === nothing
                temp = newinintervalpool!(T)
                push!(temp, int1)
                push!(temp, intsct)
            else
                push!(temp, intsct)
            end
        end
    end
    if int1 === nothing
        return EmptyInterval(T)
    elseif int1 !== nothing && (temp === nothing)
        return int1
    else
        return DisjointUnion(temp)
    end
end

function intervalintersection(a::AbstractVector{Interval{T}}, b::AbstractVector{Interval{T}})::Union{EmptyInterval{T},Interval{T},DisjointUnion{T}} where {T<:Real}
    temp = newinintervalpool!(T)
    for aint in a
        for bint in b
            intsct = intervalintersection(aint, bint)
            if !(intsct isa EmptyInterval{T})
                push!(temp, intsct)
            end
        end
    end
    if length(temp) == 0
        return EmptyInterval(T)
    else
        return DisjointUnion(temp)
    end
end

########################################################################################################################

macro intervalunionhigh(low)
    esc(quote
        if ua isa Infinity{R}
            if ub isa Infinity{R}
                return Interval($low, Infinity(R))
            else
                if ub < la
                    return DisjointUnion(a, b)
                end
                return Interval($low, Infinity(R))
            end
        else
            if ub isa Infinity{R}
                if ua < lb
                    return DisjointUnion(a, b)
                end
                return Interval($low, Infinity(R))
            else
                if ua < lb || ub < la
                    return DisjointUnion(a, b)
                end
                return Interval($low, max(ua, ub))
            end
        end
    end)
end

function intervalunion(a::Interval{R}, b::Interval{R})::Union{Interval{R},DisjointUnion{R}} where {R<:Real}
    # This method is just doing the below, but to avoid type ambiguities things have to be much more complicated
    # if upper(a) < lower(b) || upper(b) < lower(a)
    #     return DisjointUnion(a, b)
    # else
    #     low = min(lower(a), lower(b))
    #     high = max(upper(a), upper(b))
    #     return Interval(low, high)
    # end

    ua = upper(a)
    ub = upper(b)
    la = lower(a)
    lb = lower(b)
    if la isa RayOrigin{R}
        if lb isa RayOrigin{R}
            @intervalunionhigh(RayOrigin(R))
        else
            @intervalunionhigh(RayOrigin(R))
        end
    else
        if lb isa RayOrigin{R}
            @intervalunionhigh(RayOrigin(R))
        else
            low = min(la, lb)
            @intervalunionhigh(low)
        end
    end
end

intervalunion(a::Interval{T}, b::DisjointUnion{T}) where {T<:Real} = intervalunion(b, a)
function intervalunion(a::DisjointUnion{T}, b::Interval{T}) where {T<:Real}
    temp = newinintervalpool!(T)
    for intvl in intervals(a)
        push!(temp, intvl)
    end
    push!(temp, b)
    return DisjointUnion(temp)
end
function intervalunion(a::DisjointUnion{T}, b::DisjointUnion{T}) where {T<:Real}
    temp = newinintervalpool!(T)
    for intvl in intervals(a)
        push!(temp, intvl)
    end
    for intvl in intervals(b)
        push!(temp, intvl)
    end
    return DisjointUnion(temp)
end
intervalunion(::EmptyInterval{T}, ::EmptyInterval{T}) where {T<:Real} = EmptyInterval(T)
intervalunion(::EmptyInterval{T}, b::Interval{T}) where {T<:Real} = b
intervalunion(a::Interval{T}, ::EmptyInterval{T}) where {T<:Real} = a
intervalunion(::EmptyInterval{T}, b::DisjointUnion{T}) where {T<:Real} = b
intervalunion(a::DisjointUnion{T}, ::EmptyInterval{T}) where {T<:Real} = a

########################################################################################################################

function intervalcomplement(a::DisjointUnion{T}) where {T<:Real}
    res = rayorigininterval(Infinity(T))
    for intrval in intervals(a)
        @assert !(intrval isa EmptyInterval{T}) "should never have an empty interval in a DisjointUnion"
        compl = intervalcomplement(intrval)::Union{Interval{T},DisjointUnion{T}}
        if compl isa Interval{T}
            res = intervalintersection(res, compl)
        else
            res = intervalintersection(res, compl)
        end
    end
    return res
end

intervalcomplement(::EmptyInterval{T}) where {T<:Real} = Interval(RayOrigin(T), Infinity(T))
function intervalcomplement(a::Interval{T}) where {T<:Real}
    ua = upper(a)
    la = lower(a)
    if la isa RayOrigin{T}
        if ua isa Infinity{T}
            return EmptyInterval(T)
        else
            return positivehalfspace(ua) # this is not strictly correct
        end
    else
        if ua isa Infinity{T}
            return rayorigininterval(la)
        else
            return DisjointUnion(rayorigininterval(la), positivehalfspace(ua))
        end
    end
end

difference(a::DisjointUnion, b::DisjointUnion) = intervalintersection(a, intervalcomplement(b))
