using Test
using FiniteDifferences
using StaticArrays
using LinearAlgebra
using Suppressor
using Random
using Unitful

using OpticSim
# curve/surface imports
using OpticSim: findspan, makemesh, knotstoinsert, coefficients, inside, quadraticroots, tobeziersegments, evalcsg, makiemesh
# interval imports
using OpticSim: α, halfspaceintersection, positivehalfspace, lower, upper, EmptyInterval, rayorigininterval, intervalcomplement, intervalintersection, intervalunion, RayOrigin, Infinity, Intersection
using OpticSim.TestData: intersectionat
# bounding box imports
using OpticSim: doesintersect
# RBT imports
using OpticSim: rotmat, rotmatd
# lens imports
using OpticSim: reflectedray, snell, refractedray, trace, intersection, pathlength, power
#Bounding volume hierarchy imports
# using OpticSim: partition!

######################################################################################################

const COMP_TOLERANCE = 25 * eps(Float64)
const RTOLERANCE = 1e-10
const ATOLERANCE = 10 * eps(Float64)
const SEED = 12312487
const ALL_TESTS = isempty(ARGS) || "all" in ARGS || "All" in ARGS || "ALL" in ARGS

"""Optional testset: if the name of the testset is passed as an arg then it will be executed. If none are specified, or "all" is given as an arg then all optional testsets will be executed """
macro otestset(name, expr)
    quote
        if ALL_TESTS || $name in ARGS
            @testset $name begin
                $expr
            end
        end
    end
end

"""Creates a 3D vector uniformly distributed on the sphere by rejection sampling, i.e., discarding all points with norm > 1.0"""
function randunit()
    let v = rand(3)
        while (norm(v) > 1.0)
            v = rand(3)
        end
        return normalize(v)
    end
end

"""Evalute all functions not requiring arguments in a given module and test they don't throw anything"""
macro test_all_no_arg_functions(m)
    quote
        for n in names($m, all = true)
            # get all the valid functions
            if occursin("#", string(n)) || string(n) == string(nameof($m))
                continue
            end
            f = Core.eval($m, n)
            # get the methods of this function
            for meth in methods(f)
                # if this method has no args then try and evaluate it
                if meth.nargs == 1
                    # suppress STDOUT to prevent loads of stuff spamming the console
                    @test (@suppress_out begin
                        f()
                    end;
                    true)
                end
            end
        end
    end
end

######################################################################################################


# FIXING ERRONEOUS UNBOUND TYPE ERRORS THAT OCCUR WITH VARARG
#=
The set will contain something like this:
    MultiHologramInterface(interfaces::Vararg{HologramInterface{T}, N}) where {T<:Real, N}
To get the signature run:
    methods(OpticSim.MultiHologramInterface).ms[I].sig            where I is typically 1 or 2 depending on the number of methods
This gives:
    Tuple{Type{MultiHologramInterface}, Vararg{HologramInterface{T}, N}} where N where T<:Real
We then have to use which to find the method:
    which(OpticSim.MultiHologramInterface, Tuple{Vararg{HologramInterface{T}, N}} where N where T<:Real)
and pop it from the set
=#

@otestset "JuliaLang" begin
    # ensure there aren't any ambiguities or unbound args
    let unbound = Set{Method}(detect_unbound_args(OpticSim))
        # Pretty hacky way to ignore specific methods, for some weird reason the unbound args check seems to fail for some (seemingly random) methods with Vararg arguments
        # here we ignore the specific methods which we know are ok but are still failing
        if VERSION >= v"1.6.0-DEV"
            pop!(unbound, which(OpticSim.LensAssembly, Tuple{Vararg{Union{CSGTree{T},LensAssembly{T},Surface{T}},N} where N} where {T<:Real}))
            pop!(unbound, which(OpticSim.RayListSource, Tuple{Vararg{OpticalRay{T,3},N} where N} where {T<:Real}))
            pop!(unbound, which(OpticSim.OpticalSourceGroup, Tuple{Vararg{OpticalRayGenerator{T},N} where N} where {T<:Real}))
            pop!(unbound, which(OpticSim.PixelSource, Tuple{Vararg{P,N} where N} where {P<:OpticalRayGenerator{T}} where {T<:Real}))
            pop!(unbound, which(OpticSim.MultiHologramInterface, Tuple{Vararg{HologramInterface{T},N}} where {N} where {T<:Real}))
        end
        # and ignore any generate methods created due to default or keyword arguments
        for m in unbound
            if occursin("#", string(m.name))
                pop!(unbound, m)
            end
        end
        @test unbound == Set()
    end
    let ambiguous = Set{Method}(detect_ambiguities(OpticSim))
        @test ambiguous == Set()
    end
    let unbound = Set{Method}(detect_unbound_args(OpticSim.Zernike))
        @test unbound == Set()
    end
    let ambiguous = Set{Method}(detect_ambiguities(OpticSim.Zernike))
        @test ambiguous == Set()
    end
    let unbound = Set{Method}(detect_unbound_args(OpticSim.QType))
        @test unbound == Set()
    end
    let ambiguous = Set{Method}(detect_ambiguities(OpticSim.QType))
        @test ambiguous == Set()
    end
    let unbound = Set{Method}(detect_unbound_args(OpticSim.Vis))
        # Pretty hacky way to ignore specific methods, for some weird reason the unbound args check seems to fail for some (seemingly random) methods with Vararg arguments
        # here we ignore the specific methods which we know are ok but are still failing
        pop!(unbound, which(OpticSim.Vis.drawcurves, Tuple{Vararg{Spline{P,S,N,M},N1} where N1} where {M} where {N} where {S} where {P}))
        pop!(unbound, which(OpticSim.Vis.draw, Tuple{Vararg{S,N} where N} where {S<:Union{OpticSim.Surface{T},OpticSim.TriangleMesh{T}}} where {T<:Real}))
        pop!(unbound, which(OpticSim.Vis.draw!, Tuple{OpticSim.Vis.MakieLayout.LScene,Vararg{S,N} where N} where {S<:Union{OpticSim.Surface{T},OpticSim.TriangleMesh{T}}} where {T<:Real}))
        # and ignore any generate methods created due to default or keyword arguments
        for m in unbound
            if occursin("#", string(m.name))
                pop!(unbound, m)
            end
        end
        @test unbound == Set()
    end
    let ambiguous = Set{Method}(detect_ambiguities(OpticSim.Vis))
        @test ambiguous == Set()
    end
    let unbound = Set{Method}(detect_unbound_args(OpticSim.TestData))
        @test unbound == Set()
    end
    let ambiguous = Set{Method}(detect_ambiguities(OpticSim.TestData))
        @test ambiguous == Set()
    end
end # testset JuliaLang

# @otestset "BVH" begin
#     @testset "partition!" begin
#         split = 0.5

#         for i in 1:100000
#             a = rand(5)
#             b = copy(a)
#             badresult::Bool = false

#             lower, upper = partition!(a, (x) -> x, split)

#             if lower !== nothing
#                 for i in lower
#                     if i >= split
#                         badresult = true
#                         break
#                     end
#                 end
#             end

#             if upper !== nothing
#                 for i in upper
#                     if i <= split
#                         badresult = true
#                     end
#                 end
#             end

#             if badresult
#                 throw(ErrorException("array didn't partition: $(b)"))
#             end

#         end
#     end
# end # testset BVH

@otestset "TestData" begin
    @test_all_no_arg_functions TestData
end # testset TestData

# TODO may want to test this but it's very slow on azure
# @otestset "Examples" begin
#     @test_all_no_arg_functions Examples
# end # testset Examples

@otestset "General" begin
    @testset "QuadraticRoots" begin
        Random.seed!(SEED)
        similarroots(r1, r2, x1, x2) = isapprox([r1, r2], [x1, x2], rtol = 1e-9) || isapprox([r2, r1], [x1, x2], rtol = 1e-9)

        for i in 1:10000
            r1, r2, scale = rand(3) .- 0.5
            a, b, c = scale .* (1, r1 + r2, r1 * r2)
            x1, x2 = quadraticroots(a, b, c)
            @test similarroots(-r1, -r2, x1, x2)
        end

        a, b, c = 1, 0, -1
        x1, x2 = quadraticroots(a, b, c)
        @test similarroots(1, -1, x1, x2)

        a, b, c = 1, 2, 1 # (x+1)(x+1)= x^2 + 2x + 1, double root at one
        x1, x2 = quadraticroots(a, b, c)
        @test similarroots(-1, -1, x1, x2)
    end # testset QuadraticRoots

    @testset "RigidBodyTransform" begin
        Random.seed!(SEED)
        @test isapprox(rotmatd(180, 0, 0), [1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 -1.0], rtol = RTOLERANCE, atol = ATOLERANCE)
        @test isapprox(rotmatd(0.0, 180.0, 0.0), [-1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 -1.0], rtol = RTOLERANCE, atol = ATOLERANCE)
        @test isapprox(rotmatd(0, 0, 180), [-1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 1.0], rtol = RTOLERANCE, atol = ATOLERANCE)
        @test isapprox(rotmatd(0, 90, 0), [0.0 0.0 1.0; 0.0 1.0 0.0; -1.0 0.0 0.0], rtol = RTOLERANCE, atol = ATOLERANCE)
        @test isapprox(rotmatd(45, -45, 45), [0.5 -0.8535533905932737 0.1464466094067261; 0.5 0.14644660940672644 -0.8535533905932737; 0.7071067811865475 0.5 0.5], rtol = RTOLERANCE, atol = ATOLERANCE)
        x, y, z = rand(3)
        @test isapprox(rotmatd(x * 180 / π, y * 180 / π, z * 180 / π), rotmat(x, y, z), rtol = RTOLERANCE, atol = ATOLERANCE)

        @test isapprox(RigidBodyTransform(rotmatd(0, 90, 0), SVector(0.0, 0.0, 1.0)) * SVector(1.0, 0.0, 0.0), [0.0, 0.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE)

        ta = RigidBodyTransform(rotmatd(0, 90, 0), SVector(0.0, 0.0, 1.0))
        tb = RigidBodyTransform(rotmatd(90, 0, 0), SVector(1.0, 0.0, 0.0))
        @test isapprox(collect(ta * tb), collect(RigidBodyTransform(rotmatd(90, 90, 0), SVector(0.0, 0.0, 0.0))), rtol = RTOLERANCE, atol = ATOLERANCE)

        @test isapprox(collect(ta * inv(ta)), collect(identitytransform()), rtol = RTOLERANCE, atol = ATOLERANCE)
        @test isapprox(collect(tb * inv(tb)), collect(identitytransform()), rtol = RTOLERANCE, atol = ATOLERANCE)
        @test isapprox(collect(inv(ta)), collect(RigidBodyTransform(rotmatd(0, -90, 0), SVector(1.0, 0.0, 0.0))), rtol = RTOLERANCE, atol = ATOLERANCE)
    end # testset RigidBodyTransform

    @testset "Interval" begin
        pt1 = Intersection(0.3, [4.0, 5.0, 6.0], normalize(rand(3)), 0.4, 0.5, NullInterface())
        pt2 = Intersection(0.5, [1.0, 2.0, 3.0], normalize(rand(3)), 0.2, 0.3, NullInterface())
        @test pt1 < pt2 && pt1 <= pt2 && pt2 > pt1 && pt2 >= pt1 && pt1 != pt2 && pt1 <= pt1
        intvl2 = positivehalfspace(pt1)
        @test α(halfspaceintersection(intvl2)) == α(pt1)

        ### interval intersection
        ## interval/interval
        # one in another
        a = Interval(intersectionat(0.1), intersectionat(0.2))
        b = Interval(intersectionat(0.0), intersectionat(0.3))
        @test intervalintersection(a, b) == intervalintersection(b, a) == a
        # overlap
        a = Interval(intersectionat(0.1), intersectionat(0.3))
        b = Interval(intersectionat(0.2), intersectionat(0.4))
        @test intervalintersection(a, b) == intervalintersection(b, a) == Interval(intersectionat(0.2), intersectionat(0.3))
        # no overlap
        a = Interval(intersectionat(0.1), intersectionat(0.2))
        b = Interval(intersectionat(0.3), intersectionat(0.4))
        @test intervalintersection(a, b) == intervalintersection(b, a) isa EmptyInterval
        # start == end
        a = Interval(intersectionat(0.1), intersectionat(0.2))
        b = Interval(intersectionat(0.2), intersectionat(0.3))
        @test intervalintersection(a, b) == intervalintersection(b, a) isa EmptyInterval
        ## interval/du
        # encompass one
        a = Interval(intersectionat(0.0), intersectionat(0.3))
        b = DisjointUnion(Interval(intersectionat(0.1), intersectionat(0.2)), Interval(intersectionat(0.5), intersectionat(0.6)))
        @test intervalintersection(a, b) == intervalintersection(b, a) == b[1]
        # encompass two
        a = Interval(intersectionat(0.0), intersectionat(0.8))
        b = DisjointUnion(Interval(intersectionat(0.1), intersectionat(0.2)), Interval(intersectionat(0.5), intersectionat(0.6)))
        @test intervalintersection(a, b) == intervalintersection(b, a) == b
        # overlap one
        a = Interval(intersectionat(0.0), intersectionat(0.2))
        b = DisjointUnion(Interval(intersectionat(0.1), intersectionat(0.3)), Interval(intersectionat(0.5), intersectionat(0.6)))
        @test intervalintersection(a, b) == intervalintersection(b, a) == Interval(intersectionat(0.1), intersectionat(0.2))
        # overlap two
        a = Interval(intersectionat(0.2), intersectionat(0.5))
        b = DisjointUnion(Interval(intersectionat(0.1), intersectionat(0.3)), Interval(intersectionat(0.4), intersectionat(0.6)))
        res = intervalintersection(a, b)
        @test res == intervalintersection(b, a) && res[1] == Interval(intersectionat(0.2), intersectionat(0.3)) && res[2] == Interval(intersectionat(0.4), intersectionat(0.5))
        # no overlap
        a = Interval(intersectionat(0.1), intersectionat(0.2))
        b = DisjointUnion(Interval(intersectionat(0.3), intersectionat(0.4)), Interval(intersectionat(0.5), intersectionat(0.6)))
        @test intervalintersection(a, b) == intervalintersection(b, a) isa EmptyInterval
        ## du/du
        # no overlap
        a = DisjointUnion(Interval(intersectionat(0.0), intersectionat(0.1)), Interval(intersectionat(0.4), intersectionat(0.5)))
        b = DisjointUnion(Interval(intersectionat(0.2), intersectionat(0.3)), Interval(intersectionat(0.6), intersectionat(0.7)))
        @test intervalintersection(a, b) == intervalintersection(b, a) isa EmptyInterval
        # one in a overlaps one in b
        a = DisjointUnion(Interval(intersectionat(0.0), intersectionat(0.2)), Interval(intersectionat(0.4), intersectionat(0.5)))
        b = DisjointUnion(Interval(intersectionat(0.1), intersectionat(0.3)), Interval(intersectionat(0.6), intersectionat(0.7)))
        @test intervalintersection(a, b) == intervalintersection(b, a) == Interval(intersectionat(0.1), intersectionat(0.2))
        # one in a encompasses one in b
        a = DisjointUnion(Interval(intersectionat(0.0), intersectionat(0.3)), Interval(intersectionat(0.4), intersectionat(0.5)))
        b = DisjointUnion(Interval(intersectionat(0.1), intersectionat(0.2)), Interval(intersectionat(0.6), intersectionat(0.7)))
        @test intervalintersection(a, b) == intervalintersection(b, a) == b[1]
        # one in a overlaps two in b
        a = DisjointUnion(Interval(intersectionat(0.2), intersectionat(0.5)), Interval(intersectionat(0.7), intersectionat(0.8)))
        b = DisjointUnion(Interval(intersectionat(0.1), intersectionat(0.3)), Interval(intersectionat(0.4), intersectionat(0.6)))
        res = intervalintersection(a, b)
        @test res == intervalintersection(b, a) && res[1] == Interval(intersectionat(0.2), intersectionat(0.3)) && res[2] == Interval(intersectionat(0.4), intersectionat(0.5))
        # one in a encompasses two in b
        a = DisjointUnion(Interval(intersectionat(0.0), intersectionat(0.5)), Interval(intersectionat(0.7), intersectionat(0.8)))
        b = DisjointUnion(Interval(intersectionat(0.1), intersectionat(0.2)), Interval(intersectionat(0.3), intersectionat(0.4)))
        @test intervalintersection(a, b) == intervalintersection(b, a) == b
        # each in a overlap one in b
        a = DisjointUnion(Interval(intersectionat(0.0), intersectionat(0.2)), Interval(intersectionat(0.4), intersectionat(0.6)))
        b = DisjointUnion(Interval(intersectionat(0.1), intersectionat(0.3)), Interval(intersectionat(0.5), intersectionat(0.7)))
        res = intervalintersection(a, b)
        @test res == intervalintersection(b, a) && res[1] == Interval(intersectionat(0.1), intersectionat(0.2)) && res[2] == Interval(intersectionat(0.5), intersectionat(0.6))
        # each in a encompass one in b
        a = DisjointUnion(Interval(intersectionat(0.0), intersectionat(0.3)), Interval(intersectionat(0.4), intersectionat(0.7)))
        b = DisjointUnion(Interval(intersectionat(0.1), intersectionat(0.2)), Interval(intersectionat(0.5), intersectionat(0.6)))
        @test intervalintersection(a, b) == intervalintersection(b, a) == b
        ## empties
        intvl = positivehalfspace(pt1)
        du = intervalcomplement(Interval(pt1, pt2))
        @test intervalintersection(EmptyInterval(), EmptyInterval()) isa EmptyInterval
        @test intervalintersection(EmptyInterval(), intvl) isa EmptyInterval
        @test intervalintersection(intvl, EmptyInterval()) isa EmptyInterval
        @test intervalintersection(EmptyInterval(), du) isa EmptyInterval
        @test intervalintersection(du, EmptyInterval()) isa EmptyInterval

        ### interval union
        ## interval/interval
        # one in another
        a = Interval(intersectionat(0.1), intersectionat(0.2))
        b = Interval(intersectionat(0.0), intersectionat(0.3))
        @test intervalunion(a, b) == intervalunion(b, a) == b
        # overlap
        a = Interval(intersectionat(0.1), intersectionat(0.3))
        b = Interval(intersectionat(0.2), intersectionat(0.4))
        @test intervalunion(a, b) == intervalunion(b, a) == Interval(intersectionat(0.1), intersectionat(0.4))
        # no overlap
        a = Interval(intersectionat(0.1), intersectionat(0.2))
        b = Interval(intersectionat(0.3), intersectionat(0.4))
        res = intervalunion(a, b)
        @test res == intervalunion(b, a) && res[1] == a && res[2] == b
        # start == end
        a = Interval(intersectionat(0.1), intersectionat(0.2))
        b = Interval(intersectionat(0.2), intersectionat(0.3))
        @test intervalunion(a, b) == intervalunion(b, a) == Interval(intersectionat(0.1), intersectionat(0.3))
        ## interval/du
        # encompass one
        a = Interval(intersectionat(0.0), intersectionat(0.3))
        b = DisjointUnion(Interval(intersectionat(0.1), intersectionat(0.2)), Interval(intersectionat(0.5), intersectionat(0.6)))
        res = intervalunion(a, b)
        @test res == intervalunion(b, a) && res[1] == a && res[2] == b[2]
        # encompass two
        a = Interval(intersectionat(0.0), intersectionat(0.8))
        b = DisjointUnion(Interval(intersectionat(0.1), intersectionat(0.2)), Interval(intersectionat(0.5), intersectionat(0.6)))
        @test intervalunion(a, b) == intervalunion(b, a) == a
        # overlap one
        a = Interval(intersectionat(0.0), intersectionat(0.2))
        b = DisjointUnion(Interval(intersectionat(0.1), intersectionat(0.3)), Interval(intersectionat(0.5), intersectionat(0.6)))
        res = intervalunion(a, b)
        @test res == intervalunion(b, a) && res[1] == Interval(intersectionat(0.0), intersectionat(0.3)) && res[2] == b[2]
        # overlap two
        a = Interval(intersectionat(0.2), intersectionat(0.5))
        b = DisjointUnion(Interval(intersectionat(0.1), intersectionat(0.3)), Interval(intersectionat(0.4), intersectionat(0.6)))
        @test intervalunion(a, b) == intervalunion(b, a) == Interval(intersectionat(0.1), intersectionat(0.6))
        # no overlap
        a = Interval(intersectionat(0.1), intersectionat(0.2))
        b = DisjointUnion(Interval(intersectionat(0.3), intersectionat(0.4)), Interval(intersectionat(0.5), intersectionat(0.6)))
        res = intervalunion(a, b)
        @test res == intervalunion(b, a) && res[1] == a && res[2] == b[1] && res[3] == b[2]
        ## du/du
        # no overlap
        a = DisjointUnion(Interval(intersectionat(0.0), intersectionat(0.1)), Interval(intersectionat(0.4), intersectionat(0.5)))
        b = DisjointUnion(Interval(intersectionat(0.2), intersectionat(0.3)), Interval(intersectionat(0.6), intersectionat(0.7)))
        res = intervalunion(a, b)
        @test res == intervalunion(b, a) && res[1] == a[1] && res[2] == b[1] && res[3] == a[2] && res[4] == b[2]
        # one in a overlaps one in b
        a = DisjointUnion(Interval(intersectionat(0.0), intersectionat(0.2)), Interval(intersectionat(0.4), intersectionat(0.5)))
        b = DisjointUnion(Interval(intersectionat(0.1), intersectionat(0.3)), Interval(intersectionat(0.6), intersectionat(0.7)))
        res = intervalunion(a, b)
        @test res == intervalunion(b, a) && res[1] == Interval(intersectionat(0.0), intersectionat(0.3)) && res[2] == a[2] && res[3] == b[2]
        # one in a encompasses one in b
        a = DisjointUnion(Interval(intersectionat(0.0), intersectionat(0.3)), Interval(intersectionat(0.4), intersectionat(0.5)))
        b = DisjointUnion(Interval(intersectionat(0.1), intersectionat(0.2)), Interval(intersectionat(0.6), intersectionat(0.7)))
        res = intervalunion(a, b)
        @test res == intervalunion(b, a) && res[1] == a[1] && res[2] == a[2] && res[3] == b[2]
        # one in a overlaps two in b
        a = DisjointUnion(Interval(intersectionat(0.2), intersectionat(0.5)), Interval(intersectionat(0.7), intersectionat(0.8)))
        b = DisjointUnion(Interval(intersectionat(0.1), intersectionat(0.3)), Interval(intersectionat(0.4), intersectionat(0.6)))
        res = intervalunion(a, b)
        @test res == intervalunion(b, a) && res[1] == Interval(intersectionat(0.1), intersectionat(0.6)) && res[2] == a[2]
        # one in a encompasses two in b
        a = DisjointUnion(Interval(intersectionat(0.0), intersectionat(0.5)), Interval(intersectionat(0.7), intersectionat(0.8)))
        b = DisjointUnion(Interval(intersectionat(0.1), intersectionat(0.2)), Interval(intersectionat(0.3), intersectionat(0.4)))
        @test intervalunion(a, b) == intervalunion(b, a) == a
        # each in a overlap one in b
        a = DisjointUnion(Interval(intersectionat(0.0), intersectionat(0.2)), Interval(intersectionat(0.4), intersectionat(0.6)))
        b = DisjointUnion(Interval(intersectionat(0.1), intersectionat(0.3)), Interval(intersectionat(0.5), intersectionat(0.7)))
        res = intervalunion(a, b)
        @test res == intervalunion(b, a) && res[1] == Interval(intersectionat(0.0), intersectionat(0.3)) && res[2] == Interval(intersectionat(0.4), intersectionat(0.7))
        # each in a encompass one in b
        a = DisjointUnion(Interval(intersectionat(0.0), intersectionat(0.3)), Interval(intersectionat(0.4), intersectionat(0.7)))
        b = DisjointUnion(Interval(intersectionat(0.1), intersectionat(0.2)), Interval(intersectionat(0.5), intersectionat(0.6)))
        @test intervalunion(a, b) == intervalunion(b, a) == a
        ## empties
        @test intervalunion(EmptyInterval(), EmptyInterval()) isa EmptyInterval
        @test intervalunion(EmptyInterval(), intvl) == intvl
        @test intervalunion(intvl, EmptyInterval()) == intvl
        @test intervalunion(EmptyInterval(), du) == du
        @test intervalunion(du, EmptyInterval()) == du

        # interval complement
        int = intervalcomplement(EmptyInterval())
        @test lower(int) isa RayOrigin && upper(int) isa Infinity

        int = intervalcomplement(rayorigininterval(Infinity()))
        @test int isa EmptyInterval

        int = intervalcomplement(rayorigininterval(pt1))
        @test lower(int) == pt1 && upper(int) isa Infinity

        int = intervalcomplement(positivehalfspace(pt1))
        @test lower(int) isa RayOrigin && upper(int) == pt1

        int = intervalcomplement(Interval(pt1, pt2))
        @test lower(int[1]) isa RayOrigin && upper(int[1]) == pt1 && lower(int[2]) == pt2 && upper(int[2]) isa Infinity

        a = DisjointUnion(Interval(intersectionat(0.1), intersectionat(0.3)), Interval(intersectionat(0.4), intersectionat(0.7)))
        res = intervalcomplement(a)
        @test res[1] == rayorigininterval(intersectionat(0.1)) && res[2] == Interval(intersectionat(0.3), intersectionat(0.4)) && res[3] == positivehalfspace(intersectionat(0.7))

        a = DisjointUnion(rayorigininterval(intersectionat(0.2)), Interval(intersectionat(0.4), intersectionat(0.7)))
        res = intervalcomplement(a)
        @test res[1] == Interval(intersectionat(0.2), intersectionat(0.4)) && res[2] == positivehalfspace(intersectionat(0.7))
    end # testset interval
end # testset General

@otestset "SurfaceDefs" begin

    @testset "Knots" begin
        # find span
        knots = KnotVector{Int64}([0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5])
        function apply(knots::KnotVector, curveorder, u, correctindex)
            knotnum = findspan(knots, curveorder, u)
            return knotnum == correctindex
        end
        @test apply(knots, 2, 2.5, 5) && apply(knots, 2, 0, 3) && apply(knots, 2, 4, 6) && apply(knots, 2, 5, 8)

        # insert knots
        knots = [0, 0, 0, 0, 1, 2, 3, 3, 4, 4, 5, 5, 6, 7, 7, 7, 7]
        insertions = knotstoinsert(knots, 3)
        @test insertions == [(5, 2), (6, 2), (7, 1), (9, 1), (11, 1), (13, 2)]
    end

    @testset "BSpline" begin
        # bspline conversion
        correctverts = [[0.0, 0.0, 0.0], [0.0, 0.49625, 0.0], [0.6625, 0.49625, 0.65625], [0.0, 0.0, 0.0], [0.6625, 0.49625, 0.65625], [0.6625, 0.0, 0.0], [0.6625, 0.0, 0.0], [0.6625, 0.49625, 0.65625], [1.33, 0.49625, 0.0], [0.6625, 0.0, 0.0], [1.33, 0.49625, 0.0], [1.33, 0.0, 0.0], [0.0, 0.49625, 0.0], [0.0, 1.0, 0.0], [0.6625, 1.0, 0.0], [0.0, 0.49625, 0.0], [0.6625, 1.0, 0.0], [0.6625, 0.49625, 0.65625], [0.6625, 0.49625, 0.65625], [0.6625, 1.0, 0.0], [1.33, 1.0, 0.0], [0.6625, 0.49625, 0.65625], [1.33, 1.0, 0.0], [1.33, 0.49625, 0.0]]
        correcttris = [
            0x00000001 0x00000002 0x00000003
            0x00000004 0x00000005 0x00000006
            0x00000007 0x00000008 0x00000009
            0x0000000a 0x0000000b 0x0000000c
            0x0000000d 0x0000000e 0x0000000f
            0x00000010 0x00000011 0x00000012
            0x00000013 0x00000014 0x00000015
            0x00000016 0x00000017 0x00000018
        ]
        surf = TestData.bsplinesurface()
        points, triangles = makiemesh(makemesh(surf, 2))
        segments = tobeziersegments(surf) # TODO just testing that this evaluates for now
        @test triangles == correcttris
        @test all(isapprox.(points, correctverts, rtol = RTOLERANCE, atol = ATOLERANCE))
    end # testset BSpline

    @testset "PowerBasis" begin
        # polynomial is 3x^2 + 2x + 1
        coeff = [1.0 2.0 3.0]
        curve = PowerBasisCurve{OpticSim.Euclidean,Float64,1,2}(coeff)
        @test !all(isapprox.(coeff, coefficients(curve, 1), rtol = RTOLERANCE, atol = ATOLERANCE)) # TODO not sure if this is correct/what this is testing?
        correct = true
        for x in -10.0:0.1:10
            exact = 3 * x^2 + 2 * x + 1
            calculated = point(curve, x)[1]
            @test isapprox(exact, calculated, rtol = RTOLERANCE, atol = ATOLERANCE)
        end
    end # testset PowerBasis

    @testset "BezierSurface" begin
        surf = TestData.beziersurface()
        fdm = central_fdm(10, 1)
        for u in 0:0.1:1, v in 0:0.1:1
            (du, dv) = partials(surf, u, v)
            fu(u) = point(surf, u, v)
            fv(v) = point(surf, u, v)
            # compute accurate finite difference approximation to derivative
            fdu, fdv = (fdm(fu, u), fdm(fv, v))
            @test isapprox(du, fdu, rtol = 1e-12, atol = 2 * eps(Float64))
            @test isapprox(dv, fdv, rtol = 1e-12, atol = 2 * eps(Float64))
        end
    end # testset BezierSurface

    @testset "ZernikeSurface" begin
        surf = TestData.zernikesurface1a() # with normradius
        fdm = central_fdm(10, 1)
        for ρ in 0.05:0.1:0.95, ϕ in 0:(π / 10):(2π)
            (dρ, dϕ) = partials(surf, ρ, ϕ)
            fρ(ρ) = point(surf, ρ, ϕ)
            fϕ(ϕ) = point(surf, ρ, ϕ)
            # compute accurate finite difference approximation to derivative
            fdρ, fdϕ = (fdm(fρ, ρ), fdm(fϕ, ϕ))
            @test isapprox(dρ, fdρ, rtol = 1e-12, atol = 2 * eps(Float64))
            @test isapprox(dϕ, fdϕ, rtol = 1e-12, atol = 2 * eps(Float64))
        end

        @test !any(isnan.(normal(TestData.conicsurface(), 0.0, 0.0)))
    end # testset ZernikeSurface

    @testset "QTypeSurface" begin
        # test predefined values against papers
        # Forbes, G. W. "Robust, efficient computational methods for axially symmetric optical aspheres." OpticSim express 18.19 (2010): 19700-19712.
        # Forbes, G. W. "Characterizing the shape of freeform optics." OpticSim express 20.3 (2012): 2483-2499.

        f0_true = [2, sqrt(19 / 4), 4 * sqrt(10 / 19), 1 / 2 * sqrt(509 / 10), 6 * sqrt(259 / 509), 1 / 2 * sqrt(25607 / 259)]
        for i in 1:length(f0_true)
            @test isapprox(OpticSim.QType.f0(i - 1), f0_true[i], rtol = 2 * eps(Float64))
        end
        g0_true = [-1 / 2, -5 / (2 * sqrt(19)), -17 / (2 * sqrt(190)), -91 / (2 * sqrt(5090)), -473 / (2 * sqrt(131831))]
        for i in 1:length(g0_true)
            @test isapprox(OpticSim.QType.g0(i - 1), g0_true[i], rtol = 2 * eps(Float64))
        end
        h0_true = [-1 / 2, -6 / sqrt(19), -3 / 2 * sqrt(19 / 10), -20 * sqrt(10 / 509)]
        for i in 1:length(h0_true)
            @test isapprox(OpticSim.QType.h0(i - 1), h0_true[i], rtol = 2 * eps(Float64))
        end

        F_true = [1/4 1/2 27/32 5/4; 15/32 7/8 35/16 511/128; 17/72 35/36 35/16 23/6; 29/40 67/40 243/80 12287/2560]
        N, M = size(F_true)
        for m in 1:M
            for n in 0:(N - 1)
                @test isapprox(OpticSim.QType.F(m, n), F_true[n + 1, m], rtol = 2 * eps(Float64))
            end
        end
        G_true = [1/4 3/8 15/32 35/64; -1/24 -5/48 -7/32 -21/64; -7/40 -7/16 -117/160 -33/32]
        N, M = size(G_true)
        for m in 1:M
            for n in 0:(N - 1)
                @test isapprox(OpticSim.QType.G(m, n), G_true[n + 1, m], rtol = 2 * eps(Float64))
            end
        end
        f_true = [1/2 1/sqrt(2) 3 / 4*sqrt(3 / 2) sqrt(5)/2; 1 / 4*sqrt(7 / 2) 1 / 4*sqrt(19 / 2) 1 / 4*sqrt(185 / 6) 3 / 32*sqrt(427); 1 / 6*sqrt(115 / 14) 1 / 2*sqrt(145 / 38) 1 / 4*sqrt(12803 / 370) 1 / 2*sqrt(2785 / 183); 1 / 5*sqrt(3397 / 230) 1 / 4*sqrt(6841 / 290) 9 / 4*sqrt(14113 / 25606) 1 / 16*sqrt(1289057 / 1114)]
        N, M = size(f_true)
        for m in 1:M
            for n in 0:(N - 1)
                @test isapprox(OpticSim.QType.f(m, n), f_true[n + 1, m], rtol = 2 * eps(Float64))
            end
        end
        g_true = [1/2 3/(4 * sqrt(2)) 5/(4 * sqrt(6)) 7 * sqrt(5)/32; -1/(3 * sqrt(14)) -5/(6 * sqrt(38)) -7 / 4*sqrt(3 / 370) -1 / 2*sqrt(7 / 61); -21 / 10*sqrt(7 / 230) -7 / 4*sqrt(19 / 290) -117 / 4*sqrt(37 / 128030) -33 / 16*sqrt(183 / 2785)]
        N, M = size(g_true)
        for m in 1:M
            for n in 0:(N - 1)
                @test isapprox(OpticSim.QType.g(m, n), g_true[n + 1, m], rtol = 2 * eps(Float64))
            end
        end
        A_true = [2 3 5 7 9 11; -4/3 2/3 2 76/27 85/24 106/25; 9/5 26/15 2 117/50 203/75 108/35; 55/28 66/35 2 536/245 135/56 130/49; 161/81 122/63 2 515/243 143/63 161/66]
        N, M = size(A_true)
        for m in 1:M
            for n in 0:(N - 1)
                @test isapprox(OpticSim.QType.A(m, n), A_true[n + 1, m], rtol = 2 * eps(Float64))
            end
        end
        B_true = [-1 -2 -4 -6 -8 -10; -8/3 -4 -4 -40/9 -5 -28/5; -24/5 -4 -4 -21/5 -112/25 -24/5; -30/7 -4 -4 -144/35 -30/7 -220/49; -112/27 -4 -4 -110/27 -88/21 -13/3]
        N, M = size(B_true)
        for m in 1:M
            for n in 0:(N - 1)
                @test isapprox(OpticSim.QType.B(m, n), B_true[n + 1, m], rtol = 2 * eps(Float64))
            end
        end
        C_true = [NaN NaN NaN NaN NaN NaN; -11/3 -3 -5/3 -35/27 -9/8 -77/75; 0 5/9 7/15 21/50 88/225 13/35; 27/28 21/25 27/35 891/1225 39/56 33/49; 80/81 45/49 55/63 1430/1701 40/49 1105/1386]
        N, M = size(C_true)
        for m in 1:M
            for n in 0:(N - 1)
                @test isapprox(OpticSim.QType.C(m, n), C_true[n + 1, m], rtol = 2 * eps(Float64)) || (isnan(OpticSim.QType.C(m, n)) && isnan(C_true[n + 1, m]))
            end
        end

        # test deriv with finite differences
        surf = TestData.qtypesurface1()
        fdm = central_fdm(10, 1)
        for ρ in 0.05:0.1:0.95, ϕ in 0:(π / 10):(2π)
            (dρ, dϕ) = partials(surf, ρ, ϕ)
            fρ(ρ) = point(surf, ρ, ϕ)
            fϕ(ϕ) = point(surf, ρ, ϕ)
            # compute accurate finite difference approximation to derivative
            fdρ, fdϕ = (fdm(fρ, ρ), fdm(fϕ, ϕ))
            @test isapprox(dρ, fdρ, rtol = 1e-12, atol = 2 * eps(Float64))
            @test isapprox(dϕ, fdϕ, rtol = 1e-12, atol = 2 * eps(Float64))
        end
    end # testset QTypeSurface

    @testset "GridSag" begin
        surf = TestData.gridsagsurfacelinear()
        fdm = central_fdm(10, 1)
        # doesn't enforce C1 across patch boundaries meaning that finite differences won't match at all
        # just test within a patch to make sure partials() is working
        for ρ in 0.02:0.01:0.18, ϕ in 0:(π / 30):(2π)
            (dρ, dϕ) = partials(surf, ρ, ϕ)
            fρ(ρ) = point(surf, ρ, ϕ)
            fϕ(ϕ) = point(surf, ρ, ϕ)
            # compute accurate finite difference approximation to derivative
            fdρ, fdϕ = (fdm(fρ, ρ), fdm(fϕ, ϕ))
            @test isapprox(dρ, fdρ, rtol = 1e-12, atol = 2 * eps(Float64))
            @test isapprox(dϕ, fdϕ, rtol = 1e-12, atol = 2 * eps(Float64))
        end

        surf = TestData.gridsagsurfacebicubic()
        fdm = central_fdm(10, 1)
        # doesn't enforce C2 across patch boundaries meaning that finite differences won't match exactly
        # just test within a patch to make sure partials() is working
        for ρ in 0.02:0.01:0.18, ϕ in 0:(π / 30):(2π)
            (dρ, dϕ) = partials(surf, ρ, ϕ)
            fρ(ρ) = point(surf, ρ, ϕ)
            fϕ(ϕ) = point(surf, ρ, ϕ)
            # compute accurate finite difference approximation to derivative
            fdρ, fdϕ = (fdm(fρ, ρ), fdm(fϕ, ϕ))
            @test isapprox(dρ, fdρ, rtol = 1e-12, atol = 2 * eps(Float64))
            @test isapprox(dϕ, fdϕ, rtol = 1e-12, atol = 2 * eps(Float64))
        end

        surf = TestData.gridsagsurfacebicubiccheby()
        fdm = central_fdm(10, 1)
        # not valid at the very boundary of the surface as stuff gets weird outside |u| < 1 and |v| < 1
        for u in -0.3:0.01:0.3, v in -0.15:0.01:0.15
            (du, dv) = partials(surf, u, v)
            fu(u) = point(surf, u, v)
            fv(v) = point(surf, u, v)
            # compute accurate finite difference approximation to derivative
            fdu, fdv = (fdm(fu, u), fdm(fv, v))
            @test isapprox(du, fdu, rtol = 1e-12, atol = 2 * eps(Float64))
            @test isapprox(dv, fdv, rtol = 1e-12, atol = 2 * eps(Float64))
        end
    end # testset GridSag

    @testset "ChebyshevSurface" begin
        surf = TestData.chebyshevsurface1() # with normradius
        fdm = central_fdm(10, 1)
        # not valid at the very boundary of the surface as stuff gets weird outside |u| < 1 and |v| < 1
        for u in -0.99:0.1:0.99, v in -0.99:0.1:0.99
            (du, dv) = partials(surf, u, v)
            fu(u) = point(surf, u, v)
            fv(v) = point(surf, u, v)
            # compute accurate finite difference approximation to derivative
            fdu, fdv = (fdm(fu, u), fdm(fv, v))
            @test isapprox(du, fdu, rtol = 1e-11, atol = 2 * eps(Float64)) #changed rtol from 1e-12 to 1e-11. FiniteDifferences approximation to the derivative had larger than expected error.
            @test isapprox(dv, fdv, rtol = 1e-11, atol = 2 * eps(Float64)) #changed rtol from 1e-12 to 1e-11. FiniteDifferences approximation to the derivative had larger than expected error.
        end
    end # testset ChebyshevSurface

end # testset SurfaceDefs

@otestset "Intersection" begin
    function samplepoints(numsamples, lowu, highu, lowv, highv)
        samples = Array{Tuple{Float64,Float64},1}(undef, 0)

        for i in 0:numsamples
            u = lowu * i / numsamples + (1 - i / numsamples) * highu
            for j in 0:numsamples
                v = lowv * j / numsamples + (1 - j / numsamples) * highv
                push!(samples, (u, v))
            end
        end
        return samples
    end

    @testset "Cylinder" begin
        # Random.seed!(SEED)

        # # random point intersections
        # ur, vr = uvrange(Cylinder)
        # for i in 1:10000
        #     cyl = Cylinder(rand() * 3.0, rand() * 50.0)
        #     u, v = rand() * (ur[2] - ur[1]) + ur[1], rand() * (vr[2] - vr[1]) + vr[1]
        #     pointon = point(cyl, u, v)
        #     origin = SVector(4.0, 4.0, 4.0)
        #     r = Ray(origin, pointon .- origin)
        #     halfspace = surfaceintersection(cyl, r)
        #     @test halfspace !== nothing && !((lower(halfspace) isa RayOrigin) || (upper(halfspace) isa Infinity))
        #     @test isapprox(pointon, point(lower(halfspace)), rtol = RTOLERANCE, atol = ATOLERANCE) || isapprox(pointon, point(upper(halfspace)), rtol = RTOLERANCE, atol = ATOLERANCE)
        # end

        # on xy
        cyl = Cylinder(0.5)
        r = Ray([1.0, 1.0, 0.0], [-1.0, -1.0, 0.0])
        intscts = surfaceintersection(cyl, r)
        xy = sqrt(2) / 4.0
        pt1 = [xy, xy, 0.0]
        pt2 = [-xy, -xy, 0.0]
        cpt1 = point(lower(intscts))
        cpt2 = point(upper(intscts))
        @test isapprox(pt1, cpt1, rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(pt2, cpt2, rtol = RTOLERANCE, atol = ATOLERANCE)

        # starting inside
        r = Ray([0.0, 0.0, 0.0], [1.0, 1.0, 0.0])
        intscts = surfaceintersection(cyl, r)
        @test lower(intscts) isa RayOrigin && isapprox(pt1, point(upper(intscts)), rtol = RTOLERANCE, atol = ATOLERANCE)

        # inside infinite
        r = Ray([0.0, 0.0, 0.0], [0.0, 0.0, 1.0])
        intscts = surfaceintersection(cyl, r)
        @test lower(intscts) isa RayOrigin && upper(intscts) isa Infinity

        # on surface
        r = Ray([0.5, 0.0, 0.0], [0.0, 0.0, 1.0])
        intscts = surfaceintersection(cyl, r)
        @test lower(intscts) isa RayOrigin && upper(intscts) isa Infinity

        # on diagonal
        r = Ray([1.0, 1.0, 1.0], [-1.0, -1.0, -1.0])
        intscts = surfaceintersection(cyl, r)
        cpt1 = point(lower(intscts))
        cpt2 = point(upper(intscts))
        pt1 = [xy, xy, xy]
        pt2 = [-xy, -xy, -xy]
        @test isapprox(pt1, cpt1, rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(pt2, cpt2, rtol = RTOLERANCE, atol = ATOLERANCE)

        # missing
        r = Ray([5.0, 5.0, 5.0], [0.0, 0.0, -1.0])
        intscts = surfaceintersection(cyl, r)
        @test intscts isa EmptyInterval
    end # testset cylinder

    @testset "Sphere" begin
        # Random.seed!(SEED)

        # # random point intersections
        # ur, vr = uvrange(Sphere)
        # for i in 1:10000
        #     sph = Sphere(rand() * 3)
        #     u, v = rand() * (ur[2] - ur[1]) + ur[1], rand() * (vr[2] - vr[1]) + vr[1]
        #     pointon = point(sph, u, v)
        #     origin = SVector(4.0, 4.0, 4.0)
        #     r = Ray(origin, pointon .- origin)
        #     halfspace = surfaceintersection(sph, r)
        #     @test halfspace !== nothing && !((lower(halfspace) isa RayOrigin) || (upper(halfspace) isa Infinity))
        #     @test isapprox(pointon, point(lower(halfspace)), rtol = RTOLERANCE, atol = ATOLERANCE) || isapprox(pointon, point(upper(halfspace)), rtol = RTOLERANCE, atol = ATOLERANCE)
        # end

        # on xy plane
        sph = Sphere(0.5)
        r = Ray([1.0, 1.0, 0.0], [-1.0, -1.0, 0.0])
        intscts = surfaceintersection(sph, r)
        xy = sqrt(2) / 4.0
        pt1 = [xy, xy, 0.0]
        pt2 = [-xy, -xy, 0.0]
        cpt1 = point(lower(intscts))
        cpt2 = point(upper(intscts))
        @test isapprox(pt1, cpt1, rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(pt2, cpt2, rtol = RTOLERANCE, atol = ATOLERANCE)

        # starting inside
        r = Ray([0.0, 0.0, 0.0], [1.0, 1.0, 0.0])
        intscts = surfaceintersection(sph, r)
        @test lower(intscts) isa RayOrigin && isapprox(pt1, point(upper(intscts)), rtol = RTOLERANCE, atol = ATOLERANCE)

        # on diagonal
        r = Ray([1.0, 1.0, 1.0], [-1.0, -1.0, -1.0])
        intscts = surfaceintersection(sph, r)
        cpt1 = point(lower(intscts))
        cpt2 = point(upper(intscts))
        xyz = sqrt(3) / 6.0
        pt1 = [xyz, xyz, xyz]
        pt2 = [-xyz, -xyz, -xyz]
        @test isapprox(pt1, cpt1, rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(pt2, cpt2, rtol = RTOLERANCE, atol = ATOLERANCE)

        # missing
        r = Ray([5.0, 5.0, 5.0], [0.0, 0.0, -1.0])
        intscts = surfaceintersection(sph, r)
        @test intscts isa EmptyInterval

        # tangent
        r = Ray([0.5, 0.5, 0.0], [0.0, -1.0, 0.0])
        intscts = surfaceintersection(sph, r)
        @test intscts isa EmptyInterval
    end # testset sphere

    @testset "Spherical Cap" begin
        @test_throws AssertionError SphericalCap(0.0, 1.0)
        @test_throws AssertionError SphericalCap(1.0, 0.0)

        sph = SphericalCap(0.5, π / 3, SVector(1.0, 1.0, 0.0), SVector(1.0, 1.0, 1.0))

        r = Ray([10.0, 10.0, 1.0], [-1.0, -1.0, 0.0])
        int = surfaceintersection(sph, r)
        @test isapprox(point(lower(int)), [1.0, 1.0, 1.0], rtol = RTOLERANCE, atol = ATOLERANCE) && upper(int) isa Infinity

        r = Ray([0.8, 1.4, 1.2], [1.0, -1.0, 0.0])
        int = surfaceintersection(sph, r)
        @test isapprox(point(lower(int)), [0.898231126982741, 1.301768873017259, 1.2], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(int)), [1.301768873017259, 0.8982311269827411, 1.2], rtol = RTOLERANCE, atol = ATOLERANCE)

        r = Ray([2.0, 2.0, 1.5], [-1.0, -1.0, 0.0])
        @test surfaceintersection(sph, r) isa EmptyInterval
    end # testset spherical cap

    @testset "Triangle" begin
        Random.seed!(SEED)

        tri = Triangle(SVector{3}(2.0, 1.0, 1.0), SVector{3}(1.0, 2.0, 1.0), SVector{3}(1.0, 1.0, 2.0), SVector{2}(0.0, 0.0), SVector{2}(1.0, 0.0), SVector{2}(0.5, 0.5))

        # random point intersections
        for i in 1:1000
            barycentriccoords = rand(3)
            barycentriccoords /= sum(barycentriccoords)
            pointontri = point(tri, barycentriccoords...)
            origin = SVector(3.0, 3.0, 3.0)
            r = Ray(origin, pointontri .- origin)
            halfspace = surfaceintersection(tri, r)
            @test !(halfspace isa EmptyInterval) && upper(halfspace) isa Infinity
            t = α(lower(halfspace))
            @test isapprox(point(r, t), point(lower(halfspace)), rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(pointontri, point(lower(halfspace)), rtol = RTOLERANCE, atol = ATOLERANCE)
        end

        # front and back face intersections
        tri = Triangle(SVector{3}(0.0, 0.0, 0.0), SVector{3}(1.0, 0.0, 0.0), SVector{3}(0.0, 1.0, 0.0), SVector{2}(0.0, 0.0), SVector{2}(1.0, 0.0), SVector{2}(0.0, 1.0))
        r1 = Ray([0.1, 0.1, 1.0], [0.0, 0.0, -1.0])
        r2 = Ray([0.1, 0.1, -1.0], [0.0, 0.0, 1.0])

        intsct1 = halfspaceintersection(surfaceintersection(tri, r1))
        intsct2 = halfspaceintersection(surfaceintersection(tri, r2))

        @test isapprox(point(intsct1), point(intsct2), rtol = RTOLERANCE, atol = ATOLERANCE)

        # ray should miss
        r = Ray([1.0, 1.0, 1.0], [0.0, 0.0, -1.0])
        intsct = surfaceintersection(tri, r)
        @test intsct isa EmptyInterval
    end # testset triangle

    @testset "Plane" begin
        rinplane = Ray([0.0, 0.0, 1.0], [1.0, 1.0, 0.0])
        routside = Ray([0.0, 0.0, 2.0], [1.0, 1.0, 1.0])
        rinside = Ray([0.0, 0.0, -1.0], [1.0, 1.0, -1.0])
        routintersects = Ray([0.0, 0.0, 2.0], [0.0, 0.0, -1.0])
        rinintersects = Ray([0.0, 0.0, -1.0], [0.0, 0.0, 1.0])
        pln = Plane(0.0, 0.0, 1.0, 0.0, 0.0, 1.0)

        # parallel to plane and outside
        @test surfaceintersection(pln, routside) isa EmptyInterval

        # NOTE for coplanar faces to work with visualization, rays in the plane must count as being 'inside'
        # parallel and on plane
        res = surfaceintersection(pln, rinplane)
        @test lower(res) isa RayOrigin && upper(res) isa Infinity

        # inside but not hitting
        res = surfaceintersection(pln, rinside)
        @test lower(res) isa RayOrigin && upper(res) isa Infinity

        # starts outside and hits
        res = surfaceintersection(pln, routintersects)
        @test isapprox(point(lower(res)), [0.0, 0.0, 1.0], rtol = RTOLERANCE, atol = ATOLERANCE) && upper(res) isa Infinity

        # starts inside and hits
        res = surfaceintersection(pln, rinintersects)
        @test lower(res) isa RayOrigin && isapprox(point(upper(res)), [0.0, 0.0, 1.0], rtol = RTOLERANCE, atol = ATOLERANCE)
    end # testset plane

    @testset "Rectangle" begin
        @test_throws AssertionError Rectangle(0.0, 1.0)
        @test_throws AssertionError Rectangle(1.0, 0.0)

        rinplane = Ray([0.0, 0.0, 0.0], [1.0, 1.0, 0.0])
        routside = Ray([0.0, 0.0, 1.0], [1.0, 1.0, 1.0])
        rinside = Ray([0.0, 0.0, -1.0], [1.0, 1.0, -1.0])
        routintersects = Ray([0.0, 0.0, 1.0], [0.0, 0.0, -1.0])
        rinintersects = Ray([0.0, 0.0, -1.0], [0.0, 0.0, 1.0])
        routmiss = Ray([0.0, 2.0, 1.0], [0.0, 0.0, -1.0])
        rinmiss = Ray([0.0, 2.0, -1.0], [0.0, 0.0, 1.0])
        rinbounds = Ray([0.5, 0.5, -1.0], [0.0, 0.0, 1.0])
        routbounds = Ray([0.5, 0.5, 1.0], [0.0, 0.0, -1.0])

        rect = Rectangle(0.5, 0.5)

        # parallel to plane and outside
        @test surfaceintersection(rect, routside) isa EmptyInterval

        # parallel and on plane
        @test surfaceintersection(rect, rinplane) isa EmptyInterval

        # inside but not hitting
        @test surfaceintersection(rect, rinside) isa EmptyInterval

        # starts outside and hits
        res = surfaceintersection(rect, routintersects)
        @test isapprox(point(lower(res)), [0.0, 0.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE) && upper(res) isa Infinity

        # starts inside and hits
        res = surfaceintersection(rect, rinintersects)
        @test lower(res) isa RayOrigin && isapprox(point(upper(res)), [0.0, 0.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE)

        # starts inside and misses bounds
        @test surfaceintersection(rect, rinmiss) isa EmptyInterval

        # starts outside and misses bounds
        @test surfaceintersection(rect, routmiss) isa EmptyInterval

        # starts outside and hits bounds
        res = surfaceintersection(rect, routbounds)
        @test isapprox(point(lower(res)), [0.5, 0.5, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE) && upper(res) isa Infinity

        # starts inside and hits bounds
        res = surfaceintersection(rect, rinbounds)
        @test lower(res) isa RayOrigin && isapprox(point(upper(res)), [0.5, 0.5, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE)

        # test a rectangle with translation and rotation
        rect2 = Rectangle(0.3, 0.5, SVector(3.0, 1.0, 4.0), SVector(0.2, 0.3, 0.4))
        rayhit = Ray([0.6, 0.6, 0.6], [-0.3, -0.1, -0.3])
        raymiss = Ray([0.7, 0.4, 0.4], [-0.3, -0.1, -0.3])
        res = surfaceintersection(rect2, rayhit)
        @test isapprox(point(lower(res)), [0.2863636363636363, 0.4954545454545454, 0.2863636363636363], rtol = RTOLERANCE, atol = ATOLERANCE) && upper(res) isa Infinity
        @test surfaceintersection(rect2, raymiss) isa EmptyInterval
    end # testset Rectangle

    @testset "Ellipse" begin
        @test_throws AssertionError Ellipse(0.0, 1.0)
        @test_nowarn Circle(0.5)

        rinplane = Ray([0.0, 0.0, 0.0], [1.0, 1.0, 0.0])
        routside = Ray([0.0, 0.0, 1.0], [1.0, 1.0, 1.0])
        rinside = Ray([0.0, 0.0, -1.0], [1.0, 1.0, -1.0])
        routintersects = Ray([0.0, 0.0, 1.0], [0.0, 0.0, -1.0])
        rinintersects = Ray([0.0, 0.0, -1.0], [0.0, 0.0, 1.0])
        routmiss = Ray([0.0, 2.0, 1.0], [0.0, 0.0, -1.0])
        rinmiss = Ray([0.0, 2.0, -1.0], [0.0, 0.0, 1.0])
        rinbounds = Ray([0.5, 0.0, -1.0], [0.0, 0.0, 1.0])
        routbounds = Ray([0.5, 0.0, 1.0], [0.0, 0.0, -1.0])
        rasymhit = Ray([0.0, 0.9, 1.0], [0.0, 0.0, -1.0])
        rasymmiss = Ray([0.9, 0.0, 1.0], [0.0, 0.0, -1.0])

        ell = Ellipse(0.5, 1.0)

        # parallel to plane and outside
        @test surfaceintersection(ell, routside) isa EmptyInterval

        # parallel and on plane
        @test surfaceintersection(ell, rinplane) isa EmptyInterval

        # inside but not hitting
        @test surfaceintersection(ell, rinside) isa EmptyInterval

        # starts outside and hits
        res = surfaceintersection(ell, routintersects)
        @test isapprox(point(lower(res)), [0.0, 0.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE) && upper(res) isa Infinity

        # starts inside and hits
        res = surfaceintersection(ell, rinintersects)
        @test lower(res) isa RayOrigin && isapprox(point(upper(res)), [0.0, 0.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE)

        # starts inside and misses bounds
        @test surfaceintersection(ell, rinmiss) isa EmptyInterval

        # starts outside and misses bounds
        @test surfaceintersection(ell, routmiss) isa EmptyInterval

        # starts outside and hits bounds
        res = surfaceintersection(ell, routbounds)
        @test isapprox(point(lower(res)), [0.5, 0.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE) && upper(res) isa Infinity

        # starts inside and hits bounds
        res = surfaceintersection(ell, rinbounds)
        @test lower(res) isa RayOrigin && isapprox(point(upper(res)), [0.5, 0.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE)

        # asymmetric hit
        res = surfaceintersection(ell, rasymhit)
        @test isapprox(point(lower(res)), [0.0, 0.9, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE) && upper(res) isa Infinity

        # asymmetric miss
        @test surfaceintersection(ell, rasymmiss) isa EmptyInterval

        # test an ellipse with translation and rotation
        ell2 = Ellipse(0.3, 0.5, SVector(3.0, 1.0, 4.0), SVector(0.2, 0.3, 0.4))
        rayhit = Ray([0.6, 0.6, 0.6], [-0.3, -0.1, -0.3])
        raymiss = Ray([0.7, 0.4, 0.4], [-0.3, -0.1, -0.3])
        res = surfaceintersection(ell2, rayhit)
        @test isapprox(point(lower(res)), [0.2863636363636363, 0.4954545454545454, 0.2863636363636363], rtol = RTOLERANCE, atol = ATOLERANCE) && upper(res) isa Infinity
        @test surfaceintersection(ell2, raymiss) isa EmptyInterval
    end # testset Ellipse

    @testset "Hexagon" begin
        @test_nowarn Hexagon(0.5)

        rinplane = Ray([0.0, 0.0, 0.0], [1.0, 1.0, 0.0])
        routside = Ray([0.0, 0.0, 1.0], [1.0, 1.0, 1.0])
        rinside = Ray([0.0, 0.0, -1.0], [1.0, 1.0, -1.0])
        routintersects = Ray([0.0, 0.0, 1.0], [0.0, 0.0, -1.0])
        rinintersects = Ray([0.0, 0.0, -1.0], [0.0, 0.0, 1.0])
        routmiss = Ray([0.0, 2.0, 1.0], [0.0, 0.0, -1.0])
        rinmiss = Ray([0.0, 2.0, -1.0], [0.0, 0.0, 1.0])
        rinbounds = Ray([sqrt(3) / 2, 0.0, -1.0], [0.0, 0.0, 1.0])
        routbounds = Ray([sqrt(3) / 2, 0.0, 1.0], [0.0, 0.0, -1.0])
        rasymhit = Ray([0.0, 1.0, 1.0], [0.0, 0.0, -1.0])
        rasymmiss = Ray([1.0, 0.0, 1.0], [0.0, 0.0, -1.0])

        hex = Hexagon(1.0)

        # parallel to plane and outside
        @test surfaceintersection(hex, routside) isa EmptyInterval

        # parallel and on plane
        @test surfaceintersection(hex, rinplane) isa EmptyInterval

        # inside but not hitting
        @test surfaceintersection(hex, rinside) isa EmptyInterval

        # starts outside and hits
        res = surfaceintersection(hex, routintersects)
        @test isapprox(point(lower(res)), [0.0, 0.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE) && upper(res) isa Infinity

        # starts inside and hits
        res = surfaceintersection(hex, rinintersects)
        @test lower(res) isa RayOrigin && isapprox(point(upper(res)), [0.0, 0.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE)

        # starts inside and misses bounds
        @test surfaceintersection(hex, rinmiss) isa EmptyInterval

        # starts outside and misses bounds
        @test surfaceintersection(hex, routmiss) isa EmptyInterval

        # starts outside and hits bounds
        res = surfaceintersection(hex, routbounds)
        @test isapprox(point(lower(res)), [sqrt(3) / 2, 0.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE) && upper(res) isa Infinity

        # starts inside and hits bounds
        res = surfaceintersection(hex, rinbounds)
        @test lower(res) isa RayOrigin && isapprox(point(upper(res)), [sqrt(3) / 2, 0.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE)

        # asymmetric hit
        res = surfaceintersection(hex, rasymhit)
        @test isapprox(point(lower(res)), [0.0, 1.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE) && upper(res) isa Infinity

        # asymmetric miss
        @test surfaceintersection(hex, rasymmiss) isa EmptyInterval

        # test an ellipse with translation and rotation
        hex2 = Hexagon(0.4, SVector(3.0, 1.0, 4.0), SVector(0.2, 0.3, 0.4))
        rayhit = Ray([0.6, 0.6, 0.6], [-0.3, -0.1, -0.3])
        raymiss = Ray([0.7, 0.4, 0.4], [-0.3, -0.1, -0.3])
        res = surfaceintersection(hex2, rayhit)
        @test isapprox(point(lower(res)), [0.2863636363636363, 0.4954545454545454, 0.2863636363636363], rtol = RTOLERANCE, atol = ATOLERANCE) && upper(res) isa Infinity
        @test surfaceintersection(hex2, raymiss) isa EmptyInterval
    end # testset Hexagon

    @testset "Stops" begin
        infiniterect = RectangularAperture(0.4, 0.8, SVector(0.0, 1.0, 0.0), SVector(0.0, 0.0, 1.0))
        finiterect = RectangularAperture(0.4, 0.8, 1.0, 2.0, SVector(0.0, 1.0, 0.0), SVector(0.0, 0.0, 1.0))
        # through hole
        r = Ray([0.0, 1.0, 1.0], [0.0, -1.0, 0.0])
        @test surfaceintersection(infiniterect, r) isa EmptyInterval
        @test surfaceintersection(finiterect, r) isa EmptyInterval
        r = Ray([0.0, -1.0, 1.0], [0.0, 1.0, 0.0])
        @test surfaceintersection(infiniterect, r) isa EmptyInterval
        @test surfaceintersection(finiterect, r) isa EmptyInterval
        # on edge
        r = Ray([0.0, 1.0, 0.6], [0.0, -1.0, 0.0])
        res = surfaceintersection(infiniterect, r)
        @test isapprox(point(lower(res)), [0.0, 0.0, 0.6], rtol = RTOLERANCE, atol = ATOLERANCE) && upper(res) isa Infinity
        res = surfaceintersection(finiterect, r)
        @test isapprox(point(lower(res)), [0.0, 0.0, 0.6], rtol = RTOLERANCE, atol = ATOLERANCE) && upper(res) isa Infinity
        # through hole asym
        r = Ray([0.7, 1.0, 1.0], [0.0, -1.0, 0.0])
        @test surfaceintersection(infiniterect, r) isa EmptyInterval
        @test surfaceintersection(finiterect, r) isa EmptyInterval
        # on finite edge
        r = Ray([1.0, -1.0, 2.0], [0.0, 1.0, 0.0])
        res = surfaceintersection(infiniterect, r)
        @test lower(res) isa RayOrigin && isapprox(point(upper(res)), [1.0, 0.0, 2.0], rtol = RTOLERANCE, atol = ATOLERANCE)
        res = surfaceintersection(finiterect, r)
        @test lower(res) isa RayOrigin && isapprox(point(upper(res)), [1.0, 0.0, 2.0], rtol = RTOLERANCE, atol = ATOLERANCE)
        # outside finite
        r = Ray([2.0, 1.0, 3.0], [0.0, -1.0, 0.0])
        @test isapprox(point(surfaceintersection(infiniterect, r)), [2.0, 0.0, 3.0], rtol = RTOLERANCE, atol = ATOLERANCE)
        @test surfaceintersection(finiterect, r) isa EmptyInterval

        infinitecirc = CircularAperture(0.4, SVector(0.0, 1.0, 0.0), SVector(0.0, 0.0, 1.0))
        finitecirc = CircularAperture(0.4, 1.0, 2.0, SVector(0.0, 1.0, 0.0), SVector(0.0, 0.0, 1.0))
        # through hole
        r = Ray([0.0, 1.0, 1.0], [0.0, -1.0, 0.0])
        @test surfaceintersection(infinitecirc, r) isa EmptyInterval
        @test surfaceintersection(finitecirc, r) isa EmptyInterval
        r = Ray([0.0, -1.0, 1.0], [0.0, 1.0, 0.0])
        @test surfaceintersection(infinitecirc, r) isa EmptyInterval
        @test surfaceintersection(finitecirc, r) isa EmptyInterval
        # on edge
        r = Ray([0.0, 1.0, 0.6], [0.0, -1.0, 0.0])
        res = surfaceintersection(infinitecirc, r)
        @test isapprox(point(lower(res)), [0.0, 0.0, 0.6], rtol = RTOLERANCE, atol = ATOLERANCE) && upper(res) isa Infinity
        res = surfaceintersection(finitecirc, r)
        @test isapprox(point(lower(res)), [0.0, 0.0, 0.6], rtol = RTOLERANCE, atol = ATOLERANCE) && upper(res) isa Infinity
        # through hole 2
        r = Ray([0.3, 1.0, 1.0], [0.0, -1.0, 0.0])
        @test surfaceintersection(infinitecirc, r) isa EmptyInterval
        @test surfaceintersection(finitecirc, r) isa EmptyInterval
        # on finite edge
        r = Ray([1.0, -1.0, 2.0], [0.0, 1.0, 0.0])
        res = surfaceintersection(infinitecirc, r)
        @test lower(res) isa RayOrigin && isapprox(point(upper(res)), [1.0, 0.0, 2.0], rtol = RTOLERANCE, atol = ATOLERANCE)
        res = surfaceintersection(finitecirc, r)
        @test lower(res) isa RayOrigin && isapprox(point(upper(res)), [1.0, 0.0, 2.0], rtol = RTOLERANCE, atol = ATOLERANCE)
        # outside finite
        r = Ray([2.0, 1.0, 3.0], [0.0, -1.0, 0.0])
        @test isapprox(point(surfaceintersection(infinitecirc, r)), [2.0, 0.0, 3.0], rtol = RTOLERANCE, atol = ATOLERANCE)
        @test surfaceintersection(finitecirc, r) isa EmptyInterval

        annulus = Annulus(0.5, 1.0, SVector(0.0, 1.0, 0.0), SVector(0.0, 0.0, 1.0))
        # through hole
        r = Ray([0.0, 1.0, 1.0], [0.0, -1.0, 0.0])
        @test surfaceintersection(annulus, r) isa EmptyInterval
        r = Ray([0.0, -1.0, 1.0], [0.0, 1.0, 0.0])
        @test surfaceintersection(annulus, r) isa EmptyInterval
        # on edge
        r = Ray([0.0, 1.0, 0.5], [0.0, -1.0, 0.0])
        res = surfaceintersection(annulus, r)
        @test isapprox(point(lower(res)), [0.0, 0.0, 0.5], rtol = RTOLERANCE, atol = ATOLERANCE) && upper(res) isa Infinity
        # through hole 2
        r = Ray([0.3, 1.0, 1.0], [0.0, -1.0, 0.0])
        @test surfaceintersection(infiniterect, r) isa EmptyInterval
        @test surfaceintersection(annulus, r) isa EmptyInterval
        # on finite edge
        r = Ray([1.0, -1.0, 1.0], [0.0, 1.0, 0.0])
        res = surfaceintersection(annulus, r)
        @test lower(res) isa RayOrigin && isapprox(point(upper(res)), [1.0, 0.0, 1.0], rtol = RTOLERANCE, atol = ATOLERANCE)
        # outside finite
        r = Ray([0.9, 1.0, 1.9], [0.0, -1.0, 0.0])
        @test surfaceintersection(annulus, r) isa EmptyInterval
    end # testset Stops

    @testset "Bezier" begin
        surf = TestData.beziersurface()
        accelsurf = AcceleratedParametricSurface(surf)
        numsamples = 100

        samples = samplepoints(numsamples, 0.05, 0.95, 0.05, 0.95)
        missedintersections = 0

        for uv in samples
            pt = point(surf, uv[1], uv[2])
            origin = [0.5, 0.5, 5.0]
            r = Ray(origin, pt .- origin)
            allintersections = surfaceintersection(accelsurf, r)
            if allintersections isa EmptyInterval
                # sampled triangle acceleration structure doesn't guarantee all intersecting rays will be detected.
                # Have to use patch subdivision and convex hull polyhedron to do this.
                missedintersections += 1
            else
                if isa(allintersections, Interval)
                    allintersections = (allintersections,)
                end

                for intsct in allintersections
                    @test isapprox(point(halfspaceintersection(intsct)), pt, rtol = RTOLERANCE, atol = ATOLERANCE)
                    @test upper(intsct) isa Infinity
                end
            end
        end

        if missedintersections > 0
            @warn "$missedintersections out of total $(length(samples)) bezier suface intersections were missed"
        end

        # miss from outside
        r = Ray([0.5, 0.5, 5.0], [0.0, 0.0, 1.0])
        @test surfaceintersection(accelsurf, r) isa EmptyInterval
        r = Ray([5.0, 0.5, -5.0], [0.0, 0.0, -1.0])
        @test surfaceintersection(accelsurf, r) isa EmptyInterval

        # miss from inside
        r = Ray([0.5, 0.5, -5.0], [0.0, 0.0, -1.0])
        res = surfaceintersection(accelsurf, r)
        # TODO!! Fix bezier surface to create halfspace
        # @test lower(res) isa RayOrigin && upper(res) isa Infinity

        # hit from inside
        r = Ray([0.5, 0.5, 0.2], [0.0, 1.0, 0.0])
        res = surfaceintersection(accelsurf, r)
        @test lower(res) isa RayOrigin && isapprox(point(upper(res)), [0.5, 0.8996900152986361, 0.2], rtol = RTOLERANCE, atol = ATOLERANCE)

        # hit from outside
        r = Ray([0.0, 0.5, 0.2], [1.0, 0.0, -1.0])
        res = surfaceintersection(accelsurf, r)
        @test isapprox(point(lower(res)), [0.06398711204047353, 0.5, 0.13601288795952649], rtol = RTOLERANCE, atol = ATOLERANCE) && upper(res) isa Infinity

        # two hits from outside
        r = Ray([0.0, 0.5, 0.25], [1.0, 0.0, 0.0])
        res = surfaceintersection(accelsurf, r)
        @test isapprox(point(lower(res)), [0.12607777053030267, 0.5, 0.25], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(res)), [0.8705887077060419, 0.5, 0.25], rtol = RTOLERANCE, atol = ATOLERANCE)

        surf = TestData.wavybeziersurface()
        accelsurf = AcceleratedParametricSurface(surf)

        # 3 hit starting outside
        r = Ray([0.5, 0.0, 0.0], [0.0, 1.0, 0.0])
        res = surfaceintersection(accelsurf, r)
        @test isa(res, DisjointUnion) && length(res) == 2 && isapprox(point(lower(res[1])), [0.5, 0.10013694786182059, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(res[1])), [0.5, 0.49625, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(lower(res[2])), [0.5, 0.8971357794109067, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE) && (upper(res[2]) isa Infinity)

        # 3 hit starting inside
        r = Ray([0.5, 1.0, 0.0], [0.0, -1.0, 0.0])
        res = surfaceintersection(accelsurf, r)
        @test isa(res, DisjointUnion) && length(res) == 2 && (lower(res[1]) isa RayOrigin) && isapprox(point(upper(res[1])), [0.5, 0.8971357794109067, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(lower(res[2])), [0.5, 0.49625, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(res[2])), [0.5, 0.10013694786182059, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE)

        surf = TestData.verywavybeziersurface()
        accelsurf = AcceleratedParametricSurface(surf, 20)

        # five hits starting outside
        r = Ray([0.9, 0.0, -0.3], [0.0, 1.0, 0.7])
        res = surfaceintersection(accelsurf, r)
        a = isapprox(point(lower(res[1])), [0.9, 0.03172286522032046, -0.2777939943457758], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(res[1])), [0.9, 0.1733979947040411, -0.17862140370717122], rtol = RTOLERANCE, atol = ATOLERANCE)
        b = isapprox(point(lower(res[2])), [0.9, 0.5335760974594397, 0.07350326822160776], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(res[2])), [0.9, 0.7767707607392784, 0.24373953251749486], rtol = RTOLERANCE, atol = ATOLERANCE)
        c = isapprox(point(lower(res[3])), [0.9, 0.9830891958374246, 0.3881624370861975], rtol = RTOLERANCE, atol = ATOLERANCE) && (upper(res[3]) isa Infinity)
        @test isa(res, DisjointUnion) && length(res) == 3 && a && b && c

        # five hits starting inside
        r = Ray([0.9, 1.0, 0.4], [0.0, -1.0, -0.7])
        res = surfaceintersection(accelsurf, r)
        a = (lower(res[1]) isa RayOrigin) && isapprox(point(upper(res[1])), [0.9, 0.9830891958374246, 0.3881624370861975], rtol = RTOLERANCE, atol = ATOLERANCE)
        b = isapprox(point(lower(res[2])), [0.9, 0.7767707607392784, 0.24373953251749486], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(res[2])), [0.9, 0.5335760974594397, 0.07350326822160776], rtol = RTOLERANCE, atol = ATOLERANCE)
        c = isapprox(point(lower(res[3])), [0.9, 0.17339799470404108, -0.1786214037071712], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(res[3])), [0.9, 0.03172286522032046, -0.27779399434577573], rtol = RTOLERANCE, atol = ATOLERANCE)
        @test isa(res, DisjointUnion) && length(res) == 3 && a && b && c

        # 4 hits starting inside
        r = Ray([0.1, 0.0, -0.3], [0.0, 1.0, 0.7])
        res = surfaceintersection(accelsurf, r)
        a = (lower(res[1]) isa RayOrigin) && isapprox(point(upper(res[1])), [0.1, 0.2851860296285551, -0.10036977926001144], rtol = RTOLERANCE, atol = ATOLERANCE)
        b = isapprox(point(lower(res[2])), [0.1, 0.5166793625025807, 0.06167555375180668], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(res[2])), [0.1, 0.7770862508789854, 0.24396037561528983], rtol = RTOLERANCE, atol = ATOLERANCE)
        c = isapprox(point(lower(res[3])), [0.1, 0.98308919558696, 0.3881624369108719], rtol = RTOLERANCE, atol = ATOLERANCE) && (upper(res[3]) isa Infinity)
        @test isa(res, DisjointUnion) && length(res) == 3 && a && b && c

        # 4 hits starting outside
        r = Ray([0.9, 0.9, 0.4], [0.0, -0.9, -0.7])
        res = surfaceintersection(accelsurf, r)
        a = isapprox(point(lower(res[1])), [0.9, 0.736072142615238, 0.2725005553674076], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(res[1])), [0.9, 0.567439326091764, 0.141341698071372], rtol = RTOLERANCE, atol = ATOLERANCE)
        b = isapprox(point(lower(res[2])), [0.9, 0.16601081959179267, -0.1708804736508277], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(res[2])), [0.9, 0.032434058775915924, -0.274773509840954], rtol = RTOLERANCE, atol = ATOLERANCE)
        @test isa(res, DisjointUnion) && length(res) == 2 && a && b
    end # testset Bezier

    @testset "Zernike" begin
        @test_nowarn ZernikeSurface(1.5)
        @test_throws AssertionError ZernikeSurface(0.0)
        @test_throws AssertionError ZernikeSurface(1.5, aspherics = [(1, 0.1)])

        z1 = TestData.zernikesurface1()
        az1 = AcceleratedParametricSurface(z1, 20)
        numsamples = 100

        # test that an on axis ray doesn't miss
        @test !(surfaceintersection(AcceleratedParametricSurface(z1), Ray([0.0, 0.0, 10.0], [0.0, 0.0, -1.0])) isa EmptyInterval)

        samples = samplepoints(numsamples, 0.01, 0.99, 0.01, 2π - 0.01)
        missedintersections = 0

        for uv in samples
            pt = point(z1, uv[1], uv[2])
            origin = [0.0, 0.0, 5.0]
            r = Ray(origin, pt .- origin)
            allintersections = surfaceintersection(az1, r)
            if (allintersections isa EmptyInterval)
                # sampled triangle acceleration structure doesn't guarantee all intersecting rays will be detected.
                # Have to use patch subdivision and convex hull polyhedron to do this.
                missedintersections += 1
            else
                if isa(allintersections, Interval)
                    allintersections = (allintersections,)
                end

                for intsct in allintersections
                    if lower(intsct) isa RayOrigin
                        missedintersections += 1
                    else
                        # closest point should be on the surface, furthest should be on bounding prism
                        @test isapprox(point(lower(intsct)), pt, rtol = RTOLERANCE, atol = ATOLERANCE)
                        @test isa(upper(intsct), Intersection{Float64,3}) || (OpticSim.direction(r)[3] == -1 && isa(upper(intsct), Infinity{Float64}))
                    end
                end
            end
        end

        if missedintersections > 0
            @warn "$missedintersections out of total $(length(samples)) zernike suface intersections were missed"
        end

        # hit from inside
        r = Ray([0.0, 0.0, -0.5], [0.1, 0.2, 0.5])
        intvl = surfaceintersection(az1, r)
        @test (lower(intvl) isa RayOrigin) && isapprox(point(upper(intvl)), [0.12363711269619936, 0.24727422539239863, 0.11818556348099664], rtol = RTOLERANCE, atol = ATOLERANCE)

        # hit from outside
        r = Ray([0.0, 0.0, 0.5], [0.1, 0.2, -0.5])
        intvl = surfaceintersection(az1, r)
        @test isapprox(point(lower(intvl)), [0.07118311042701463, 0.1423662208540292, 0.144084447864927], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(intvl)), [0.6708203932499369, 1.3416407864998738, -2.8541019662496843], rtol = RTOLERANCE, atol = ATOLERANCE)

        # miss from inside
        r = Ray([0.2, 0.2, -0.5], [0.0, 0.0, -1.0])
        intvl = surfaceintersection(az1, r)
        @test (lower(intvl) isa RayOrigin) && (upper(intvl) isa Infinity)

        # miss from outside
        r = Ray([0.2, 2.0, -0.5], [0.0, 0.0, 1.0])
        @test surfaceintersection(az1, r) isa EmptyInterval
        r = Ray([0.2, 2.0, 0.5], [0.0, 0.0, -1.0])
        @test surfaceintersection(az1, r) isa EmptyInterval
        r = Ray([0.2, 0.2, 0.5], [0.0, 0.0, 1.0])
        @test surfaceintersection(az1, r) isa EmptyInterval
        r = Ray([0.2, 2.0, 2.0], [0.0, -1.0, 0.0])
        @test surfaceintersection(az1, r) isa EmptyInterval

        # two hits from outside
        r = Ray([2.0, 0.0, 0.25], [-1.0, 0.0, 0.0])
        hit = surfaceintersection(az1, r)
        a = isapprox(point(lower(hit[1])), [1.5, 0.0, 0.25], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(hit[1])), [1.3571758210851095, 0.0, 0.25], rtol = RTOLERANCE, atol = ATOLERANCE)
        b = isapprox(point(lower(hit[2])), [-1.2665828947521165, 0.0, 0.25], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(hit[2])), [-1.5, 0.0, 0.25], rtol = RTOLERANCE, atol = ATOLERANCE)
        @test isa(hit, DisjointUnion) && length(hit) == 2 && a && b

        # hit cyl from outside
        r = Ray([0.0, 2.0, -0.5], [0.0, -1.0, 0.0])
        intvl = surfaceintersection(az1, r)
        @test isapprox(point(lower(intvl)), [0.0, 1.5, -0.5], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(intvl)), [0.0, -1.5, -0.5], rtol = RTOLERANCE, atol = ATOLERANCE)

        # hit cyl from inside
        r = Ray([0.0, 0.0, -0.5], [0.0, -1.0, 0.0])
        intvl = surfaceintersection(az1, r)
        @test (lower(intvl) isa RayOrigin) && isapprox(point(upper(intvl)), [0.0, -1.5, -0.5], rtol = RTOLERANCE, atol = ATOLERANCE)

        z2 = TestData.zernikesurface2()
        az2 = AcceleratedParametricSurface(z2, 20)

        # two hits
        r = Ray([2.0, -1.0, 0.2], [-1.0, 0.0, 0.0])
        intvl = surfaceintersection(az2, r)
        @test isapprox(point(lower(intvl)), [0.238496383738369, -1.0, 0.2], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(intvl)), [-0.8790518227874484, -1.0, 0.2], rtol = RTOLERANCE, atol = ATOLERANCE)

        # three hits
        r = Ray([-0.7, 1.0, 0.2], [0.0, -1.0, 0.0])
        hit = surfaceintersection(az2, r)
        a = (lower(hit[1]) isa RayOrigin) && isapprox(point(upper(hit[1])), [-0.7, 0.8732489598020176, 0.2], rtol = RTOLERANCE, atol = ATOLERANCE)
        b = isapprox(point(lower(hit[2])), [-0.7, -0.34532502965048606, 0.2], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(hit[2])), [-0.7, -1.1615928003236047, 0.2], rtol = RTOLERANCE, atol = ATOLERANCE)
        @test isa(hit, DisjointUnion) && length(hit) == 2 && a && b

        # four hits
        r = Ray([1.5, 0.1, 0.35], [-1.0, -0.5, 0.0])
        hit = surfaceintersection(az2, r)
        a = isapprox(point(lower(hit[1])), [1.4371467631969683, 0.06857338159848421, 0.35], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(hit[1])), [1.1428338578611368, -0.07858307106943167, 0.35], rtol = RTOLERANCE, atol = ATOLERANCE)
        b = isapprox(point(lower(hit[2])), [0.12916355683491865, -0.5854182215825409, 0.35], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(hit[2])), [-0.7290713614371003, -1.0145356807185504, 0.35], rtol = RTOLERANCE, atol = ATOLERANCE)
        @test isa(hit, DisjointUnion) && length(hit) == 2 && a && b

        # failure case, had a bug where rounding error would cause o2 to fail
        o1 = leaf(AcceleratedParametricSurface(ZernikeSurface(1.4 * 1.15)), translation(0.0, 0.0, 3.0))()
        o2 = leaf(AcceleratedParametricSurface(ZernikeSurface(1.4 * 1.15)), translation(2.0, 0.0, 3.0))()
        r = Ray([-5.0, 0.0, 0.0], [1.0, 0.0, 0.0])
        @test !(surfaceintersection(o1, r) isa EmptyInterval)
        @test !(surfaceintersection(o2, r) isa EmptyInterval)
    end # testset Zernike

    @testset "QType" begin
        @test_nowarn QTypeSurface(1.5)
        @test_throws AssertionError QTypeSurface(0.0)

        q1 = TestData.qtypesurface1()
        aq1 = AcceleratedParametricSurface(q1, 20)
        numsamples = 100

        # test that an on axis ray doesn't miss
        @test !(surfaceintersection(AcceleratedParametricSurface(q1), Ray([0.0, 0.0, 10.0], [0.0, 0.0, -1.0])) isa EmptyInterval)

        samples = samplepoints(numsamples, 0.01, 0.99, 0.01, 2π - 0.01)
        missedintersections = 0

        for uv in samples
            pt = point(q1, uv[1], uv[2])
            origin = [0.0, 0.0, 5.0]
            r = Ray(origin, pt .- origin)
            allintersections = surfaceintersection(aq1, r)
            if (allintersections isa EmptyInterval)
                # sampled triangle acceleration structure doesn't guarantee all intersecting rays will be detected.
                # Have to use patch subdivision and convex hull polyhedron to do this.
                missedintersections += 1
            else
                if isa(allintersections, Interval)
                    allintersections = (allintersections,)
                end

                for intsct in allintersections
                    if lower(intsct) isa RayOrigin
                        missedintersections += 1
                    else
                        # closest point should be on the surface, furthest should be on bounding prism
                        @test isapprox(point(lower(intsct)), pt, rtol = RTOLERANCE, atol = ATOLERANCE)
                        @test isa(upper(intsct), Intersection{Float64,3}) || (OpticSim.direction(r)[3] == -1 && isa(upper(intsct), Infinity{Float64}))
                    end
                end
            end
        end

        if missedintersections > 0
            @warn "$missedintersections out of total $(length(samples)) qtype suface intersections were missed"
        end

        # hit from inside
        r = Ray([0.0, 0.0, -0.5], [0.1, 0.2, 0.5])
        intvl = surfaceintersection(aq1, r)
        @test (lower(intvl) isa RayOrigin) && isapprox(point(upper(intvl)), [0.09758258750208074, 0.19516517500416142, -0.01208706248959643], rtol = RTOLERANCE, atol = ATOLERANCE)

        # hit from outside
        r = Ray([0.0, 0.0, 0.5], [0.1, 0.2, -0.5])
        intvl = surfaceintersection(aq1, r)
        @test isapprox(point(lower(intvl)), [0.10265857811957124, 0.20531715623914257, -0.013292890597856405], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(intvl)), [0.6708203932499369, 1.3416407864998738, -2.8541019662496843], rtol = RTOLERANCE, atol = ATOLERANCE)

        # miss from inside
        r = Ray([0.2, 0.2, -0.5], [0.0, 0.0, -1.0])
        intvl = surfaceintersection(aq1, r)
        @test (lower(intvl) isa RayOrigin) && (upper(intvl) isa Infinity)

        # miss from outside
        r = Ray([0.2, 2.0, -0.5], [0.0, 0.0, 1.0])
        @test surfaceintersection(aq1, r) isa EmptyInterval
        r = Ray([0.2, 2.0, 0.5], [0.0, 0.0, -1.0])
        @test surfaceintersection(aq1, r) isa EmptyInterval
        r = Ray([0.2, 0.2, 0.5], [0.0, 0.0, 1.0])
        @test surfaceintersection(aq1, r) isa EmptyInterval
        r = Ray([0.2, 2.0, 2.0], [0.0, -1.0, 0.0])
        @test surfaceintersection(aq1, r) isa EmptyInterval
    end # testset QType

    @testset "BoundingBox" begin
        function boundingval(a::BoundingBox{T}, axis::Int, plane::Bool) where {T<:Real}
            if axis === 1
                return plane ? a.xmax : a.xmin
            elseif axis === 2
                return plane ? a.ymax : a.ymin
            elseif axis == 3
                return plane ? a.zmax : a.zmin
            else
                throw(ErrorException("Invalid axis: $axis"))
            end
        end

        Random.seed!(SEED)

        axis(x) = (x + 1) ÷ 2
        face(x) = mod(x, 2)
        facevalue(a::BoundingBox, x) = boundingval(a, axis(x), Bool(face(x)))
        boundingindices = ((2, 3), (1, 3), (1, 2))

        function facepoint(a::BoundingBox, facenumber)
            axisindex = axis(facenumber)
            indices = boundingindices[axisindex]
            b1min, b1max = boundingval(a, indices[1], false), boundingval(a, indices[1], true)
            b2min, b2max = boundingval(a, indices[2], false), boundingval(a, indices[2], true)

            step = 0.00001
            pt1val = rand((b1min + step):step:(b1max - step))
            pt2val = rand((b2min + step):step:(b2max - step))

            result = Array{Float64,1}(undef, 3)
            result[indices[1]] = pt1val
            result[indices[2]] = pt2val
            result[axisindex] = facevalue(a, facenumber)
            return result
        end

        bbox = BoundingBox(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0)

        for i in 1:10000
            # randomly pick two planes
            face1 = rand(1:6)
            face2 = rand(1:6)
            while face2 == face1
                face2 = rand(1:6)
            end
            # pick a point on each face
            pt1 = facepoint(bbox, face1)
            pt2 = facepoint(bbox, face2)

            direction = normalize(pt1 - pt2)
            origin = pt2 - 3 * direction

            r = Ray(origin, direction)

            @test doesintersect(bbox, r)
        end

        # should miss
        @test !doesintersect(bbox, Ray([-2.0, -2.0, 0.0], [1.0, 0.0, 0.0]))
        @test !doesintersect(bbox, Ray([-2.0, -2.0, 0.0], [0.0, 1.0, 0.0]))
        @test !doesintersect(bbox, Ray([-2.0, -2.0, 0.0], [0.0, 0.0, 1.0]))
        @test !doesintersect(bbox, Ray([-2.0, 0.0, 0.0], [-1.0, 0.0, 0.0]))
        @test !doesintersect(bbox, Ray([0.0, -2.0, 0.0], [0.0, -1.0, 0.0]))
        @test !doesintersect(bbox, Ray([0.0, 0.0, -2.0], [0.0, 0.0, -1.0]))

        bbox = BoundingBox(-0.5, 0.5, -0.75, 0.75, typemin(Float64), typemax(Float64))

        # starting outside
        r = Ray([-1.0, -1.1, 0.0], [1.0, 1.0, 0.0])
        intscts = surfaceintersection(bbox, r)
        @test isapprox(point(lower(intscts)), [-0.5, -0.6, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(intscts)), [0.5, 0.4, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE)

        # starting inside
        r = Ray([0.0, 0.0, 0.0], [1.0, 1.0, 0.0])
        intscts = surfaceintersection(bbox, r)
        @test lower(intscts) isa RayOrigin && isapprox(point(upper(intscts)), [0.5, 0.5, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE)

        # inside infinite
        r = Ray([0.0, 0.0, 0.0], [0.0, 0.0, 1.0])
        intscts = surfaceintersection(bbox, r)
        @test lower(intscts) isa RayOrigin && upper(intscts) isa Infinity

        # on surface
        r = Ray([0.5, 0.0, 0.0], [0.0, 0.0, 1.0])
        intscts = surfaceintersection(bbox, r)
        @test lower(intscts) isa RayOrigin && upper(intscts) isa Infinity

        # missing
        r = Ray([5.0, 5.0, 5.0], [0.0, 0.0, -1.0])
        intscts = surfaceintersection(bbox, r)
        @test intscts isa EmptyInterval

        r = Ray([5.0, 5.0, 5.0], [0.0, -1.0, 0.0])
        intscts = surfaceintersection(bbox, r)
        @test intscts isa EmptyInterval
    end # testset BoundingBox

    @testset "CSG" begin
        pln1 = Plane([0.0, 0.0, -1.0], [0.0, 0.0, -1.0])
        pln2 = Plane([0.0, 0.0, 1.0], [0.0, 0.0, 1.0])
        r = Ray([0.0, 0, 2.0], [0.0, 0.0, -1.0])
        gen1 = csgintersection(pln1, pln2)
        gen2 = csgunion(pln1, pln2)

        csg = gen1(identitytransform())
        intsct = evalcsg(csg, r)

        intvlpt1 = lower(intsct)
        intvlpt2 = upper(intsct)

        @test isapprox(point(intvlpt1), SVector{3}(0.0, 0.0, 1.0)) && isapprox(point(intvlpt2), SVector{3}(0.0, 0.0, -1.0))

        csg = gen2(identitytransform())
        intsct = evalcsg(csg, r)
        @test (lower(intsct) isa RayOrigin) && (upper(intsct) isa Infinity)

        # test simple rays for intersection, union and difference operations
        # INTERSECTION
        intersection_obj = csgintersection(leaf(Cylinder(0.5, 3.0), OpticSim.rotationd(90.0, 0.0, 0.0)), (leaf(Sphere(1.0))))()
        r = Ray([0.7, 0.0, 0.0], [-1.0, 0.0, 0.0])
        int = OpticSim.evalcsg(intersection_obj, r)
        @test isapprox(point(lower(int)), [0.5, 0.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(int)), [-0.5, 0.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE)

        r = Ray([5.0, 0.0, 0.0], [-1.0, 0.0, 0.0])
        int = OpticSim.evalcsg(intersection_obj, r)
        @test isapprox(point(lower(int)), [0.5, 0.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(int)), [-0.5, 0.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE)

        r = Ray([0.7, 0.0, 0.0], [1.0, 0.0, 0.0])
        int = OpticSim.evalcsg(intersection_obj, r)
        @test (int isa EmptyInterval)

        r = Ray([0.0, 0.0, 0.0], [1.0, 0.0, 0.0])
        int = OpticSim.evalcsg(intersection_obj, r)
        @test (lower(int) isa RayOrigin) && isapprox(point(upper(int)), [0.5, 0.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE)

        r = Ray([0.0, 0.0, 0.0], [0.0, 1.0, 0.0])
        int = OpticSim.evalcsg(intersection_obj, r)
        @test (lower(int) isa RayOrigin) && isapprox(point(upper(int)), [0.0, 1.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE)

        r = Ray([0.0, 2.0, 0.0], [0.0, -1.0, 0.0])
        int = OpticSim.evalcsg(intersection_obj, r)
        @test isapprox(point(lower(int)), [0.0, 1.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(int)), [0.0, -1.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE)

        # UNION
        union_obj = csgunion(leaf(Cylinder(0.5, 3.0)), (leaf(Sphere(1.0))))(OpticSim.translation(1.0, 0.0, 0.0))
        r = Ray([1.0, 0.0, 0.0], [0.0, 0.0, 1.0])
        int = OpticSim.evalcsg(union_obj, r)
        @test (lower(int) isa RayOrigin) && (upper(int) isa Infinity)

        r = Ray([1.0, 0.0, 0.0], [1.0, 0.0, 0.0])
        int = OpticSim.evalcsg(union_obj, r)
        @test (lower(int) isa RayOrigin) && isapprox(point(upper(int)), [2.0, 0.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE)

        r = Ray([1.6, 0.0, 0.0], [1.0, 0.0, 0.0])
        int = OpticSim.evalcsg(union_obj, r)
        @test (lower(int) isa RayOrigin) && isapprox(point(upper(int)), [2.0, 0.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE)

        r = Ray([1.6, 0.0, 0.0], [-1.0, 0.0, 0.0])
        int = OpticSim.evalcsg(union_obj, r)
        @test (lower(int) isa RayOrigin) && isapprox(point(upper(int)), [0.0, 0.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE)

        r = Ray([-1.0, 0.0, 0.0], [1.0, 0.0, 0.0])
        int = OpticSim.evalcsg(union_obj, r)
        @test isapprox(point(lower(int)), [0.0, 0.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(int)), [2.0, 0.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE)

        r = Ray([-1.0, 0.0, 4.0], [1.0, 0.0, 0.0])
        int = OpticSim.evalcsg(union_obj, r)
        @test isapprox(point(lower(int)), [0.5, 0.0, 4.0], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(int)), [1.5, 0.0, 4.0], rtol = RTOLERANCE, atol = ATOLERANCE)

        # DIFFERENCE
        difference_obj = csgdifference(leaf(Cylinder(0.5, 3.0)), leaf(Sphere(1.0), OpticSim.translation(0.75, 0.0, 0.2)))()
        r = Ray([0.25, 0.0, 0.0], [-1.0, 0.0, 0.0])
        int = OpticSim.evalcsg(difference_obj, r)
        @test isapprox(point(lower(int)), [-0.2297958971132712, 0.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(int)), [-0.5, 0.0, 0.0], rtol = RTOLERANCE, atol = ATOLERANCE)

        r = Ray([0.25, 0.0, 0.0], [1.0, 0.0, 0.0])
        int = OpticSim.evalcsg(difference_obj, r)
        @test int isa EmptyInterval

        r = Ray([0.25, 0.0, 0.0], [0.0, 0.0, 1.0])
        int = OpticSim.evalcsg(difference_obj, r)
        @test isapprox(point(lower(int)), [0.25, 0.0, 1.0660254037844386], rtol = RTOLERANCE, atol = ATOLERANCE) && (upper(int) isa Infinity)

        r = Ray([0.25, 0.0, 1.5], [0.0, 0.0, -1.0])
        int = OpticSim.evalcsg(difference_obj, r)
        @test (lower(int[1]) isa RayOrigin) && isapprox(point(upper(int[1])), [0.25, 0.0, 1.0660254037844386], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(lower(int[2])), [0.25, 0.0, -0.6660254037844386], rtol = RTOLERANCE, atol = ATOLERANCE) && (upper(int[2]) isa Infinity)

        r = Ray([0.25, 0.0, 1.5], [1.0, 0.0, 0.0])
        int = OpticSim.evalcsg(difference_obj, r)
        @test (lower(int) isa RayOrigin) && isapprox(point(upper(int)), [0.5, 0.0, 1.5], rtol = RTOLERANCE, atol = ATOLERANCE)

        r = Ray([0.25, 0.0, 1.5], [0.0, 0.0, 1.0])
        int = OpticSim.evalcsg(difference_obj, r)
        @test (lower(int) isa RayOrigin) && (upper(int) isa Infinity)

        # DisjointUnion result on CSG
        surf = TestData.verywavybeziersurface()
        accelsurf = leaf(AcceleratedParametricSurface(surf, 20))()

        # five hits starting outside
        r = Ray([0.9, 0.0, -0.3], [0.0, 1.0, 0.7])
        res = surfaceintersection(accelsurf, r)
        a = isapprox(point(lower(res[1])), [0.9, 0.03172286522032046, -0.2777939943457758], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(res[1])), [0.9, 0.1733979947040411, -0.17862140370717122], rtol = RTOLERANCE, atol = ATOLERANCE)
        b = isapprox(point(lower(res[2])), [0.9, 0.5335760974594397, 0.07350326822160776], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(res[2])), [0.9, 0.7767707607392784, 0.24373953251749486], rtol = RTOLERANCE, atol = ATOLERANCE)
        c = isapprox(point(lower(res[3])), [0.9, 0.9830891958374246, 0.3881624370861975], rtol = RTOLERANCE, atol = ATOLERANCE) && (upper(res[3]) isa Infinity)
        @test isa(res, DisjointUnion) && length(res) == 3 && a && b && c

        # five hits starting inside
        r = Ray([0.9, 1.0, 0.4], [0.0, -1.0, -0.7])
        res = surfaceintersection(accelsurf, r)
        a = (lower(res[1]) isa RayOrigin) && isapprox(point(upper(res[1])), [0.9, 0.9830891958374246, 0.3881624370861975], rtol = RTOLERANCE, atol = ATOLERANCE)
        b = isapprox(point(lower(res[2])), [0.9, 0.7767707607392784, 0.24373953251749486], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(res[2])), [0.9, 0.5335760974594397, 0.07350326822160776], rtol = RTOLERANCE, atol = ATOLERANCE)
        c = isapprox(point(lower(res[3])), [0.9, 0.17339799470404108, -0.1786214037071712], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(res[3])), [0.9, 0.03172286522032046, -0.27779399434577573], rtol = RTOLERANCE, atol = ATOLERANCE)
        @test isa(res, DisjointUnion) && length(res) == 3 && a && b && c
    end # testset CSG

    @testset "ThinGrating" begin
        angle_from_ray(raydirection) = 90 + atand(raydirection[3], raydirection[2])
        true_diff(order, λ, period, θi) = asind((order * λ / period + sind(θi)))
        period = 3.0
        int = TestData.transmissivethingrating(period, 2)
        for k in 1:50
            for θi in [-5.0, 0.0, 5.0, 10.0]
                for λ in [0.35, 0.55, 1.0]
                    ray = OpticalRay(SVector(0.0, 0.0, 2.0), SVector(0.0, sind(θi), -cosd(θi)), 1.0, λ)
                    raydir, _, _ = OpticSim.processintersection(int, SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 1.0), ray, 20.0, 1.0, true, true)
                    # order is random so just check that it is correct for one of them
                    @test any([isapprox(x, angle_from_ray(raydir), rtol = RTOLERANCE, atol = ATOLERANCE) for x in [true_diff(m, λ, period, θi) for m in -2:2]])
                end
            end
        end
    end # testset ThinGrating

    # TODO Hologram intersection tests

    @testset "Chebyshev" begin
        @test_throws AssertionError ChebyshevSurface(0.0, 1.0, [(1, 2, 1.0)])
        @test_throws AssertionError ChebyshevSurface(1.0, 0.0, [(1, 2, 1.0)])

        z1 = TestData.chebyshevsurface2()
        az1 = AcceleratedParametricSurface(z1, 20)
        numsamples = 100

        # test that an on axis ray doesn't miss
        @test !(surfaceintersection(AcceleratedParametricSurface(z1), Ray([0.0, 0.0, 10.0], [0.0, 0.0, -1.0])) isa EmptyInterval)

        samples = samplepoints(numsamples, -0.95, 0.95, -0.95, 0.95)
        missedintersections = 0

        for uv in samples
            pt = point(z1, uv[1], uv[2])
            origin = [0.0, 0.0, 5.0]
            r = Ray(origin, pt .- origin)
            allintersections = surfaceintersection(az1, r)
            if (allintersections isa EmptyInterval)
                # sampled triangle acceleration structure doesn't guarantee all intersecting rays will be detected.
                # Have to use patch subdivision and convex hull polyhedron to do this.
                missedintersections += 1
            else
                if isa(allintersections, Interval)
                    allintersections = (allintersections,)
                end

                for intsct in allintersections
                    if lower(intsct) isa RayOrigin
                        missedintersections += 1
                    else
                        # closest point should be on the surface, furthest should be on bounding prism
                        @test isapprox(point(lower(intsct)), pt, rtol = RTOLERANCE, atol = ATOLERANCE)
                        @test isa(upper(intsct), Intersection{Float64,3}) || (OpticSim.direction(r)[3] == -1 && isa(upper(intsct), Infinity{Float64}))
                    end
                end
            end
        end

        if missedintersections > 0
            @warn "$missedintersections out of total $(length(samples)) chebychev suface intersections were missed"
        end

        z1 = TestData.chebyshevsurface1()
        az1 = AcceleratedParametricSurface(z1, 20)

        # hit from inside
        r = Ray([0.0, 0.0, -0.5], [0.1, 0.2, 0.5])
        intvl = surfaceintersection(az1, r)
        @test (lower(intvl) isa RayOrigin) && isapprox(point(upper(intvl)), [0.07870079995991423, 0.15740159991982847, -0.10649600020042879], rtol = RTOLERANCE, atol = ATOLERANCE)

        # hit from outside
        r = Ray([0.0, 0.0, 0.5], [0.1, 0.2, -0.5])
        intvl = surfaceintersection(az1, r)
        @test isapprox(point(lower(intvl)), [0.13584702907541094, 0.2716940581508219, -0.1792351453770547], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(intvl)), [1.0, 2.0, -4.5], rtol = RTOLERANCE, atol = ATOLERANCE)

        # miss from inside
        r = Ray([0.2, 0.2, -0.5], [0.0, 0.0, -1.0])
        intvl = surfaceintersection(az1, r)
        @test (lower(intvl) isa RayOrigin) && (upper(intvl) isa Infinity)

        # miss from outside
        r = Ray([0.2, 3.0, -0.5], [0.0, 0.0, 1.0])
        @test surfaceintersection(az1, r) isa EmptyInterval
        r = Ray([0.2, 3.0, 0.5], [0.0, 0.0, -1.0])
        @test surfaceintersection(az1, r) isa EmptyInterval
        r = Ray([0.2, 0.2, 0.5], [0.0, 0.0, 1.0])
        @test surfaceintersection(az1, r) isa EmptyInterval
        r = Ray([0.2, 3.0, 2.0], [0.0, -1.0, 0.0])
        @test surfaceintersection(az1, r) isa EmptyInterval

        # two hits from outside
        r = Ray([0.0, -1.5, 0.25], [-1.0, 0.0, 0.0])
        intvl = surfaceintersection(az1, r)
        @test isapprox(point(lower(intvl)), [-1.030882920068228, -1.5, 0.25], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(intvl)), [-1.786071230727613, -1.5, 0.25], rtol = RTOLERANCE, atol = ATOLERANCE)

        r = Ray([1.5, -0.8, 0.25], [-1.0, 0.0, 0.0])
        hit = surfaceintersection(az1, r)
        a = isapprox(point(lower(hit[1])), [0.21526832427065165, -0.8, 0.25], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(hit[1])), [-0.25523748548950076, -0.8, 0.25], rtol = RTOLERANCE, atol = ATOLERANCE)
        b = isapprox(point(lower(hit[2])), [-1.9273057166986296, -0.8, 0.25], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(hit[2])), [-2.0, -0.8, 0.25], rtol = RTOLERANCE, atol = ATOLERANCE)
        @test isa(hit, DisjointUnion) && length(hit) == 2 && a && b

        # hit cyl from outside
        r = Ray([0.0, 3.0, -0.5], [0.0, -1.0, 0.0])
        intvl = surfaceintersection(az1, r)
        @test isapprox(point(lower(intvl)), [0.0, 2.0, -0.5], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(intvl)), [0.0, -2.0, -0.5], rtol = RTOLERANCE, atol = ATOLERANCE)

        # hit cyl from inside
        r = Ray([0.0, 0.0, -0.5], [0.0, -1.0, 0.0])
        intvl = surfaceintersection(az1, r)
        @test (lower(intvl) isa RayOrigin) && isapprox(point(upper(intvl)), [0.0, -2.0, -0.5], rtol = RTOLERANCE, atol = ATOLERANCE)

    end # testset Chebyshev

    @testset "GridSag" begin
        z1 = TestData.gridsagsurfacebicubic()
        az1 = AcceleratedParametricSurface(z1, 20)

        numsamples = 100

        # test that an on axis ray doesn't miss
        @test !(surfaceintersection(AcceleratedParametricSurface(z1), Ray([0.0, 0.0, 10.0], [0.0, 0.0, -1.0])) isa EmptyInterval)

        samples1 = samplepoints(numsamples, 0.01, 0.99, 0.01, 2π - 0.01)
        missedintersections = 0

        for uv in samples1
            pt = point(z1, uv[1], uv[2])
            origin = [0.0, 0.0, 5.0]
            r = Ray(origin, pt .- origin)
            allintersections = surfaceintersection(az1, r)
            if (allintersections isa EmptyInterval)
                # sampled triangle acceleration structure doesn't guarantee all intersecting rays will be detected.
                # Have to use patch subdivision and convex hull polyhedron to do this.
                missedintersections += 1
            else
                if isa(allintersections, Interval)
                    allintersections = (allintersections,)
                end

                for intsct in allintersections
                    if lower(intsct) isa RayOrigin
                        missedintersections += 1
                    else
                        # closest point should be on the surface, furthest should be on bounding prism
                        @test isapprox(point(lower(intsct)), pt, rtol = RTOLERANCE, atol = ATOLERANCE)
                        @test isa(upper(intsct), Intersection{Float64,3}) || (OpticSim.direction(r)[3] == -1 && isa(upper(intsct), Infinity{Float64}))
                    end
                end
            end
        end

        z2 = TestData.gridsagsurfacebicubiccheby()
        az2 = AcceleratedParametricSurface(z2, 20)

        samples2 = samplepoints(numsamples, -0.95, 0.95, -0.95, 0.95)
        for uv in samples2
            pt = point(z2, uv[1], uv[2])
            origin = [0.0, 0.0, 5.0]
            r = Ray(origin, pt .- origin)
            allintersections = surfaceintersection(az2, r)
            if (allintersections isa EmptyInterval)
                # sampled triangle acceleration structure doesn't guarantee all intersecting rays will be detected.
                # Have to use patch subdivision and convex hull polyhedron to do this.
                missedintersections += 1
            else
                if isa(allintersections, Interval)
                    allintersections = (allintersections,)
                end

                for intsct in allintersections
                    if lower(intsct) isa RayOrigin
                        missedintersections += 1
                    else
                        # closest point should be on the surface, furthest should be on bounding prism
                        @test isapprox(point(lower(intsct)), pt, rtol = RTOLERANCE, atol = ATOLERANCE)
                        @test isa(upper(intsct), Intersection{Float64,3}) || (OpticSim.direction(r)[3] == -1 && isa(upper(intsct), Infinity{Float64}))
                    end
                end
            end
        end

        if missedintersections > 0
            @warn "$missedintersections out of total $(length(samples1) + length(samples2)) gridsag suface intersections were missed"
        end

        # hit from inside
        r = Ray([0.0, 0.0, -0.5], [0.1, 0.2, 0.5])
        intvl = surfaceintersection(az1, r)
        @test (lower(intvl) isa RayOrigin) && isapprox(point(upper(intvl)), [0.13038780645120746, 0.26077561290241486, 0.15193903225603717], rtol = RTOLERANCE, atol = ATOLERANCE)

        # hit from outside
        r = Ray([0.0, 0.0, 0.5], [0.1, 0.2, -0.5])
        intvl = surfaceintersection(az1, r)
        @test isapprox(point(lower(intvl)), [0.046677740288643785, 0.0933554805772876, 0.26661129855678095], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(intvl)), [0.6708203932499369, 1.3416407864998738, -2.8541019662496843], rtol = RTOLERANCE, atol = ATOLERANCE)

        # miss from inside
        r = Ray([0.2, 0.2, -0.5], [0.0, 0.0, -1.0])
        intvl = surfaceintersection(az1, r)
        @test (lower(intvl) isa RayOrigin) && (upper(intvl) isa Infinity)

        # miss from outside
        r = Ray([0.2, 2.0, -0.5], [0.0, 0.0, 1.0])
        @test surfaceintersection(az1, r) isa EmptyInterval
        r = Ray([0.2, 2.0, 0.5], [0.0, 0.0, -1.0])
        @test surfaceintersection(az1, r) isa EmptyInterval
        r = Ray([0.2, 0.2, 0.5], [0.0, 0.0, 1.0])
        @test surfaceintersection(az1, r) isa EmptyInterval
        r = Ray([0.2, 2.0, 2.0], [0.0, -1.0, 0.0])
        @test surfaceintersection(az1, r) isa EmptyInterval

        # two hits from outside
        r = Ray([2.0, 0.0, 0.25], [-1.0, 0.0, 0.0])
        hit = surfaceintersection(az1, r)
        a = isapprox(point(lower(hit[1])), [1.5, 0.0, 0.25], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(hit[1])), [1.394132887629247, 0.0, 0.25], rtol = RTOLERANCE, atol = ATOLERANCE)
        b = isapprox(point(lower(hit[2])), [0.30184444700168467, 0.0, 0.25], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(hit[2])), [-0.6437764440680096, 0.0, 0.25], rtol = RTOLERANCE, atol = ATOLERANCE)
        @test isa(hit, DisjointUnion) && length(hit) == 2 && a && b

        # hit cyl from outside
        r = Ray([0.0, 2.0, -0.5], [0.0, -1.0, 0.0])
        intvl = surfaceintersection(az1, r)
        @test isapprox(point(lower(intvl)), [0.0, 1.5, -0.5], rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(point(upper(intvl)), [0.0, -1.5, -0.5], rtol = RTOLERANCE, atol = ATOLERANCE)

        # hit cyl from inside
        r = Ray([0.0, 0.0, -0.5], [0.0, -1.0, 0.0])
        intvl = surfaceintersection(az1, r)
        @test (lower(intvl) isa RayOrigin) && isapprox(point(upper(intvl)), [0.0, -1.5, -0.5], rtol = RTOLERANCE, atol = ATOLERANCE)
    end # testset GridSag
end # testset intersection

@otestset "Lenses" begin
    test_n = 5000

    @testset "Refraction" begin
        Random.seed!(SEED)
        ηᵢ = 1.4
        ηₜ = 1.2
        for i in 1:test_n
            nₛ = randunit()
            rᵢ = randunit()
            rₜ = refractedray(ηᵢ, ηₜ, nₛ, rᵢ)
            if !(rₜ === nothing)
                sinθₜ = norm(cross(-nₛ, rₜ))
                sinθᵢ = norm(cross(nₛ, rᵢ))
                a = isapprox(ηₜ * sinθₜ, ηᵢ * sinθᵢ, rtol = RTOLERANCE, atol = ATOLERANCE) #verify snell's law for the the transmitted and incident ray
                b = isapprox(1.0, norm(rₜ), rtol = RTOLERANCE, atol = ATOLERANCE) #verify transmitted ray is unit
                perp = normalize(cross(nₛ, rᵢ))
                c = isapprox(0.0, rₜ ⋅ perp, rtol = RTOLERANCE, atol = ATOLERANCE) #verify incident and transmitted ray are in the same plane
                @test a && b && c
            end
        end
    end # testset refraction

    # snell
    @testset "Snell" begin
        Random.seed!(SEED)
        for i in 1:test_n
            r = randunit()
            nₛ = randunit()
            n1 = 1.0
            n2 = 1.4
            sθ1, sθ2 = snell(nₛ, r, n1, n2)
            sθ3, sθ4 = snell(nₛ, r, n2, n1)
            @test isapprox(sθ1 * n1, sθ2 * n2, rtol = RTOLERANCE, atol = ATOLERANCE) && isapprox(sθ3 * n2, sθ4 * n1, rtol = RTOLERANCE, atol = ATOLERANCE)
        end
    end # testset snell

    @testset "Reflection" begin
        Random.seed!(SEED)
        # reflection
        for i in 1:test_n
            r = randunit()
            nₛ = randunit()
            if abs(r ⋅ nₛ) < 0.9999 # if r and n are exactly aligned then their sum will be zero, not a vector pointing in the direction of the normal
                reflected = reflectedray(nₛ, r)
                #ensure reflected and refracted rays are in the same plane
                a = isapprox(norm(reflected ⋅ cross(r, nₛ)), 0.0, rtol = RTOLERANCE, atol = ATOLERANCE)
                b = isapprox(norm(reflected), 1.0, rtol = RTOLERANCE, atol = ATOLERANCE)
                c = isapprox(0.0, reflected ⋅ nₛ + r ⋅ nₛ, rtol = RTOLERANCE, atol = ATOLERANCE)
                #ensure sum of reflected and origin ray are in the direction of the normal
                d = isapprox(0.0, norm(nₛ - sign(reflected ⋅ nₛ) * (normalize(reflected - r))), rtol = RTOLERANCE, atol = ATOLERANCE)
                @test a & b & c & d
            end
        end
    end # testset reflection

    @testset "Paraxial" begin
        # check that normally incident rays are focussed across the lens
        l = ParaxialLensEllipse(100.0, 10.0, 10.0, [1.0, 1.0, 0.0], [3.0, 3.0, 3.0])
        r = OpticalRay([2.0, 2.0, 3.0], [1.0, 1.0, 0.0], 1.0, 0.55)
        intsct = halfspaceintersection(surfaceintersection(l, r))
        ref, _, _ = OpticSim.processintersection(OpticSim.interface(intsct), OpticSim.point(intsct), OpticSim.normal(intsct), r, 20.0, 1.0, true, true)
        @test isapprox(ref, normalize([1.0, 1.0, 0.0]), rtol = RTOLERANCE, atol = ATOLERANCE)
        l = ParaxialLensRect(100.0, 10.0, 10.0, [1.0, 1.0, 0.0], [3.0, 3.0, 3.0])
        r = OpticalRay([2.3, 2.4, 3.1], [1.0, 1.0, 0.0], 1.0, 0.55)
        intsct = halfspaceintersection(surfaceintersection(l, r))
        fp = [3.0, 3.0, 3.0] + 100 * normalize([1.0, 1.0, 0.0])
        ref, _, _ = OpticSim.processintersection(OpticSim.interface(intsct), OpticSim.point(intsct), OpticSim.normal(intsct), r, 20.0, 1.0, true, true)
        @test isapprox(ref, normalize(fp - OpticSim.point(intsct)), rtol = RTOLERANCE, atol = ATOLERANCE)
    end # testset Paraxial
end # testset lenses

@otestset "Emitters" begin
    # TODO!! once emitters are finalised
end # testset Emitters

@otestset "Comparison" begin
    @testset "Refraction" begin
        λ = 0.550
        r1 = OpticalRay([0.0, 0.0, 1.0], [0.0, 0.0, -1.0], 1.0, λ)
        r2 = OpticalRay([2.0, 2.0, 1.0], [0.0, 0.0, -1.0], 1.0, λ)
        r3 = OpticalRay([5.0, 5.0, 1.0], [0.0, 0.0, -1.0], 1.0, λ)
        r4 = OpticalRay([0.0, -5.0, 1.0], [0.0, 0.08715574274765818, -0.9961946980917454], 1.0, λ)
        r5 = OpticalRay([-5.0, -5.0, 1.0], [0.08715574274765818, -0.01738599476176408, -0.9960429728140486], 1.0, λ)

        # Test a number of rays under standard environmental conditions
        a = TestData.planoplano()
        res = trace(a, r1, test = true)
        @test isapprox(point(res), [0.0, 0.0, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 73.9852238762079, rtol = COMP_TOLERANCE)
        res = trace(a, r2, test = true)
        @test isapprox(point(res), [2.0, 2.0, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 73.9852238762079, rtol = COMP_TOLERANCE)
        res = trace(a, r3, test = true)
        @test isapprox(point(res), [5.0, 5.0, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 73.9852238762079, rtol = COMP_TOLERANCE)
        res = trace(a, r4, test = true)
        @test isapprox(point(res), [0.0, 0.7192321011619232, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 74.234903851062, rtol = COMP_TOLERANCE)
        res = trace(a, r5, test = true)
        @test isapprox(point(res), [0.7200535362928893, -6.141047252697188, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 74.2448952769065, rtol = COMP_TOLERANCE)

        a = TestData.concaveplano()
        res = trace(a, r1, test = true)
        @test isapprox(point(res), [0.0, 0.0, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 73.9852238762079, rtol = COMP_TOLERANCE)
        res = trace(a, r2, test = true)
        @test isapprox(point(res), [3.633723967570758, 3.633723967570758, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 74.0772722321026, rtol = COMP_TOLERANCE)
        res = trace(a, r3, test = true)
        @test isapprox(point(res), [9.153654757938119, 9.153654757938119, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 74.5695599263722, rtol = COMP_TOLERANCE)
        res = trace(a, r4, test = true)
        @test isapprox(point(res), [0.0, -3.401152468994745, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 74.1609946836496, rtol = COMP_TOLERANCE)
        res = trace(a, r5, test = true)
        @test isapprox(point(res), [-3.457199906282556, -10.32836827735406, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 74.5393056652617, rtol = COMP_TOLERANCE)

        a = TestData.doubleconcave()
        res = trace(a, r1, test = true)
        @test isapprox(point(res), [0.0, 0.0, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 73.9852238762079, rtol = COMP_TOLERANCE)
        res = trace(a, r2, test = true)
        @test isapprox(point(res), [5.235579681684886, 5.235579681684886, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 74.2624909340242, rtol = COMP_TOLERANCE)
        res = trace(a, r3, test = true)
        @test isapprox(point(res), [13.42351905137007, 13.42351905137007, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 75.8131310483928, rtol = COMP_TOLERANCE)
        res = trace(a, r4, test = true)
        @test isapprox(point(res), [0.0, -6.904467939352202, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 74.3286918571586, rtol = COMP_TOLERANCE)
        res = trace(a, r5, test = true)
        @test isapprox(point(res), [-7.089764102262856, -14.61837033417989, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 75.4286280772859, rtol = COMP_TOLERANCE)

        a = TestData.convexplano()
        res = trace(a, r1, test = true)
        @test isapprox(point(res), [0.0, 0.0, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 73.9852238762079, rtol = COMP_TOLERANCE)
        res = trace(a, r2, test = true)
        @test isapprox(point(res), [0.8867891519368289, 0.8867891519368289, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 73.9698941263224, rtol = COMP_TOLERANCE)
        res = trace(a, r3, test = true)
        @test isapprox(point(res), [2.212067233969831, 2.212067233969831, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 73.8895474037682, rtol = COMP_TOLERANCE)
        res = trace(a, r4, test = true)
        @test isapprox(point(res), [0.0, 3.489759589087602, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 74.4310036252394, rtol = COMP_TOLERANCE)
        res = trace(a, r5, test = true)
        @test isapprox(point(res), [3.491858246934857, -3.331802834115089, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 74.344370622227, rtol = COMP_TOLERANCE)

        a = TestData.doubleconvex()
        res = trace(a, r1, test = true)
        @test isapprox(point(res), [0.0, 0.0, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 73.9852238762079, rtol = COMP_TOLERANCE)
        res = trace(a, r2, test = true)
        @test isapprox(point(res), [-0.06191521590711035, -0.06191521590711035, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 73.987421926598, rtol = COMP_TOLERANCE)
        res = trace(a, r3, test = true)
        @test isapprox(point(res), [-0.2491105067897657, -0.2491105067897657, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 74.0110871422915, rtol = COMP_TOLERANCE)
        res = trace(a, r4, test = true)
        @test isapprox(point(res), [0.0, 5.639876913179362, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 74.6758349889372, rtol = COMP_TOLERANCE)
        res = trace(a, r5, test = true)
        @test isapprox(point(res), [5.705452331562673, -0.7679713587894854, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 74.6175821043825, rtol = COMP_TOLERANCE)
    end # testset Refraction

    @testset "Temperature/Pressure" begin
        λ = 0.550
        r1 = OpticalRay([0.0, 0.0, 1.0], [0.0, 0.0, -1.0], 1.0, λ)
        r2 = OpticalRay([2.0, 2.0, 1.0], [0.0, 0.0, -1.0], 1.0, λ)
        r3 = OpticalRay([5.0, 5.0, 1.0], [0.0, 0.0, -1.0], 1.0, λ)
        r4 = OpticalRay([0.0, -5.0, 1.0], [0.0, 0.08715574274765818, -0.9961946980917454], 1.0, λ)
        r5 = OpticalRay([-5.0, -5.0, 1.0], [0.08715574274765818, -0.01738599476176408, -0.9960429728140486], 1.0, λ)

        # Test a component whose material does have a ΔT component
        a = TestData.doubleconvex(temperature = 40 * u"°C")
        @test isapprox(point(intersection(trace(a, r1, test = true))), [0.0, 0.0, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(point(intersection(trace(a, r2, test = true))), [-0.06214282132053373, -0.06214282132053373, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(point(intersection(trace(a, r3, test = true))), [-0.2497005807710933, -0.2497005807710933, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(point(intersection(trace(a, r4, test = true))), [0.0, 5.640416741927821, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(point(intersection(trace(a, r5, test = true))), [5.706006107804734, -0.7673622008537766, -67.8], rtol = COMP_TOLERANCE)

        # note that this material doesn't have a ΔT component
        λ2 = 0.533
        r1 = OpticalRay([0.0, 0.0, 1.0], [0.0, 0.0, -1.0], 1.0, λ2)
        r2 = OpticalRay([2.0, 2.0, 1.0], [0.0, 0.0, -1.0], 1.0, λ2)
        r3 = OpticalRay([5.0, 5.0, 1.0], [0.0, 0.0, -1.0], 1.0, λ2)
        r4 = OpticalRay([0.0, -5.0, 1.0], [0.0, 0.08715574274765818, -0.9961946980917454], 1.0, λ2)
        r5 = OpticalRay([-5.0, -5.0, 1.0], [0.08715574274765818, -0.01738599476176408, -0.9960429728140486], 1.0, λ2)
    end # testset Temperature/Pressure

    @testset "Reflection" begin
        λ = 0.550
        r1 = OpticalRay([0.0, 0.0, 1.0], [0.0, 0.0, -1.0], 1.0, λ)
        r2 = OpticalRay([2.0, 2.0, 1.0], [0.0, 0.0, -1.0], 1.0, λ)
        r3 = OpticalRay([5.0, 5.0, 1.0], [0.0, 0.0, -1.0], 1.0, λ)
        r4 = OpticalRay([0.0, -5.0, 1.0], [0.0, 0.08715574274765818, -0.9961946980917454], 1.0, λ)
        r5 = OpticalRay([-5.0, -5.0, 1.0], [0.08715574274765818, -0.01738599476176408, -0.9960429728140486], 1.0, λ)

        a = TestData.planoconcaverefl()
        res = trace(a, r1, test = true)
        @test isapprox(point(res), [0.0, 0.0, 47.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 79.1704477524159, rtol = COMP_TOLERANCE)
        res = trace(a, r2, test = true)
        @test isapprox(point(res), [-6.19723306038804, -6.19723306038804, 47.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 80.0979229417755, rtol = COMP_TOLERANCE)
        res = trace(a, r3, test = true)
        @test isapprox(point(res), [-17.7546255175372, -17.7546255175372, 47.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 86.2063094599075, rtol = COMP_TOLERANCE)
        res = trace(a, r4, test = true)
        @test isapprox(point(res), [0.0, 19.5349836031196, 47.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 83.5364646376395, rtol = COMP_TOLERANCE)
        res = trace(a, r5, test = true)
        @test isapprox(point(res), [21.2477706541112, 17.3463704045055, 47.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 87.4033932087711, rtol = COMP_TOLERANCE)
    end # testset Reflection

    @testset "Complex Lenses" begin
        λ = 0.550
        r1 = OpticalRay([0.0, 0.0, 1.0], [0.0, 0.0, -1.0], 1.0, λ)
        r2 = OpticalRay([2.0, 2.0, 1.0], [0.0, 0.0, -1.0], 1.0, λ)
        r3 = OpticalRay([5.0, 5.0, 1.0], [0.0, 0.0, -1.0], 1.0, λ)
        r4 = OpticalRay([0.0, -5.0, 1.0], [0.0, 0.08715574274765818, -0.9961946980917454], 1.0, λ)
        r5 = OpticalRay([-5.0, -5.0, 1.0], [0.08715574274765818, -0.01738599476176408, -0.9960429728140486], 1.0, λ)

        a = TestData.conicsystemZ()
        res = trace(a, r1, test = true)
        @test isapprox(point(res), [0.0, 0.0, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 73.9852238762079, rtol = COMP_TOLERANCE)
        res = trace(a, r2, test = true)
        @test isapprox(point(res), [-1.80279270185495, -1.80279270185495, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 74.1000975813922, rtol = COMP_TOLERANCE)
        res = trace(a, r3, test = true)
        @test isapprox(point(res), [-2.8229241607807, -2.8229241607807, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 74.282094499601, rtol = COMP_TOLERANCE)
        res = trace(a, r4, test = true)
        @test isapprox(point(res), [0.0, 9.01421545841289, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 75.2211078071893, rtol = COMP_TOLERANCE)
        res = trace(a, r5, test = true)
        @test isapprox(point(res), [8.13946939328266, 2.23006981338816, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 75.1025133542546, rtol = COMP_TOLERANCE)

        a = TestData.asphericsystem()
        res = trace(a, r1, test = true)
        @test isapprox(point(res), [0.0, 0.0, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 73.9852238762079, rtol = COMP_TOLERANCE)
        res = trace(a, r2, test = true)
        @test isapprox(point(res), [0.0735282671574837, 0.0735282671574837, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 74.0161411455851, rtol = COMP_TOLERANCE)
        track = Vector{OpticSim.LensTrace{Float64,3}}(undef, 0)
        res = trace(a, r3, test = true, trackrays = track)
        @test (res === nothing) # TIR
        @test isapprox(point(track[end]), [-5.0, -5.0, 0.0], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(track[end]), 46.5556716286238, rtol = COMP_TOLERANCE)
        res = trace(a, r4, test = true)
        @test isapprox(point(res), [0.0, 11.9748998399885, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 76.0760286320348, rtol = COMP_TOLERANCE)
        track = Vector{OpticSim.LensTrace{Float64,3}}(undef, 0)
        res = trace(a, r5, test = true, trackrays = track)
        @test (res === nothing) # TIR
        @test isapprox(point(track[end]), [5.49905367197174, 5.66882664623822, 0.0], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(track[end]), 47.6825931025333, rtol = COMP_TOLERANCE)

        a = TestData.zernikesystem()
        res = trace(a, r1, test = true)
        @test isapprox(point(res), [0.0, -8.9787010034042, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 75.5977029819277, rtol = COMP_TOLERANCE)
        res = trace(a, r2, test = true)
        @test isapprox(point(res), [-1.58696749235066, -9.71313213852721, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 75.7603266537023, rtol = COMP_TOLERANCE)
        res = trace(a, r3, test = true)
        @test isapprox(point(res), [0.790348081563859, -0.762155123619682, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 74.1933359669848, rtol = COMP_TOLERANCE)
        res = trace(a, r4, test = true, trackrays = track)
        @test isapprox(point(res), [0.0, 10.0037899692172, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 76.8041236301678, rtol = COMP_TOLERANCE)
        res = trace(a, r5, test = true)
        @test isapprox(point(res), [35.5152289731456, 28.3819941055557, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 94.3059947823954, rtol = COMP_TOLERANCE)

        a = TestData.conicsystemQ()
        res = trace(a, r1, test = true)
        @test isapprox(point(res), [0.0, 0.0, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 73.9852238762079, rtol = COMP_TOLERANCE)
        res = trace(a, r2, test = true)
        @test isapprox(point(res), [-1.80279270185495, -1.80279270185495, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 74.1000975813922, rtol = COMP_TOLERANCE)
        res = trace(a, r3, test = true)
        @test isapprox(point(res), [-2.8229241607807, -2.8229241607807, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 74.282094499601, rtol = COMP_TOLERANCE)
        res = trace(a, r4, test = true)
        @test isapprox(point(res), [0.0, 9.01421545841289, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 75.2211078071893, rtol = COMP_TOLERANCE)
        res = trace(a, r5, test = true)
        @test isapprox(point(res), [8.13946939328266, 2.23006981338816, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 75.1025133542546, rtol = COMP_TOLERANCE)

        # No NSC qtype so have less precise SC values and no angled rays
        a = TestData.qtypesystem()
        res = trace(a, r1, test = true)
        @test isapprox(point(res), [-4.3551074306, -0.318112811010, -67.8], rtol = 1e-12)
        # TODO these are the values that I get for comparison - I'm pretty certain that they are just wrong...
        # res = trace(a, r2, test = true)
        # @test isapprox(point(res), [-0.091773202667, 3.9362695851, -67.8], rtol = 1e-12)
        # res = trace(a, r3, test = true)
        # @test isapprox(point(res), [-7.0993098547, -2.6848726242, -67.8], rtol = 1e-12)

        # No NSC chebyshev so have less precise SC values and no angled rays
        a = TestData.chebyshevsystem()
        res = trace(a, r1, test = true)
        @test isapprox(point(res), [-1.07963031980, 0.53981515992, -67.8], rtol = 1e-12)
        res = trace(a, r2, test = true)
        @test isapprox(point(res), [0.00851939011, 1.25870229260, -67.8], rtol = 1e-12)
        res = trace(a, r3, test = true)
        @test isapprox(point(res), [-1.20441411240, -0.07554605390, -67.8], rtol = 1e-12)

        a = TestData.gridsagsystem()
        res = trace(a, r1, test = true)
        @test isapprox(point(res), [21.0407756733608, 21.724638830759, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 82.0777022031378, rtol = COMP_TOLERANCE)
        res = trace(a, r2, test = true)
        @test isapprox(point(res), [-0.489183765274452, 0.405160352533666, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 75.536902213118, rtol = COMP_TOLERANCE)
        res = trace(a, r3, test = true)
        @test isapprox(point(res), [-12.0770793886528, -8.7705340321259, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 77.5651917266687, rtol = COMP_TOLERANCE)
        res = trace(a, r4, test = true)
        @test isapprox(point(res), [-1.20535362019229, 6.07973526939659, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 74.7349148204077, rtol = COMP_TOLERANCE)
        res = trace(a, r5, test = true)
        @test isapprox(point(res), [6.5112407340601, -0.22055440245024, -67.8], rtol = COMP_TOLERANCE)
        @test isapprox(pathlength(res), 74.6955888563913, rtol = COMP_TOLERANCE)
    end #testset complex lenses

    # @testset "Power" begin
    #     λ = 0.550
    #     r1 = OpticalRay([0.0, 0.0, 1.0], [0.0, 0.0, -1.0], 1.0, λ)
    #     r2 = OpticalRay([2.0, 2.0, 1.0], [0.0, 0.0, -1.0], 1.0, λ)
    #     r3 = OpticalRay([5.0, 5.0, 1.0], [0.0, 0.0, -1.0], 1.0, λ)
    #     a = TestData.doubleconvex()
    #     res = trace(a, r1, test = true)
    #     @test isapprox(power(res), 0.915508396)
    #     res = trace(a, r2, test = true)
    #     @test isapprox(power(res), 0.915526077)
    #     res = trace(a, r3, test = true)
    #     @test isapprox(power(res), 0.915577586)
    # end
end # testset Comparison

@otestset "Visualization" begin
    # test that this all at least runs
    surf1 = AcceleratedParametricSurface(TestData.beziersurface(), 15)
    surf2 = AcceleratedParametricSurface(TestData.upsidedownbeziersurface(), 15)
    m = csgintersection(leaf(surf1, translation(-0.5, -0.5, 0.0)), csgintersection(leaf(Cylinder(0.3, 5.0)), leaf(surf1, RigidBodyTransform{Float64}(0.0, Float64(π), 0.0, 0.5, -0.5, 0.0))))()
    @test_nowarn Vis.draw(m)
    m = csgintersection(leaf(surf1, translation(-0.5, -0.5, 0.0)), csgintersection(leaf(Cylinder(0.3, 5.0)), leaf(surf2, translation(-0.5, -0.5, 0.0))))()
    @test_nowarn Vis.draw(m)
    m = csgunion(leaf(Cylinder(0.5, 3.0)), (leaf(Sphere(1.0))))()
    @test_nowarn Vis.draw(m)
    m = csgintersection(leaf(Cylinder(0.5, 3.0)), (leaf(Sphere(1.0))))()
    @test_nowarn Vis.draw(m)
    m = csgdifference(leaf(Cylinder(0.5, 3.0)), leaf(Sphere(1.0)))()
    @test_nowarn Vis.draw(m)
end # testset Visualization

@otestset "Allocations" begin
    # ensure that there are 0 allocations for all the benchmarks
    for b in Benchmarks.all_benchmarks()
        @test Benchmarks.runbenchmark(b, samples = 1, evals = 1).allocs == 0
    end
end # testset Allocations
