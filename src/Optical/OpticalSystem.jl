# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

"""
    AbstractOpticalSystem{T<:Real}

Abstract type for any optical system, must parameterized by the datatype of entities within the system `T`.
"""
abstract type AbstractOpticalSystem{T<:Real} end
export AbstractOpticalSystem

"""
    CSGOpticalSystem{T,D<:Real,S<:Surface{T},L<:LensAssembly{T}} <: AbstractOpticalSystem{T}

An optical system containing a lens assembly with all optical elements and a detector surface with associated image. The system can be at a specified temperature and pressure.

There are two number types in the type signature. The `T` type parameter is the numeric type for geometry in the optical system, the `D` type parameter is the numeric type of the pixels in the detector image.
This way you can have `Float64` geometry, where high precision is essential, but the pixels in the detector can be `Float32` since precision is much less critical for image data.

The detector can be any [`Surface`](@ref) which implements [`uv`](@ref), [`uvtopix`](@ref) and [`onsurface`](@ref), typically this is one of [`Rectangle`](@ref), [`Ellipse`](@ref) or [`SphericalCap`](@ref).

```julia
CSGOpticalSystem(assembly::LensAssembly, detector::Surface, detectorpixelsx = 1000, detectorpixelsy = 1000, ::Type{D} = Float32; temperature = OpticSim.GlassCat.TEMP_REF, pressure = OpticSim.GlassCat.PRESSURE_REF)
```
"""
struct CSGOpticalSystem{T,D<:Number,S<:Surface{T},L<:LensAssembly{T}} <: AbstractOpticalSystem{T}
    assembly::L
    detector::S
    detectorimage::HierarchicalImage{D}
    temperature::T
    pressure::T
    function CSGOpticalSystem(assembly::L, detector::S, detectorpixelsx::Int = 1000, detectorpixelsy::Int = 1000, ::Type{D} = Float32; temperature::Union{T,Unitful.Temperature} = T(OpticSim.GlassCat.TEMP_REF), pressure::T = T(OpticSim.GlassCat.PRESSURE_REF)) where {T<:Real,S<:Surface{T},L<:LensAssembly{T},D<:Number}
        @assert hasmethod(uv, (S, SVector{3,T})) "Detector must implement uv()"
        @assert hasmethod(uvtopix, (S, SVector{2,T}, Tuple{Int,Int})) "Detector must implement uvtopix()"
        @assert hasmethod(onsurface, (S, SVector{3,T})) "Detector must implement onsurface()"
        opticalinterface = interface(detector)
        @assert insidematerialid(opticalinterface) == outsidematerialid(opticalinterface) "Detector must have same material either side"
        @assert interface(detector) !== NullInterface(T) "Detector can't have null interface"
        image = HierarchicalImage{D}(detectorpixelsy, detectorpixelsx)
        if temperature isa Unitful.Temperature
            temperature = Unitful.ustrip(T, Unitful.u"°C", temperature)
        end
        return new{T,D,S,L}(assembly, detector, image, temperature, T(pressure))
    end
end
export CSGOpticalSystem

Base.copy(a::CSGOpticalSystem) = CSGOpticalSystem(a.assembly, a.detector, size(a.detectorimage)..., temperature = a.temperature * Unitful.u"°C", pressure = a.pressure)

#added this show method because the type of CSGOpticalSystem is gigantic and printing it in the REPL can crash the system.
Base.show(io::IO, a::CSGOpticalSystem{T}) where {T} = print(io, "CSGOpticalSystem{$T}($(temperature(a)), $(pressure(a)), $(assembly(a)), $(detector(a)))")

"""
    assembly(system::AbstractOpticalSystem{T}) -> LensAssembly{T}

Get the [`LensAssembly`](@ref) of `system`.
"""
assembly(system::CSGOpticalSystem{T}) where {T<:Real} = system.assembly
detector(system::CSGOpticalSystem) = system.detector
"""
    detectorimage(system::AbstractOpticalSystem{T}) -> HierarchicalImage{D}

Get the detector image of `system`.
`D` is the datatype of the detector image and is not necessarily the same as the datatype of the system `T`.
"""
detectorimage(system::CSGOpticalSystem) = system.detectorimage
detectorsize(system::CSGOpticalSystem) = size(system.detectorimage)
"""
    temperature(system::AbstractOpticalSystem{T}) -> T

Get the temperature of `system` in °C.
"""
temperature(system::CSGOpticalSystem{T}) where {T<:Real} = system.temperature
"""
    pressure(system::AbstractOpticalSystem{T}) -> T

Get the pressure of `system` in Atm.
"""
pressure(system::CSGOpticalSystem{T}) where {T<:Real} = system.pressure
"""
    resetdetector!(system::AbstractOpticalSystem{T})

Reset the deterctor image of `system` to zero.
"""
resetdetector!(system::CSGOpticalSystem{T}) where {T<:Real} = reset!(system.detectorimage)
export temperature, pressure, detectorimage, resetdetector!, assembly


Base.Float32(a::T) where {T<:ForwardDiff.Dual} = Float32(ForwardDiff.value(a))

"""
    trace(system::AbstractOpticalSystem{T}, ray::OpticalRay{T}; trackrays = nothing, test = false)

Traces `system` with `ray`, if `test` is enabled then fresnel reflections are disabled and the power distribution will not be correct.
Returns either a [`LensTrace`](@ref) if the ray hits the detector or `nothing` otherwise.

`trackrays` can be passed an empty vector to accumulate the `LensTrace` objects at each intersection of `ray` with a surface in the system.
"""
function trace(system::CSGOpticalSystem{T,D}, r::OpticalRay{T,N}; trackrays::Union{Nothing,Vector{LensTrace{T,N}}} = nothing, test::Bool = false) where {T<:Real,N,D<:Number}
    if power(r) < POWER_THRESHOLD
        return nothing
    end

    result = trace(system.assembly, r, temperature(system), pressure(system), trackrays = trackrays, test = test)

    if result === nothing || result === nopower
        emptyintervalpool!(T)
        return nothing
    else #ray intersected lens assembly so continue to see if ray intersects detector
        intsct = surfaceintersection(detector(system), ray(result))
        if intsct === nothing # no intersection of final ray with detector
            emptyintervalpool!(T)
            return nothing
        end

        detintsct = closestintersection(intsct)
        if detintsct === nothing
            emptyintervalpool!(T)
            return nothing
        else
            # need to modify power and path length accordingly for the intersection with the detector
            surfintsct = point(detintsct)
            nml = normal(detintsct)
            opticalinterface = interface(detintsct)::FresnelInterface{T}
            λ = wavelength(r)

            # Optical path length is measured in mm
            # in this case the result ray is exact so no correction for RAY_OFFSET is needed
            geometricpathlength = norm(surfintsct - origin(ray(result)))
            opticalpathlength = geometricpathlength
            pow = power(result)

            m = outsidematerialid(opticalinterface)
            # compute updated power based on absorption coefficient of material using Beer's law
            # this will almost always not apply as the detector will be in air, but it's possible that the detector is not in air, in which case this is necessary
            if !isair(m)
                mat = glassforid(m)::OpticSim.GlassCat.Glass
                nᵢ = index(mat, λ, temperature = temperature(system), pressure = pressure(system))::T
                α = absorption(mat, λ, temperature = temperature(system), pressure = pressure(system))::T
                if α > zero(T)
                    internal_trans = exp(-α * geometricpathlength)
                    if rand() >= internal_trans
                        return nothing
                    end
                    pow = pow * internal_trans
                end
                opticalpathlength = nᵢ * geometricpathlength
            end

            temp = LensTrace{T,N}(OpticalRay(ray(ray(result)), pow, wavelength(result), opl = pathlength(result) + opticalpathlength, nhits = nhits(result) + 1, sourcenum = sourcenum(r), sourcepower = sourcepower(r)), detintsct)
            if trackrays !== nothing
                push!(trackrays, temp)
            end

            # incerement the detector image
            pixu, pixv = uvtopix(detector(system), uv(detintsct), size(system.detectorimage))
            system.detectorimage[pixv, pixu] += D(sourcepower(r)) # TODO will need to handle different detector image types a bit better than this

            # should be okay to assume intersection will not be a DisjointUnion for all the types of detectors we will be using.
            emptyintervalpool!(T)
            return temp
        end
    end
end

######################################################################################################################

"""
    AxisymmetricOpticalSystem{T,C<:CSGOpticalSystem{T}} <: AbstractOpticalSystem{T}

Optical system which has lens elements and an image detector, created from a `DataFrame` containing prescription data.

These tags are supported for columns: `:Radius`, `:SemiDiameter`, `:Surface`, `:Thickness`, `:Conic`, `:Aspherics`, `:Reflectance`, `:Material`, `:OptimizeRadius`, `:OptimizeThickness`, `:OptimizeConic`.
These tags are supported for entries in a `:Surface` column: `:Object`, `:Image`, `:Stop`
Assumes the `:Image` row will be the last row in the `DataFrame`.

In practice a [`CSGOpticalSystem`](@ref) is generated automatically and stored within this system.

```julia
AxisymmetricOpticalSystem{T}(prescription::DataFrame, detectorpixelsx = 1000, detectorpixelsy:: = 1000, ::Type{D} = Float32; temperature = OpticSim.GlassCat.TEMP_REF, pressure = OpticSim.GlassCat.PRESSURE_REF)
```
"""
struct AxisymmetricOpticalSystem{T,C<:CSGOpticalSystem{T}} <: AbstractOpticalSystem{T}
    system::C
    prescription::DataFrame
    semidiameter::T

    function AxisymmetricOpticalSystem{T}(prescription::DataFrame, detectorpixelsx::Int = 1000, detectorpixelsy::Int = 1000, ::Type{D} = Float32; temperature::Union{T,Unitful.Temperature} = T(OpticSim.GlassCat.TEMP_REF), pressure::T = T(OpticSim.GlassCat.PRESSURE_REF)) where {T<:Real,D<:Number}
        pr = prescription
        elements = Array{Union{Surface{T},CSGTree{T}},1}(undef, 0)
        (rows, cols) = size(pr)
        systemsemidiameter = zero(T)

        let sumstart = 1
            if pr[1, :Surface] == :Object
                sumstart = 2
            end
            firstelement = true

            for i in 1:(rows - 1)
                if pr[i, :Surface] == :Stop
                    stop = CircularAperture(T(pr[i, :SemiDiameter]), SVector{3,T}(0.0, 0.0, 1.0), SVector{3,T}(0.0, 0.0, T(-sum(pr[sumstart:(i - 1), :Thickness]))))
                    push!(elements, stop)
                else
                    material = pr[i, :Material]
                    frontreflection = zero(T)
                    backreflection = zero(T)
                    frontconic = zero(T)
                    backconic = zero(T)
                    frontaspherics = nothing
                    backaspherics = nothing
                    if material != OpticSim.GlassCat.Air && !isequal(material, missing)
                        semidiameter::T = max(pr[i, :SemiDiameter], pr[i + 1, :SemiDiameter])
                        if firstelement
                            systemsemidiameter = semidiameter
                            firstelement = false
                        end

                        vertex::T = -sum(pr[sumstart:(i - 1), :Thickness])
                        frontradius::T = pr[i, :Radius]
                        backradius::T = pr[i + 1, :Radius]
                        if "Reflectance" in names(prescription)
                            temp = pr[i, :Reflectance]
                            if !ismissing(temp)
                                frontreflection = temp
                            end
                            temp = pr[i + 1, :Reflectance]
                            if !ismissing(temp)
                                backreflection = temp
                            end
                        end
                        if "Conic" in names(prescription)
                            temp = pr[i, :Conic]
                            if !ismissing(temp)
                                frontconic = temp
                            end
                            temp = pr[i + 1, :Conic]
                            if !ismissing(temp)
                                backconic = temp
                            end
                        end
                        if "Aspherics" in names(prescription)
                            temp = pr[i, :Aspherics]
                            if !ismissing(temp)
                                frontaspherics = temp
                            end
                            temp = pr[i + 1, :Aspherics]
                            if !ismissing(temp)
                                backaspherics = temp
                            end
                        end

                        lastmaterial = pr[i - 1, :Material]
                        nextmaterial = pr[i + 1, :Material]

                        thickness = T(pr[i, :Thickness])

                        if frontaspherics !== nothing || backaspherics !== nothing
                            newelt = AsphericLens(material, vertex, frontradius, frontconic, frontaspherics, backradius, backconic, backaspherics, thickness, semidiameter, lastmaterial = lastmaterial, nextmaterial = nextmaterial, frontsurfacereflectance = frontreflection, backsurfacereflectance = backreflection)
                        elseif frontconic != zero(T) || backconic != zero(T)
                            newelt = ConicLens(material, vertex, frontradius, frontconic, backradius, backconic, thickness, semidiameter, lastmaterial = lastmaterial, nextmaterial = nextmaterial, frontsurfacereflectance = frontreflection, backsurfacereflectance = backreflection)
                        else
                            newelt = SphericalLens(material, vertex, frontradius, backradius, thickness, semidiameter, lastmaterial = lastmaterial, nextmaterial = nextmaterial, frontsurfacereflectance = frontreflection, backsurfacereflectance = backreflection)
                        end

                        push!(elements, newelt())
                    end
                end
            end

            indexofimage = findfirst(isequal(:Image), prescription[!, :Surface])
            imagesize = prescription[indexofimage, :SemiDiameter]
            vertextoimage = T(-sum(pr[sumstart:(indexofimage - 1), :Thickness]))
            imagerad = prescription[indexofimage, :Radius]

            if imagerad != zero(T) && imagerad != typemax(T)
                det = SphericalCap(abs(imagerad), NaNsafeasin(imagesize / abs(imagerad)), imagerad < 0 ? SVector{3,T}(0, 0, 1) : SVector{3,T}(0, 0, -1), SVector{3,T}(0, 0, vertextoimage), interface = opaqueinterface(T))
            else
                det = Rectangle(T(imagesize), T(imagesize), SVector{3,T}(0, 0, 1), SVector{3,T}(0, 0, vertextoimage), interface = opaqueinterface(T))
            end

            system = CSGOpticalSystem(OpticSim.LensAssembly(elements...), det, detectorpixelsx, detectorpixelsy, D, temperature = temperature, pressure = pressure)

            return new{T,typeof(system)}(system, prescription, systemsemidiameter)
        end
    end

    AxisymmetricOpticalSystem(prescription::DataFrame) = AxisymmetricOpticalSystem{Float64}(prescription)
end
export AxisymmetricOpticalSystem

Base.show(io::IO, a::AxisymmetricOpticalSystem) = print(io, a.prescription)
Base.copy(a::AxisymmetricOpticalSystem) = AxisymmetricOpticalSystem(a.prescription, size(detectorimage(a))..., temperature = temperature(a), pressure = pressure(a))

trace(system::AxisymmetricOpticalSystem{T,C}, r::OpticalRay{T,N}; trackrays::Union{Nothing,Vector{LensTrace{T,N}}} = nothing, test::Bool = false) where {T<:Real,N,C<:CSGOpticalSystem{T}} = trace(system.system, r, trackrays = trackrays, test = test)

"""
    semidiameter(system::AxisymmetricOpticalSystem{T}) -> T

Get the semidiameter of `system`, that is the semidiameter of the entrance pupil (i.e. first surface) of the system.
"""
semidiameter(a::AxisymmetricOpticalSystem) = a.semidiameter
assembly(system::AxisymmetricOpticalSystem) = assembly(system.system)
detector(system::AxisymmetricOpticalSystem) = detector(system.system)
detectorimage(system::AxisymmetricOpticalSystem) = detectorimage(system.system)
detectorsize(system::AxisymmetricOpticalSystem) = detectorsize(system.system)
resetdetector!(system::AxisymmetricOpticalSystem) = resetdetector!(system.system)
temperature(system::AxisymmetricOpticalSystem) = temperature(system.system)
pressure(system::AxisymmetricOpticalSystem) = pressure(system.system)
export semidiameter

######################################################################################################################

trace(system::AxisymmetricOpticalSystem{T}, raygenerator::OpticalRayGenerator{T}; printprog::Bool = true, test::Bool = false, outpath::Union{Nothing,String} = nothing) where {T<:Real} = trace(system.system, raygenerator, printprog = printprog, test = test, outpath = outpath)

"""
    trace(system::AbstractOpticalSystem{T}, raygenerator::OpticalRayGenerator{T}; printprog = true, test = false)

Traces `system` with rays generated by `raygenerator` on a single thread.
Optionally the progress can be printed to the REPL.
If `test` is enabled then fresnel reflections are disabled and the power distribution will not be correct.
If `outpath` is specified then the result will be saved to this path.

Returns the detector image of the system.
"""
function trace(system::CSGOpticalSystem{T}, raygenerator::OpticalRayGenerator{T}; printprog::Bool = true, test::Bool = false, outpath::Union{Nothing,String} = nothing) where {T<:Real}
    start_time = time()
    update_timesteps = 1000
    total_traced = 0
    for (k, r) in enumerate(raygenerator)
        if k % update_timesteps == 0
            total_traced += update_timesteps
            dif = round(time() - start_time, digits = 1)
            left = round((time() - start_time) * (length(raygenerator) / total_traced - 1), digits = 1)
            if printprog
                print("\rTraced: ~ $t / $(length(raygenerator))        Elapsed: $(dif)s        Left: $(left)s           ")
            end
        end
        trace(system, r, test = test)
    end
    numrays = length(raygenerator)
    tracetime = round(time() - start_time)
    if printprog
        print("\rFinished tracing $numrays rays in $(tracetime)s")
        if tracetime != 0.0
            print(",  $(Int32(round(numrays/tracetime))) rays per second")
        end
        println()
    end

    det = detectorimage(system)

    if outpath !== nothing
        save(outpath, colorview(Gray, det ./ maximum(det)))
    end

    return det
end


function traceMT(system::AxisymmetricOpticalSystem{T}, raygenerator::OpticalRayGenerator{T}; printprog::Bool = true, test::Bool = false, outpath::Union{Nothing,String} = nothing) where {T<:Real}
    traceMT(system.system, raygenerator, printprog = printprog, test = test, outpath = outpath)
end

"""
    traceMT(system::AbstractOpticalSystem{T}, raygenerator::OpticalRayGenerator{T}; printprog = true, test = false)

Traces `system` with rays generated by `raygenerator` using as many threads as possible.
Optionally the progress can be printed to the REPL.
If `test` is enabled then fresnel reflections are disabled and the power distribution will not be correct.
If `outpath` is specified then the result will be saved to this path.

Returns the accumulated detector image from all threads.
"""
function traceMT(system::CSGOpticalSystem{T,S}, raygenerator::OpticalRayGenerator{T}; printprog::Bool = true, test::Bool = false, outpath::Union{Nothing,String} = nothing) where {T<:Real,S<:Number}
    if printprog
        println("Initialising...")
    end

    all_start_time = time()
    N = Threads.nthreads()
    copies = Vector{Tuple{CSGOpticalSystem{T},OpticalRayGenerator{T}}}(undef, N)
    for i in eachindex(copies)
        copies[i] = (copy(system), copy(raygenerator))
    end

    if printprog
        println("Tracing...")
    end
    total_traced = Threads.Atomic{Int}(0)
    update_timesteps = 1000
    Threads.@threads for sys in copies
        (system, sources) = sys
        i = Threads.threadid()
        len, rem = divrem(length(sources), N)
        # not enough iterations for all the threads
        if len == 0
            if i > rem
                return
            end
            len, rem = 1, 0
        end
        # compute this thread's iterations
        f = firstindex(sources) + ((i - 1) * len)
        l = f + len - 1
        # distribute remaining iterations evenly
        if rem > 0
            if i <= rem
                f = f + (i - 1)
                l = l + i
            else
                f = f + rem
                l = l + rem
            end
        end
        if i == 1
            start_time = time()
            for k in f:l
                if k % update_timesteps == 0
                    Threads.atomic_add!(total_traced, update_timesteps)
                    dif = round(time() - start_time, digits = 1)
                    t = total_traced[]
                    left = round((time() - start_time) * (length(sources) / t - 1), digits = 1)
                    if printprog
                        print("\rTraced: ~ $t / $(length(sources))        Elapsed: $(dif)s        Left: $(left)s           ")
                    end
                end
                trace(system, sources[k]; test)
            end
        else
            for k in f:l
                if k % update_timesteps == 0
                    Threads.atomic_add!(total_traced, update_timesteps)
                end
                trace(system, sources[k]; test)
            end
        end
    end

    if printprog
        println("\rAccumulating images...                                                                       ")
    end

    det = OpticSim.detectorimage(copies[1][1])
    # sum all the copies
    @inbounds for i in 2:N
        OpticSim.sum!(det, OpticSim.detectorimage(copies[i][1]))
    end

    numrays = length(raygenerator)
    tracetime = round(time() - all_start_time)
    if printprog
        print("\rFinished tracing $numrays rays in $(tracetime)s,  ")
        if tracetime != 0.0
            print("$(Int32(round(numrays/tracetime))) rays per second")
        end
        println()
    end

    if outpath !== nothing
        save(outpath, colorview(Gray, det ./ maximum(det)))
    end

    return det
end

function tracehitsMT(system::AxisymmetricOpticalSystem{T}, raygenerator::OpticalRayGenerator{T}; printprog::Bool = true, test::Bool = false) where {T<:Real}
    tracehitsMT(system.system, raygenerator, printprog = printprog, test = test)
end

"""
    tracehitsMT(system::AbstractOpticalSystem{T}, raygenerator::OpticalRayGenerator{T}; printprog = true, test = false)

Traces `system` with rays generated by `raygenerator` using as many threads as possible.
Optionally the progress can be printed to the REPL.
If `test` is enabled then fresnel reflections are disabled and the power distribution will not be correct.

Returns a list of [`LensTrace`](@ref)s which hit the detector, accumulated from all threads.
"""
function tracehitsMT(system::CSGOpticalSystem{T}, raygenerator::OpticalRayGenerator{T}; printprog::Bool = true, test::Bool = false) where {T<:Real}
    if printprog
        println("Initialising...")
    end

    all_start_time = time()
    N = Threads.nthreads()
    copies = Vector{Tuple{CSGOpticalSystem{T},OpticalRayGenerator{T}}}(undef, N)
    results = Vector{Vector{LensTrace{T,3}}}(undef, N)
    for i in eachindex(copies)
        copies[i] = (copy(system), copy(raygenerator))
        results[i] = Vector{LensTrace{T,3}}(undef, 0)
    end

    if printprog
        println("Tracing...")
    end
    total_traced = Threads.Atomic{Int}(0)
    update_timesteps = 1000
    Threads.@threads for sys in copies
        (system, sources) = sys
        i = Threads.threadid()
        len, rem = divrem(length(sources), N)
        # not enough iterations for all the threads
        if len == 0
            if i > rem
                return
            end
            len, rem = 1, 0
        end
        # compute this thread's iterations
        f = firstindex(sources) + ((i - 1) * len)
        l = f + len - 1
        # distribute remaining iterations evenly
        if rem > 0
            if i <= rem
                f = f + (i - 1)
                l = l + i
            else
                f = f + rem
                l = l + rem
            end
        end
        if i == 1
            start_time = time()
            for k in f:l
                if k % update_timesteps == 0
                    Threads.atomic_add!(total_traced, update_timesteps)
                    dif = round(time() - start_time, digits = 1)
                    t = total_traced[]
                    left = round((time() - start_time) * (length(sources) / t - 1), digits = 1)
                    if printprog
                        print("\rTraced: ~ $t / $(length(sources))        Elapsed: $(dif)s        Left: $(left)s           ")
                    end
                end
                lt = trace(system, sources[k]; test)
                if lt !== nothing
                    push!(results[i], lt)
                end
            end
        else
            for k in f:l
                if k % update_timesteps == 0
                    Threads.atomic_add!(total_traced, update_timesteps)
                end
                lt = trace(system, sources[k]; test)
                if lt !== nothing
                    push!(results[i], lt)
                end
            end
        end
    end

    if printprog
        println("\rAccumulating results...                                                                       ")
    end

    accumres = results[1]
    # sum all the copies
    @inbounds for i in 2:N
        accumres = vcat(accumres, results[i])
    end

    numrays = length(raygenerator)
    tracetime = round(time() - all_start_time)
    if printprog
        print("\rFinished tracing $numrays rays in $(tracetime)s")
        if tracetime != 0.0
            print(",  $(Int32(round(numrays/tracetime))) rays per second")
        end
        println()
    end

    return accumres
end


tracehits(system::AxisymmetricOpticalSystem{T}, raygenerator::OpticalRayGenerator{T}; printprog::Bool = true, test::Bool = false) where {T<:Real} = tracehits(system.system, raygenerator, printprog = printprog, test = test)

"""
    tracehits(system::AbstractOpticalSystem{T}, raygenerator::OpticalRayGenerator{T}; printprog = true, test = false)

Traces `system` with rays generated by `raygenerator` on a single thread.
Optionally the progress can be printed to the REPL.
If `test` is enabled then fresnel reflections are disabled and the power distribution will not be correct.

Returns a list of [`LensTrace`](@ref)s which hit the detector.
"""
function tracehits(system::CSGOpticalSystem{T}, raygenerator::OpticalRayGenerator{T}; printprog::Bool = true, test::Bool = false) where {T<:Real}
    start_time = time()
    update_timesteps = 1000
    total_traced = 0
    res = Vector{LensTrace{T,3}}(undef, 0)
    for (k, r) in enumerate(raygenerator)
        if k % update_timesteps == 0
            total_traced += update_timesteps
            dif = round(time() - start_time, digits = 1)
            left = round((time() - start_time) * (length(sources) / total_traced - 1), digits = 1)
            if printprog
                print("\rTraced: ~ $t / $(length(sources))        Elapsed: $(dif)s        Left: $(left)s           ")
            end
        end
        lt = trace(system, r, test = test)
        if lt !== nothing
            push!(res, lt)
        end
    end
    numrays = length(raygenerator)
    tracetime = round(time() - all_start_time)
    if printprog
        print("\rFinished tracing $numrays rays in $(tracetime)s")
        if tracetime != 0.0
            print(",  $(Int32(round(numrays/tracetime))) rays per second")
        end
        println()
    end

    return res
end


export trace, traceMT, tracehits, tracehitsMT
