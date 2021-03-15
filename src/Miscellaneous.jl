module Junk
# TODO clean this up

using ..Optics
using BenchmarkTools
import Unitful

# space for junk code which should ultimately be deleted or cleaned up and moved elsewhere

function drawbezier()
    curve = TestData.onespanspline()
    drawcurve(curve, 100, 1000)
end

function drawbsplinecurve()
    curve = TestData.bsplinecurve()
    drawcurve(curve, 100, 1000)
end

function drawbeziersurface()
    vx = 0.0:0.1:1
    vy = 0.0:0.1:1
    surface = TestData.beziersurface()

    f(x, y) = point(surface, x, y)[3]
    scene = Makie.Scene(resolution = (1000, 1000))
    # One way to style the axis is to pass a nested dictionary / named tuple to it.

    # for i in 1:10000
    #     surface.controlpolygon.controlpoints[2,2].point[3] = i/100.0
    #     Makie.wireframe!(vx,vy,zvals)
    # end
    # record(scene, "surface.mp4", 1:255; framerate = 60) do i
    # surface.controlpolygon[2,2][3] = Float64(i)/40.0
    zvals = [f(x, y) for x in vx, y in vy]
    # Makie.wireframe!(vx,vy,zvals,axis = (frame = (linewidth = 2.0,),))
    # Makie.surface!(scene, vx, vy, f, axis = (frame = (linewidth = 2.0,),))
    Makie.wireframe!(vx, vy, zvals)
    # end
end

function plotcurveintersections()
    curve = TestData.curvecoefficients()
    maxtheta = 1.01

    let p
        curvepts = collect((evaluatecurve(curve, theta) for theta in 0.0:0.01:maxtheta))
        let x = Array{Real,1}(undef, 0), y = Array{Real,1}(undef, 0)
            for pt in curvepts
                push!(x, pt[1])
                push!(y, pt[2])
            end
            p = Plots.plot(x, y)
        end

        origin = [0, 0]

        for theta in 0.01:0.2:maxtheta
            endpoint = evaluatecurve(curve, theta)
            line = Ray([0.0, 0.0], endpoint)
            points = intersections(line, curve)
            maxalpha = 0

            for int in points
                pt = int[1]
                theta = int[3]

                if 0 <= theta <= maxtheta # only plot intersection points in the curve parameter range from 0..maxtheta
                    if maxalpha < int[2]
                        maxalpha = int[2]
                    end

                    ptx = [pt[1]]
                    pty = [pt[2]]
                    Plots.scatter!(p, ptx, pty, legend = false)
                end
            end

            line = hcat([0, 0], endpoint .* maxalpha)
            x = line'[:, 1]
            y = line'[:, 2]
            Plots.plot!(p, x, y)
        end
        Plots.display(p)
    end
end

function makiemeshexample5()
    points = [(0, 0, 0), (1, 0, 2), (2, 1, 0), (0, 2, 1)]
    color = [:red, :green, :blue, :yellow]

    # indices interpreted as triangles (every 3 sequential indices)
    indices = [1 2 3; 1 3 4; 1 4 2; 2 3 4]
    for i in 1:1
        Makie.mesh!(points, indices, color = color)
    end
end

function makiemeshexample4()
    points = [(0, 0, 0), (1, 0, 2), (2, 1, 0), (0, 2, 1)]
    color = [:red, :green, :blue, :yellow]

    # indices interpreted as triangles (every 3 sequential indices)
    indices = [1 2 3; 1 3 4; 1 4 2; 2 3 4]
    # replprint(indices)
    # for i in 1:1
    scene = Makie.Scene(resolution = (1000, 1000))
    Makie.mesh!(points, indices, color = color)
    scene
    # end
end

function makiemeshexample3()
    points = [(0, 0, 0), (1, 0, 2), (2, 1, 0), (0, 2, 1)]
    color = [:red, :green, :blue, :yellow]

    # indices interpreted as triangles (every 3 sequential indices)
    indices = [1, 2, 3, 1, 3, 4, 1, 4, 2, 2, 3, 4]
    Makie.mesh(points, indices, color = color)
end

function makiemeshexample2()
    x = [0, 1, 2, 0]
    y = [0, 0, 1, 2]
    z = [0, 2, 0, 1]
    color = [:red, :green, :blue, :yellow]
    i = [0, 0, 0, 1]
    j = [1, 2, 3, 2]
    k = [2, 3, 1, 3]
    # indices interpreted as triangles (every 3 sequential indices)
    indices = [1, 2, 3, 1, 3, 4, 1, 4, 2, 2, 3, 4]
    Makie.mesh(x, y, z, indices, color = color)
end

function makiemeshexample()
    coordinates = [
        0.0 0.0
        0.5 0.0
        1.0 0.0
        0.0 0.5
        0.5 0.5
        1.0 0.5
        0.0 1.0
        0.5 1.0
        1.0 1.0
    ]
    connectivity = [
        1 2 5
        1 4 5
        2 3 6
        2 5 6
        4 5 8
        4 7 8
        5 6 9
        5 8 9
    ]
    color = [0.0, 0.0, 0.0, 0.0, -0.375, 0.0, 0.0, 0.0, 0.0]
    scene = Makie.mesh(coordinates, connectivity, color = color, shading = false)

    Makie.wireframe!(scene[end][1], color = (:black, 0.6), linewidth = 3)
end




# function simplerender(samples::Int = 100)
#     cyl = Cylinder(0.3, 2.0)
#     cyl2 = Cylinder(0.2, 2.0)
#     rbt = translation(0.5, 0.5, 0.0)
#     # generator = csgintersection(leaf(cyl),leaf(cyl2,rbt))
#     generator = csgunion(leaf(cyl), leaf(cyl2, rbt))
#     # generator = leaf(cyl2,rbt)

#     # surf = rbt*AcceleratedParametricSurface(TestData.beziersurface(),10)
#     # generator = csgintersection(leaf(cyl),leaf(surf))

#     csg = generator(identitytransform())
#     view = Vis.CSGViewer([3.0, 3.0, 20.0], 2.0, samples, (csg,))
#     image = Vis.render(view)
#     Vis.show(image)
# end

# function csgrender(samples::Int = 100)
#     # generator = TestData.planecylindercsg()
#     # csg = generator(identitytransform())

#     # generator = TestData.trianglemeshcsg()
#     # csg = generator(RigidBodyTransform(-.5, -.5, 0.0))
#     # generator = TestData.csgobject()
#     # csg = generator(RigidBodyTransform(-.5, -.5, 0.0))

#     generator = TestData.lensarray(2)
#     csg = generator(identitytransform())
#     # view = CSGViewer(SVector{3}(0.0, 60.0, 80.0), .5, samples, (csg,))
#     # generator = leaf(AcceleratedParametricSurface(beziersurface(),1

#     # generator = leaf(Cylinder(.3,10.0))
#     # csg = generator(RigidBodyTransform(0.0,0.0,0.0))

#     # csg = AcceleratedParametricSurface(beziersurface(), 10)
#     view = Vis.CSGViewer([0.0, 4.0, 8.0], 1.4, samples, (csg,))
#     image = Vis.render(view)

#     # Vis.show(image)
# end

function timebezierpoint()
    surf = TestData.beziersurface()
    res = @benchmark point($surf, 0.5, 0.5)
    show(res)
end

function timereflectedray()
    res = @benchmark TestData.reflectedray(rand(3), rand(3))
    println(" Array{Float64,1} time")
    show(res)
    res = @benchmark TestData.reflectedray(SVector{3,Float64}(rand(3)), SVector{3,Float64}(rand(3)))
    println("SVector time")
    show(res)
end

function timebezierrayintersection()
    u = 0.5
    v = 0.5
    surf, r = TestData.beziersurfaceandray(u, v)
    surf = AcceleratedParametricSurface(surf, 5)
    bnch = @benchmark surfaceintersection($surf, $r)
    println(surfaceintersection(surf, r), " =? ", point(surf, u, v))
    println()
    show(bnch)
end

function timebboxrayintersection()
    bbox, r = TestData.bboxray()
    println(doesintersect(bbox, r))
    res = @benchmark doesintersect($bbox, $r)
    show(res)
end

function timecylinderrayintersection()
    cyl = Cylinder(0.5)
    r = Ray([1.0, 1.0, 1.0], [-1.0, -1.0, -1.0])
    bnch = @benchmark surfaceintersection($cyl, $r)
    println(surfaceintersection(cyl, r))
    println()
    show(bnch)
end

function timesphererayintersection()
    sph = Sphere(0.5)
    r = Ray([1.0, 1.0, 1.0], [-1.0, -1.0, -1.0])
    bnch = @benchmark surfaceintersection($sph, $r)
    println(surfaceintersection(sph, r))
    println()
    show(bnch)
end

function timetrianglerayintersection()
    tri, r = TestData.triangleray()
    res = @benchmark surfaceintersection($tri, $r)
    println(surfaceintersection(tri, r))
    println()
    show(res)
end

function bboxrayintersection(numits)
    bbox, r = TestData.bbox()

    for i in 1:numits
        res = surfaceintersection(bbox, r)
    end
end

function timeplanerayintersection()
    pln, r = TestData.planeray()
    res = @benchmark surfaceintersection($pln, $r)
    println(surfaceintersection(pln, r))
    println()
    show(res)
end

function onaxissphereintersection()
    println(surfaceintersection(Sphere(1.0), Ray([0.0, 0.0, -4.0], [0.0, 0.0, 1.0])))
    println(surfaceintersection(Sphere(1.0), Ray([0.0, 0.0, 4.0], [0.0, 0.0, -1.0])))
    println(surfaceintersection(Sphere(1.0), Ray([0.0, 4.0, 0.0], [0.0, -1.0, 0.0])))
    println(surfaceintersection(Sphere(1.0), Ray([0.0, -4.0, 0.0], [0.0, 1.0, 0.0])))
    println(surfaceintersection(Sphere(1.0), Ray([4.0, 0.0, 0.0], [-1.0, 0.0, 0.0])))
    println(surfaceintersection(Sphere(1.0), Ray([-4.0, 0.0, 0.0], [1.0, 0.0, 0.0])))
end

function makeimagedata(samples)
    zeropix = (0.0f0, 0.0f0, 0.0f0)
    return (fill(zeropix, samples, samples), fill(zeropix, samples, samples))
end

function threadingtest(samples)
    image, image2 = makeimagedata(samples)

    function innerloop(x, image, image2, samples, convolutionsize)
        for y in (convolutionsize + 1):samples
            for l in 1:convolutionsize
                for m in 1:convolutionsize
                    pix = image[y - l, x - m]
                    oldpix = image2[y, x]
                    image2[y, x] = (oldpix[1] + pix[1] * 0.8, oldpix[2] + pix[2] * 0.8, oldpix[3] + pix[3] * 0.8)
                end
            end
        end
    end

    convolutionsize = 30
    # println("multithreaded")
    # @sync
    for x in (convolutionsize + 1):samples
        #  for x in convolutionsize + 1:samples
        #    Threads.@spawn
        innerloop(x, image, image2, samples, convolutionsize)
    end
    return image
end

function timeconvolution()
    result = @benchmark threadingtest(1000)
    show(result)
end

function testdrawmesh()
    Vis.draw(TestData.beziersurface())
end

function testintervalsforcsg()
    surf1 = AcceleratedParametricSurface(TestData.beziersurface(), 5)
    surf2 = AcceleratedParametricSurface(TestData.upsidedownbeziersurface(), 5)
    r = Ray([0.5, 0.5, 2.0], [0.0, 0.0, -1.0])
    int1 = DisjointUnion(surfaceintersection(surf1, r)...)
    int2 = DisjointUnion(surfaceintersection(surf2, r)...)
    intsct = intervalintersection(int1, int2)
    if !(intsct isa EmptyInterval)
        s = Vis.scene()
        Vis.draw!(s, surf1)
        Vis.draw!(s, surf2)
        Vis.draw!(s, intsct)
        Vis.display(s)
    else
        throw(ErrorException("This should intersect (testintervalsforcsg)"))
    end
end

function testclippedbeziersurface()
    r = Ray([0.5, 0.5, 3.0], [0.0, 0.0, -1.0])
    gen = TestData.clippedbeziersurface()
    csg = gen(identitytransform())
    intsct = Optics.evalcsg(csg, r)
    # println(intsct)
    view = Vis.CSGViewer(SVector{3}(2.0, 2.0, 5.0), 2.0, 1000, (csg,))
    Vis.imshow(Vis.render(view))
end

function testdrawcylinder()
    surf = Cylinder(1.0)
    s = Vis.scene()
    Vis.draw!(s, surf, numdivisions = 100)
    r = Ray([0.0, 2.0, 0.0], [0.0, -1.0, 0.5])
    int = surfaceintersection(surf, r)
    Vis.draw!(s, int)
    Vis.display(s)
end

function testdrawsphere()
    surf = Sphere(1.0)
    s = Vis.scene()
    Vis.draw!(s, surf, numdivisions = 100)
    r = Ray([0.0, 2.0, 0.0], [0.0, -1.0, 0.5])
    int = surfaceintersection(surf, r)
    Vis.draw!(s, int)
    Vis.display(s)
end

function testcsgeval()
    surf1 = AcceleratedParametricSurface(TestData.beziersurface(), 5)
    surf2 = AcceleratedParametricSurface(TestData.upsidedownbeziersurface(), 5)
    r = Ray([0.5, 0.5, 2.0], [0.0, 0.0, -1.0])
    gen = csgintersection(leaf(surf1), leaf(surf2))
    csgtree = gen(identitytransform())
    intsct = surfaceintersection(csgtree, r)
    if !(intsct isa EmptyInterval)
        s = Vis.scene()
        Vis.draw!(s, surf1)
        Vis.draw!(s, surf2)
        Vis.draw!(s, intsct)
        Vis.display(s)
    else
        throw(ErrorException("This should intersect (testcsgeval)"))
    end
end

function testcsgwithmultipleobjects()
    generator, _, surf1, surf2 = TestData.csgobject()
    csg = generator(identitytransform())
    r = Ray{Float64,3}(SVector(0.5, 0.5, 4.0), SVector(0.0, 0.0, -1.0))

    intscts = Optics.evalcsg(csg, r)[1]
    # println(intscts)
    s = Vis.scene()
    Vis.draw!(s, surf1)
    Vis.draw!(s, surf2)
    Vis.draw!(s, intscts)
    Vis.display(s)
end

function testmeshcsgsurf1(subdiv = 20)
    _, _, surf, _ = TestData.csgobject()
    m = csgintersection(leaf(surf, translation(-0.5, -0.5, 0.0)), csgintersection(leaf(Cylinder(0.3, 5.0)), leaf(surf, RigidBodyTransform(rotmatd(0, 180, 0), SVector(0.5, -0.5, 0.0)))))()
    b = @benchmark makemesh($m, $subdiv)
    Vis.draw(m)
    return b
end

function testmeshcsgsurf2(subdiv = 20)
    _, _, surf1, surf2 = TestData.csgobject()
    m = csgintersection(leaf(surf1, translation(-0.5, -0.5, 0.0)), csgintersection(leaf(Cylinder(0.3, 5.0)), leaf(surf2, translation(-0.5, -0.5, 0.0))))()
    b = @benchmark makemesh($m, $subdiv)
    Vis.draw(m)
    return b
end

function testmeshcsguni(subdiv = 20)
    m = csgunion(leaf(Cylinder(0.5, 3.0)), (leaf(Sphere(1.0))))()
    b = @benchmark makemesh($m, $subdiv)
    Vis.draw(m)
    return b
end

function testmeshcsgint(subdiv = 20)
    m = csgintersection(leaf(Cylinder(0.5, 3.0)), (leaf(Sphere(1.0))))()
    b = @benchmark makemesh($m, $subdiv)
    Vis.draw(m)
    return b
end

function testmeshcsgsub(subdiv = 20)
    m = csgdifference(leaf(Cylinder(0.5, 3.0)), leaf(Sphere(1.0)))()
    b = @benchmark makemesh($m, $subdiv)
    Vis.draw(m)
    return b
end

function testmeshcsgmla(subdiv = 20)
    topsphere = leaf(Sphere(2.0), translation(0.05, 0.0, -1.5))
    botsphere = leaf(Sphere(3.0), translation(0.0, 0.1, 1.5))
    lens = csgintersection(leaf(Cylinder(0.7, 5.0)), csgintersection(botsphere, topsphere))
    d = 0.8
    MLA = lens
    for x in 0:2
        for y in 0:2
            if x == 0 && y == 0
                continue
            end
            MLA = csgunion(MLA, leaf(lens, RigidBodyTransform(rotmatd(0, 0, rand(Int) % 180), SVector{3}(x * d, y * d, 0.0))))
        end
    end
    m = MLA(identitytransform())
    b = @benchmark makemesh($m, $subdiv)
    Vis.draw(m)
    return b
end

function testorthogonalitymatrix()
    x = [0 4 -4.5 4.5; 0 11.25 -31.5 20.25]
    correctCmatrix = [
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0
        4.0 0.0 0.0 0.0 11.25 0.0 0.0 0.0 0.0 1.0 0.0 0.0
        -4.5 4.0 0.0 0.0 -31.5 11.25 0.0 0.0 0.0 0.0 1.0 0.0
        4.5 -4.5 4.0 0.0 20.25 -31.5 11.25 0.0 0.0 0.0 0.0 1.0
        0.0 4.5 -4.5 4.0 0.0 20.25 -31.5 11.25 0.0 0.0 0.0 0.0
        0.0 0.0 4.5 -4.5 0.0 0.0 20.25 -31.5 0.0 0.0 0.0 0.0
        0.0 0.0 0.0 4.5 0.0 0.0 0.0 20.25 0.0 0.0 0.0 0.0
    ]
    cmat = Optics.orthogonalitymatrix(x, 3)
    show(IOContext(stdout), "text/plain", cmat)
    for i in CartesianIndices(cmat)
        if cmat[i] != correctCmatrix[i]
            println("did not pass test")
            break
        end
    end
end

function testtriangulate()
    surf = TestData.beziersurface()
    numtris = 20 # number of triangles per row or column
    tris = triangulate(surf, numtris)
    tmesh = TriangleMesh(tris)
    Vis.draw(tmesh)
end

function testcurvetobeziersegments()
    orig = TestData.homogeneousbsplinecurve()
    segments = tobeziersegments(orig)
    allsegments = [BezierCurve{Optics.Rational,Float64,3,3}(segment) for segment in segments]
    Vis.draw(orig, allsegments...)
end

function testsurfacetobeziersegments()
    orig = TestData.bsplinesurface()
    segments = tobeziersegments(orig)
    Vis.draw(segments...)
end

function testinsertknot()
    orig = TestData.bsplinecurve()
    curve1 = TestData.bsplinecurve()
    let curve2 = curve1
        for i in 1:2
            curve2 = insertknot(curve1, 5 + i - 1)
            curve1 = curve2
        end
        Vis.draw(orig, curve2)
    end
end

# function timelinalg(vec::Bool = true, mat::Bool = true)
#     # pf(q) = print(BenchmarkTools.prettytime(time(median(q))), "\n")

#     if vec
#         print("VECTOR OPS\n")
#         vec3a = @SVector rand(3)
#         vec3b = @SVector rand(3)
#         vec1000a = @SVector rand(1000)
#         vec1000b = @SVector rand(1000)

#         for (vec1, vec2) in [(vec3a, vec3b), (vec1000a, vec1000b)]
#             print("SIZE: $(length(vec1))\n")
#             # NORMS
#             print("NORMALIZE\n")
#             @btime (Optics.mnormalize($vec1))
#             @btime (LinearAlgebra.normalize($vec1))

#             # DOT
#             print("DOT\n")
#             @btime (Optics.dot($vec1, $vec2))
#             @btime (LinearAlgebra.dot($vec1, $vec2))

#             # CROSS
#             if length(vec1) == 3
#                 print("CROSS\n")
#                 @btime (Optics.cross($vec1, $vec2))
#                 @btime (LinearAlgebra.cross($vec1, $vec2))
#             end
#         end
#     end

#     if mat
#         print("MATRIX OPS\n")
#         mat3a = @SMatrix rand(3, 3)
#         mat3b = @SMatrix rand(3, 3)
#         mat100a = @SMatrix rand(12, 12)
#         mat100b = @SMatrix rand(12, 12)
#         vec100 = @SVector rand(12)

#         for (mat1, mat2) in [(mat3a, mat3b), (mat100a, mat100b)]
#             print("SIZE: $(size(mat1))\n")

#             # DETERMINANT
#             print("DET\n")
#             @btime (Optics.det($mat1))
#             @btime (LinearAlgebra.det($mat1))

#             # INVERSE
#             print("INV\n")
#             @btime (Optics.minv($mat1))
#             @btime (LinearAlgebra.inv($mat1))

#             # DIVISON
#             print("DIV\n")
#             @btime (Optics.mdiv($mat1, $mat2))
#             @btime ($mat1 \ $mat2)
#         end
#     end
# end

function testniandnt()
    r = Ray([0.0, 0.0, 7.0], [0.0, 0.0, -1.0])
    normal = SVector{3,Float64}(0.0, 0.0, 1.0)


end

function testrefractedray()
    elt = TestData.planoconvexelement(Optics.GlassCat.SCHOTT.N_BK7, 0.0, 5.0, 60.0)
    lens = Optics.LensAssembly(elt.lens)
    r = Ray([0.0, 0.0, 7.0], [0.0, 0.0, -1.0])

    green = 500 * Unitful.u"nm"
    glass = Optics.GlassCat.SCHOTT.BAK50

    intsct = Optics.closestintersection(lens, r)

    pt = point(intsct[1])
    nml = normal(intsct[1])

    rdir = direction(r)
    (nᵢ, nₜ) = Optics.nᵢandnₜ(index(Optics.GlassCat.Air, green), index(glass, green), nml, r)

    refracted = Optics.refractedray(nᵢ, nₜ, nml, rdir)
end

badray() = Ray([5.000000042187876, 5.000000042187876, 0.6143579891640621], [0.042187876699796796, 0.042187876699796796, -0.9982185963600987])

function printtrace(raytrace)
    for trace in raytrace
        println("RAY***************")
        println(trace)
        println("END RAY***********")
    end
end
function testdoubleconcave()
    system = TestData.doubleconcave()
    green = 500 * Unitful.u"nm"
    r1 = Ray([5.0, 5.0, 2.0], [0.0, 0.0, -1.0])
    r2 = badray()
    elt = system.system.assembly.elements[1]
    # intsct = Optics.surfaceintersection(elt,r)
    # println(intsct)
    for r in (r1,)
        allrays = Array{Optics.LensTrace{Float64,3},1}(undef, 0)

        res = Optics.trace(system, r, green, trackrays = allrays)
        printtrace(allrays)
    end
end

function testplane()
    rparallel = Ray([0.0, 0.0, 2.0], [1.0, 1.0, 0.0])
    rinplane = Ray([0.0, 0.0, 1.0], [1.0, 1.0, 0.0])
    routside = Ray([0.0, 0.0, 2.0], [1.0, 1.0, 1.0])
    rinside = Ray([0.0, 0.0, -1.0], [1.0, 1.0, -1.0])
    routintersects = Ray([0.0, 0.0, 2.0], [0.0, 0.0, -1.0])
    rinintersects = Ray([0.0, 0.0, -1.0], [0.0, 0.0, 1.0])
    pln = Plane(0.0, 0.0, 1.0, 0.0, 0.0, 1.0)

    surfaceintersection(pln, routside)
    if !(nothing === surfaceintersection(pln, rparallel))
        @error "should have returned nothing for ray parallel to plane and outside"
    else
        if !(nothing === surfaceintersection(pln, routside))
            @error "shouldhave returned nothing for entirely outside plane but not parallel"
        else
            # NOTE for coplanar faces to work with visualization, rays in the plane must count as being 'inside'
            res = surfaceintersection(pln, rinplane)
            if !(lower(res) isa RayOrigin && upper(res) isa Infinity)
                @error "should have returned (rayorigin,∞)"
                return false
            else
                res = surfaceintersection(pln, rinside)
                if !(lower(res) isa RayOrigin && upper(res) isa Infinity)
                    @error "should have returned (rayorigin,∞)"
                    return false
                else
                    res = surfaceintersection(pln, routintersects)
                    if !(samepoint(point(lower(res)), [0.0, 0.0, 1.0]) && upper(res) isa Infinity)
                        @error "should have returned ([0.0,0.0,0.0],∞)"
                        return false
                    else
                        res = surfaceintersection(pln, rinintersects)
                        if !(lower(res) isa RayOrigin && samepoint(point(upper(res)), [0.0, 0.0, 1.0]))
                            @error "should have returned(rayorigin,[0.0,0.0,0.0])"
                            return false
                        else
                            return true
                        end
                    end
                end
            end
        end
    end
end

function testintervalcomplement()
    intvl = surfaceintersection((leaf(Sphere(1.0)))(identitytransform()), Ray([2.0, 0.0, 0.0], [-1.0, 0.0, 0.0]))
    comp = Optics.intervalcomplement(intvl)
    println(comp)
end

function testcsgoperations()
    pln1 = Plane(SVector{3}(0.0, 0.0, -1.0), SVector{3}(0.0, 0.0, -1.0))
    pln2 = Plane(SVector{3}(0.0, 0.0, 1.0), SVector{3}(0.0, 0.0, 1.0))
    r = Ray{Float64,3}(SVector{3}(0.0, 0, 2.0), SVector{3}(0.0, 0.0, -1.0))
    gen1 = csgintersection(pln1, pln2)
    gen2 = csgunion(pln1, pln2)

    csg = gen1(identitytransform())
    intsct = surfaceintersection(csg, r)

    intvlpt1 = lower(intsct)
    intvlpt2 = upper(intsct)
    println(intvlpt1, "  ", intvlpt2)
    if !(isapprox(point(intvlpt1), SVector{3}(0.0, 0.0, 1.0)) && isapprox(point(intvlpt2), SVector{3}(0.0, 0.0, -1.0)))
        @error "Should have two intersections with the two planes"
        return false
    end
    csg = gen2(identitytransform())
    intsct = surfaceintersection(csg, r)
    if !((lower(intsct) isa RayOrigin) && (upper(intsct) isa Infinity))
        @error "Union space should be half infinite"
        return false
    end
    return true
end


function testinsideplane()
    r = Ray([0.0, 0.0, 1.0], [0.0, 0.0, -1.0])
    p = Plane(0.0, 0.0, 1.0, 0.0, 0.0, 1.5)
    res = surfaceintersection(p, r)
    println(res)
end

function testsinglenegsurface()
    frontradius = -60.0
    frontvertex = 0.0
    semidiameter = 9.0
    r = Ray([5.0, 5.0, 1.0], [0.0, 0.0, -1.0])

    temp = leaf(Sphere(abs(frontradius)), Optics.translation(0.0, 0.0, frontvertex - frontradius))
    d = -frontradius - sqrt(frontradius^2 - semidiameter^2) #offset from vertex to cutting plane. Plane is necessary to prevent parts of the complement from showing up in the final CSG
    p₀ = frontvertex + d
    plane = Plane(0.0, 0.0, 1.0, 0.0, 0.0, p₀)
    gen = csgdifference(leaf(plane), temp)

    concave = gen()
    intsct = surfaceintersection(concave, r)
    println(intsct)
end

# function testsinglepossurface()
#     backradius = 60.0
#     frontvertex = 0.0
#     thickness = 10.0
#     semidiameter = 9.0
#     r = Ray([5.000000042187876, 5.000000042187876, 0.6143579891640621], [0.042187876699796796, 0.042187876699796796, -0.9982185963600987])

#     temp = leaf(Sphere(abs(backradius)), Optics.translation(0.0, 0.0, (frontvertex - thickness) - backradius))
#     d = backradius - sqrt(backradius^2 - semidiameter^2)
#     p₀ = frontvertex - thickness - d
#     plane = Plane(0.0, 0.0, -1.0, 0.0, 0.0, p₀)
#     # gen = csgintersection(leaf(plane),csgcomplement(temp))
#     gen = csgcomplement(temp)
#     concave = gen()
#     intsct = closestintersection(surfaceintersection(concave, r))
#     println(intsct)
# end


# function testnegcurveintersection()
#     gen = Optics.csgcomplement(leaf(Sphere(60.0), Optics.translation(0.0, 0.0, 60.0)))
#     concave = gen()

#     r = Ray([0.0, 0.0, 1.0], [0.0, 0.0, -1.0])
#     intsct = surfaceintersection(concave, r)

#     println(intsct)
# end

function testdoubleconvex()
    system = Examples.doubleconvex()
    green = 500 * Unitful.u"nm"
    r = Ray([0.0, 0.0, 10.0], [0.0, 0.0, -1.0])

    res = Optics.trace(system, r, green)
    println(res)
end

function testcooketriplet()
    system = Examples.cooketriplet()

    green = 500 * Unitful.u"nm"
    temperature = 20 * Unitful.u"°C"

    r = Ray([0.0, 0.0, 1.0], [0.0, 0.0, -1.0])

    res = Optics.trace(system, r, green, temperature)
end

function testopticalsystem(lenssystem::T) where {T<:OpticalSystem}
    # lenssystem = Examples.cooketriplet()
    # lenssystem = Examples.planoconvex()
    lenssystem = Examples.doubleconvex()
    green = 500 * Unitful.u"nm"

    rad = Optics.semidiameter(lenssystem)
    points = Array{SVector{3,Float64},1}(undef, 0)

    count = 0

    for x in (-rad):0.5:rad, y in (-rad):0.5:rad
        r = Ray([x, y, 1.0], [0.0, 0.0, -1.0])
        if mod(count, 100) == 0
            println(count)
        end


        count += 1
        res = Optics.trace(lenssystem, r, green)

        if !(nothing === res)
            push!(points, point(res.intersection))
        end
    end

    twodpoints = map(x -> (x[1], x[2]), points)

    # scatter(twodpoints, markersize = 3)
end

function testelementsurfaceintersection()
    elt = Optics.SphericalLens(Optics.GlassCat.SCHOTT.BK6, 0.0, 60.0, Inf64, 5.0, 9.0)
    r = Ray([0.0, 0.0, 10.0], [0.0, 0.0, -1.0])
    res = surfaceintersection(elt, r)
    println(res)
end

function testsurfaceintersection()
    # hit from inside
    surf = TestData.beziersurface()
    accelsurf = AcceleratedParametricSurface(surf)
    TOLERANCE = 1e-9

    r = Ray([0.5, 0.5, 0.2], [0.0, 1.0, 0.0])
    println("1")
    res = surfaceintersection(accelsurf, r)
    lower(res) isa RayOrigin && isapprox(point(upper(res)), [0.5, 0.8996900152986361, 0.2], atol = TOLERANCE)

    # hit from outside
    r = Ray([0.0, 0.5, 0.2], [1.0, 0.0, -1.0])
    println("2")
    res = surfaceintersection(accelsurf, r)
    isapprox(point(lower(res)), [0.06398711204047353, 0.5, 0.13601288795952649], atol = TOLERANCE) && upper(res) isa Infinity

    # two hits from outside
    r = Ray([0.0, 0.5, 0.25], [1.0, 0.0, 0.0])
    println("3")
    res = surfaceintersection(accelsurf, r)
    isapprox(point(lower(res)), [0.12607777053030267, 0.5, 0.25], atol = TOLERANCE) && isapprox(point(upper(res)), [0.8705887077060419, 0.5, 0.25], atol = TOLERANCE)

    surf = TestData.wavybeziersurface()
    accelsurf = AcceleratedParametricSurface(surf)

    # 3 hit starting outside
    r = Ray([0.5, 0.0, 0.0], [0.0, 1.0, 0.0])
    println("4")
    res = surfaceintersection(accelsurf, r)
    isa(res, DisjointUnion) && length(res) == 2 && isapprox(point(lower(res[1])), [0.5, 0.10013694786182059, 0.0], atol = TOLERANCE) && isapprox(point(upper(res[1])), [0.5, 0.49625, 0.0], atol = TOLERANCE) && isapprox(point(lower(res[2])), [0.5, 0.8971357794109067, 0.0], atol = TOLERANCE) && (upper(res[2]) isa Infinity)

    # 3 hit starting inside
    #     r = Ray([0.5, 1.0, 0.0], [0.0, -1.0, 0.0])
    #     println("5")
    #     res = surfaceintersection(accelsurf, r)
    #     isa(res, DisjointUnion) && length(res) == 2 && (lower(res[1]) isa RayOrigin) && isapprox(point(upper(res[1])), [0.5, 0.8971357794109067, 0.0], atol = TOLERANCE) && isapprox(point(lower(res[2])), [0.5, 0.49625, 0.0], atol = TOLERANCE) && isapprox(point(upper(res[2])), [0.5, 0.10013694786182059, 0.0], atol = TOLERANCE)
    # println("after 5")

    surf = TestData.verywavybeziersurface()
    accelsurf = AcceleratedParametricSurface(surf)

    # five hits starting outside
    r = Ray([0.9, 0.0, -0.3], [0.0, 1.0, 0.7])
    println("6")
    res = surfaceintersection(accelsurf, r)
    a = isapprox(point(lower(res[1])), [0.9, 0.03172286522032046, -0.2777939943457758], atol = TOLERANCE) && isapprox(point(upper(res[1])), [0.9, 0.1733979947040411, -0.17862140370717122], atol = TOLERANCE)
    b = isapprox(point(lower(res[2])), [0.9, 0.5335760974594397, 0.07350326822160776], atol = TOLERANCE) && isapprox(point(upper(res[2])), [0.9, 0.7767707607392784, 0.24373953251749486], atol = TOLERANCE)
    c = isapprox(point(lower(res[3])), [0.9, 0.9830891958374246, 0.3881624370861975], atol = TOLERANCE) && (upper(res[3]) isa Infinity)
    isa(res, DisjointUnion) && length(res) == 3 && a && b && c

    # five hits starting inside
    r = Ray([0.9, 1.0, 0.4], [0.0, -1.0, -0.7])
    println("7")
    res = surfaceintersection(accelsurf, r)
    a = (lower(res[1]) isa RayOrigin) && isapprox(point(upper(res[1])), [0.9, 0.9830891958374246, 0.3881624370861975], atol = TOLERANCE)
    b = isapprox(point(lower(res[2])), [0.9, 0.7767707607392784, 0.24373953251749486], atol = TOLERANCE) && isapprox(point(upper(res[2])), [0.9, 0.5335760974594397, 0.07350326822160776], atol = TOLERANCE)
    c = isapprox(point(lower(res[3])), [0.9, 0.17339799470404108, -0.1786214037071712], atol = TOLERANCE) && isapprox(point(upper(res[3])), [0.9, 0.03172286522032046, -0.27779399434577573], atol = TOLERANCE)
    isa(res, DisjointUnion) && length(res) == 3 && a && b && c

    # 4 hits starting inside
    r = Ray([0.1, 0.0, -0.3], [0.0, 1.0, 0.7])
    println("8")
    res = surfaceintersection(accelsurf, r)
    a = (lower(res[1]) isa RayOrigin) && isapprox(point(upper(res[1])), [0.1, 0.2851860296285551, -0.10036977926001144], atol = TOLERANCE)
    b = isapprox(point(lower(res[2])), [0.1, 0.5166793625025807, 0.06167555375180668], atol = TOLERANCE) && isapprox(point(upper(res[2])), [0.1, 0.7770862508789854, 0.24396037561528983], atol = TOLERANCE)
    c = isapprox(point(lower(res[3])), [0.1, 0.98308919558696, 0.3881624369108719], atol = TOLERANCE) && (upper(res[3]) isa Infinity)
    isa(res, DisjointUnion) && length(res) == 3 && a && b && c

    # 4 hits starting outside
    r = Ray([0.9, 0.9, 0.4], [0.0, -0.9, -0.7])
    println("9")
    res = surfaceintersection(accelsurf, r)
    a = isapprox(point(lower(res[1])), [0.9, 0.736072142615238, 0.2725005553674076], atol = TOLERANCE) && isapprox(point(upper(res[1])), [0.9, 0.567439326091764, 0.141341698071372], atol = TOLERANCE)
    b = isapprox(point(lower(res[2])), [0.9, 0.16601081959179267, -0.1708804736508277], atol = TOLERANCE) && isapprox(point(upper(res[2])), [0.9, 0.032434058775915924, -0.274773509840954], atol = TOLERANCE)
    isa(res, DisjointUnion) && length(res) == 2 && a && b

end

function errortest()
    surf = TestData.verywavybeziersurface()
    accelsurf = AcceleratedParametricSurface(surf)

    TOLERANCE = 1e-9
    r = Ray([0.5, 1.0, 0.0], [0.0, -1.0, 0.0])

    res = surfaceintersection(accelsurf, r)
    isa(res, DisjointUnion) && length(res) == 2 && (lower(res[1]) isa RayOrigin) && isapprox(point(upper(res[1])), [0.5, 0.8971357794109067, 0.0], atol = TOLERANCE) && isapprox(point(lower(res[2])), [0.5, 0.49625, 0.0], atol = TOLERANCE) && isapprox(point(upper(res[2])), [0.5, 0.10013694786182059, 0.0], atol = TOLERANCE)
end

function testintervalunion()
    alphasequal(int1, int2) = α(lower(int1)) == α(lower(int2)) && α(upper(int1)) == α(upper(int2))
    intersectionat(a) = TestData.intersectionat(a)
    a = Interval(intersectionat(0.0), intersectionat(0.2))
    aa = Interval(intersectionat(0.1), intersectionat(0.3))

    b = DisjointUnion(aa, Interval(intersectionat(0.5), intersectionat(0.6)))
    res = intervalunion(a, b)

    @assert res == intervalunion(b, a)
    println("alphas equal $(alphasequal(res[1],Interval(intersectionat(0.0), intersectionat(0.3)))) ")


    # @assert res[1] == Interval(intersectionat(0.0), intersectionat(0.3)) this equality test does not work
    @assert res[2] == b[2]
end



function fixequalitybugintests()
    #conclusion: putting BezierSurface as a field in Intersection causes simple value semantics == to fail. Maybe BezierSurface field becomes a ref and is heap allocated?
    a = DisjointUnion(Interval(TestData.intersectionat(0.2), TestData.intersectionat(0.5)), Interval(TestData.intersectionat(0.7), TestData.intersectionat(0.8)))
    b = DisjointUnion(Interval(TestData.intersectionat(0.1), TestData.intersectionat(0.3)), Interval(TestData.intersectionat(0.4), TestData.intersectionat(0.6)))
    res = intervalunion(a, b)

    # println("a∪b $(intervalunion(a,b))")
    # println("identity $(intervalunion(a,b) == intervalunion(a,b))")
    # println("b∪a $(intervalunion(b,a))")
    @assert res == intervalunion(b, a)
    temp = Interval(TestData.intersectionat(0.1), TestData.intersectionat(0.6))
    println(temp)
    println(TestData.intersectionat(0.1))
    @assert res[1] == temp

    @assert res[2] == a[2]

end

function testcylinder()
    sph = Sphere(1.0)
    cyl = Cylinder(0.2)

    csg = csgdifference(leaf(sph), leaf(cyl))
    Vis.draw(csg)
end

function testtransmission()
    lens = TestData.planoplano()
    green = 500 * Unitful.u"nm"

    # result = trace(lens,r,green)
    # println(result)

    r = Ray([10.0, 0.0, 0.0], [-1.0, 0.0, -1.0])
    result = trace(lens, r, green)
    println(result)

    lenselement = lens.system.assembly.elements[1].objecttree
    println("Direct result of surface intersection on csg tree $(Optics.evalcsg(lenselement,r))")
    # pln1 = Plane(0.0,0.0,1.0,0.0,0.0,0.0)
    # pln2 = Plane(0.0,0.0,-1.0,0.0,0.0,-10.0)

    pln1 = Rectangle([0.0, 0.0, 1.0], [0.0, 0.0, 0.0], 9.0, 9.0)
    pln2 = Rectangle([0.0, 0.0, -1.0], [0.0, 0.0, -10.0], 9.0, 9.0)

    println("intersection top rectangle $(surfaceintersection(pln1,r)) \n\n\n intersection bottom rectangle $(surfaceintersection(pln2,r))")

    cyl = Cylinder(9.0, 10.0)
    csg = csgintersection(leaf(pln1), csgintersection(leaf(pln2), leaf(cyl)))
    intsct = Optics.evalcsg(csg(), r)
    println("rectangle top and bottom CSG intsct $intsct")
    # result = surfaceintersection(cyl, r)
    # println(result)
end

function testgenerateray()
    gen = Optics.RayGeneratorUniform(10, 2.0, 10.0)

    for i in gen
        println(i)
    end
end

function testrectangle()
    r = Rectangle([0.0, 0.0, 1.0], [0.0, 0.0, 0.0], 1.0, 1.0)
    ry = Ray([2.0, 0.0, 1.0], [0.0, 0.0, -1.0])
    println(surfaceintersection(r, ry))

end

function randunit()
    let v = rand(SVector{3,Float64})
        while (norm(v) > 1.0)
            v = rand(SVector{3,Float64})
        end
        return normalize(v)
    end
end

function makeeye()
    eye = CSGOpticalSystem(Optics.ModelEye(), Rectangle(10.0, 10.0, SVector{3,Float64}(0.0, 0.0, 1.0), SVector{3,Float64}(0.0, 0.0, -26.0)))
    # Vis.drawtracerays(eye)
    ray = OpticalRay([0.0, 0.0, 10.0], [0.0, 0.001, -1.0], 1.0, 0.55)

    return ray, eye
end

function dorays(ray, eye)
    # Threads.@threads for i in 1:1000000
    for i in 1:1000
        trace(eye, ray)
    end
end

function testnearestsquareroots()
    for i in 1:10:1000
        a, b = Optics.nearestsqrts(i)
        approx = a * b
        println(approx - i)
    end
end

function hierarchicalimage()
    a = Optics.HierarchicalImage{Float32}(100, 200)

    for i in CartesianIndices(a)
        a[i] = 2.0
    end
    return a
end



struct SpecialNumber{T<:Real} <: Real
    val::T

    function SpecialNumber(num::T) where {T<:Real}
        if isnan(num)
            println("here")
        else
            return new{T}(num)
        end
    end
end

# const operators = (:+,:-,:*,:/,:sqrt,:cos,:sin,:exp,:ln)
for op in (:sin, :cos, :tan, :log, :exp, :sqrt, :-)
    eval(quote
        Base.$op(a::SpecialNumber) = SpecialNumber($op(a.val))
    end)
end

for op in (:+, :-, :*, :/)
    eval(quote
        Base.$op(a::SpecialNumber, b::SpecialNumber) = SpecialNumber($op(a.val, b.val))
    end)
end

Base.promote_rule(::Type{T}, ::Type{SpecialNumber{T1}}) where {T<:Real,T1<:Real} = SpecialNumber{promote_type(T, T1)}




function singleemitter(numrays)
    T = Float64
    temp = SVector{3,T}(0.0, 0.0, 0.0)
    # λ = T(0)
    a = Source(UniformSpectrum{T}(), PointOrigin{T}(SVector{3,T}(0.0, 0.0, 0.0)), ConeDistribution(π / 5, numrays), Lambertian{T}())

    # function tracerays()
    #   temp = T(0)
    #         for ray in a
    #             temp = ray.power
    #         end
    #         return temp
    # end

    for ray in a
        temp = Optics.direction(a, 1)
        #  temp  = Optics.spectrumsample(a)
        # temp = Optics.spectrumpower(a.spectrum,.6)
    end

    return temp
end
export singleemitter

printrays(a::Source) =
    for ray in a
        println(ray)
    end

function diagdistributions(numsamples)
    temp = 0.0

    for i in 1:numsamples
        temp = rand(Distributions.Uniform{Float64}(400.0, 680.0))
    end
    return temp
end
export diagdistributions


end # module Junk
