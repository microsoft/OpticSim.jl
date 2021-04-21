# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

using Test
using OpticSim.GlassCat

using Unitful
using Base.Filesystem
using DataFrames
using StaticArrays

@testset "GlassCat" begin
    @testset "sources.jl" begin
        @testset "add_agf" begin
            tmpdir = mktempdir()
            agfdir = mktempdir(tmpdir)
            sourcefile, _ = mktemp(tmpdir)

            @testset "bad implicit catalog names" begin
                for agffile in split("a 1.agf a.Agf")
                    @test_logs (:error, "invalid implicit catalog name \"$agffile\". Should be purely alphabetical with a .agf/.AGF extension.") add_agf(agffile; agfdir, sourcefile)
                end
            end

            @testset "file not found" begin
                @test_logs (:error, "file not found at nonexistentfile.agf") add_agf("nonexistentfile.agf"; agfdir, sourcefile)
            end

            @testset "add to source file" begin
                for name in split("a b")
                    open(joinpath(tmpdir, "$name.agf"), "w") do io
                        write(io, "")
                    end
                end
                empty_sha = "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"

                @test isempty(readlines(sourcefile))

                add_agf(joinpath(tmpdir, "a.agf"); agfdir, sourcefile, rebuild=false)
                @test length(readlines(sourcefile)) === 1
                @test readlines(sourcefile)[1] === "a $empty_sha"

                add_agf(joinpath(tmpdir, "b.agf"); agfdir, sourcefile, rebuild=false)
                @test length(readlines(sourcefile)) === 2
                @test readlines(sourcefile)[1] === "a $empty_sha"
                @test readlines(sourcefile)[2] === "b $empty_sha"

                @test_logs (:error, "adding the catalog name \"a\" would create a duplicate entry in source file $sourcefile") add_agf(joinpath(tmpdir, "a.agf"); agfdir, sourcefile)
            end

            # TODO rebuild=true
        end

        # integration test
        @testset "verify_sources!" begin
            tmpdir = mktempdir()
            agfdir = mktempdir(tmpdir)

            sources = split.(readlines(GlassCat.SOURCES_PATH))
            GlassCat.verify_sources!(sources, agfdir)

            @test first.(sources) == first.(split.(readlines(GlassCat.SOURCES_PATH)))

            # TODO missing_sources
        end

        # TODO unit tests
    end

    @testset "generate.jl" begin
        CATALOG_NAME = "TEST_CAT"
        SOURCE_DIR = joinpath(@__DIR__, "..", "..", "test")
        SOURCE_FILE = joinpath(SOURCE_DIR, "$(CATALOG_NAME).agf")

        TEST_CAT_VALUES = DataFrame(
            name = split("MG463_ MG523 _5MG448 MG452_STAR"),
            raw_name = split("MG463& MG523 5MG448 MG452*"),
            dispform = repeat([2], 4),
            Nd = repeat([1.0], 4),
            Vd = repeat([0.0], 4),
            exclude_sub = [5, 0, 0, 0],
            status = [6, 0, 0, 0],
            meltfreq = repeat([-1], 4),
            TCE = [19.8, 18.3, 30.2, 25.7],
            p = [4.63, 5.23, 4.48, 4.52],
            ΔPgF = repeat([0.0], 4),
            ignore_thermal_exp = [3, 0, 0, 0],
            C1 = [6.781409960E+00, 4.401393677E+00, 8.427948454E+00, 4.656568861E+00],
            C2 = [9.868538020E-02, 2.234722960E-02, 6.325268390E-02, 1.170534894E-01],
            C3 = [3.469529960E-03, 4.562935070E+00, -2.441576539E+00, 1.409322578E+00],
            C4 = [2.277823700E+00, 3.239318260E-01, 1.997453950E-03, 2.445302672E-04],
            C5 = [2.126055280E+00, 6.863473749E-01, 4.975938104E-01, 6.352738622E-01],
            C6 = [3.227490100E+03, 1.625888344E+03, 1.443307391E+03, 1.397305161E+03],
            C7 = repeat([0.0], 4),
            C8 = repeat([0.0], 4),
            C9 = repeat([0.0], 4),
            C10 = repeat([0.0], 4),
            D₀ = [2.226000000E-05, 1.122401430E-04, -1.718689957E-05, -2.545621409E-07],
            D₁ = [5.150100000E-08, 2.868995096E-13, -6.499841076E-13, 1.258087492E-10],
            D₂ = [1.309700000E-10, -4.157332734E-16, -4.153311088E-16, 1.119703940E-11],
            E₀ = [4.382400000E-05, 9.532969159E-05, 8.220254442E-06, 2.184273853E-05],
            E₁ = [4.966000000E-07, -1.045247704E-11, 6.086168578E-12, -2.788674395E-09],
            λₜₖ = [3.728000000E-01, 6.338475915E-01, 3.610058073E-01, -1.948154058E+00],
            temp = [2.002000000E+01, 2.000000000E+01, 2.000000000E+01, 2.000000000E+01],
            relcost = [1.1, -1.0, -1.0, -1.0],
            CR = [1.2, -1.0, -1.0, -1.0],
            FR = [1.3, -1.0, -1.0, -1.0],
            SR = [1.4, -1.0, -1.0, -1.0],
            AR = [1.5, -1.0, -1.0, -1.0],
            PR = [1.6, -1.0, -1.0, -1.0],
            λmin = repeat([2.0], 4),
            λmax = repeat([14.0], 4),
            transmission = [
                [[2.0, 1.0, 3.0], [14.0, 4.0, 5.0]],
                [[2.0, 1.0, 1.0], [14.0, 1.0, 1.0]],
                [[2.0, 1.0, 1.0], [14.0, 1.0, 1.0]],
                [[2.0, 1.0, 1.0], [14.0, 1.0, 1.0]],
            ],
        )
        FIELDS = names(TEST_CAT_VALUES)[2:end]

        @testset "sourcefile_to_catalog" begin
            catalog = GlassCat.sourcefile_to_catalog(SOURCE_FILE)

            for row in eachrow(TEST_CAT_VALUES)
                name = row["name"]
                @test name ∈ keys(catalog)
                for field in FIELDS
                    @test row[field] == catalog[name][field]
                end
            end
        end

        # integration test
        @testset "generate_jls" begin
            jldir = mktempdir()
            mainfile = joinpath(jldir, "AGF_TEST_CAT.jl")

            GlassCat.generate_jls([CATALOG_NAME], mainfile, jldir, SOURCE_DIR)

            include(mainfile)

            for row in eachrow(TEST_CAT_VALUES)
                name = row["name"]
                glass = getfield(getfield(Main, Symbol(CATALOG_NAME)), Symbol(name))
                @test "$CATALOG_NAME.$name" ∈ AGF_GLASS_NAMES
                @test Symbol("$CATALOG_NAME.$name") ∈ AGF_GLASSES
                for field in FIELDS
                    if field === "raw_name"
                    elseif field === "transmission"
                        @test row[field] == getfield(glass, Symbol(field))[1:length(row[field])]
                        @test zero(SVector{100-length(row[field]),SVector{3,Float64}}) == getfield(glass, Symbol(field))[length(row[field])+1:end]
                    else
                        @test row[field] == getfield(glass, Symbol(field))
                    end
                end
            end

            # these used to be in the "Glass Tests" testset, but they rely on the generated AGF_TEST_CAT.jl file
            g = TEST_CAT.MG523
            @test GlassCat.index(g, ((g.λmin + g.λmax) / 2)u"μm") ≈ 3.1560980389455593 atol=1e-14
            @test_throws ErrorException GlassCat.index(g, (g.λmin - 1)u"μm")
            @test_throws ErrorException GlassCat.index(g, (g.λmax + 1)u"μm")

            # [hacky] AGF_TEST_CAT.jl redefines consts AGF_GLASS_NAMES and AGF_GLASSES - this restores them
            include(GlassCat.AGFGLASSCAT_PATH)
        end
    end

    @testset "Glass Tests" begin
        @test GlassCat.absairindex(500u"nm") ≈ 1.0002741948670688 atol=1e-14
        @test GlassCat.absairindex(600u"nm", temperature=35u"°C", pressure=2.0) ≈ 1.0005179096900811 atol=1e-14
        @test GlassCat.index(GlassCat.Air, 500u"nm") == 1.0
        @test GlassCat.index(GlassCat.Air, 600u"nm") == 1.0

        # test against true values
        g = GlassCat.SCHOTT.N_BK7
        @test GlassCat.index(g, 533u"nm") ≈ 1.519417351519283 atol=1e-14
        @test GlassCat.index(g, 533u"nm", temperature=35u"°C") ≈ 1.519462486258311 atol=1e-14
        @test GlassCat.index(g, 533u"nm", pressure=2.0) ≈ 1.518994119690216 atol=1e-14
        @test GlassCat.index(g, 533u"nm", temperature=35u"°C", pressure=2.0) ≈ 1.519059871499476 atol=1e-14

        # test transmission
        @test GlassCat.absorption(GlassCat.Air, 500u"nm") == 0.0
        @test GlassCat.absorption(GlassCat.Air, 600u"nm") == 0.0
        # TODO these are currently taken from this package (i.e. regression tests), ideally we would get true values somehow
        @test GlassCat.absorption(g, 500u"nm") == 0.0002407228930225164
        @test GlassCat.absorption(g, 600u"nm") == 0.00022060722752440445
        @test GlassCat.absorption(g, 3000u"nm") == 0.04086604990127926
        @test GlassCat.absorption(g, 600u"nm", temperature = 35u"°C", pressure=2.0) == 0.00022075540719494738

        # test that everything is alloc-less
        @test (@allocated GlassCat.absorption(g, 600u"nm")) == 0
        @test (@allocated GlassCat.index(g, 533u"nm")) == 0
        @test (@allocated GlassCat.SCHOTT.N_BK7) == 0

        @test GlassCat.glassforid(GlassCat.glassid(GlassCat.SCHOTT.N_BK7)) == GlassCat.SCHOTT.N_BK7

        @test GlassCat.isair(GlassCat.Air) == true
        @test GlassCat.isair(GlassCat.SCHOTT.N_BK7) == false

        # test MIL and model glasses
        fit_acc = 0.001
        bk7 = GlassCat.glassfromMIL(517642)
        @test GlassCat.glassforid(GlassCat.glassid(bk7)) == bk7
        @test GlassCat.index(bk7, 0.5875618) ≈ 1.517 atol=fit_acc
        @test GlassCat.index(bk7, 0.533) ≈ 1.519417351519283 atol=fit_acc
        @test GlassCat.index(bk7, 0.743) ≈ 1.511997032563557 atol=fit_acc
        @test GlassCat.index(bk7, 533u"nm", temperature = 35u"°C", pressure = 2.0) ≈ 1.519059871499476 atol=fit_acc

        t2 = GlassCat.glassfromMIL(1.135635)
        @test GlassCat.index(t2, 0.5875618) ≈ 2.135 atol=fit_acc

        fit_acc = 0.0001
        bk7 = GlassCat.modelglass(1.5168, 64.167336, -0.0009)
        @test GlassCat.glassforid(GlassCat.glassid(bk7)) == bk7
        @test GlassCat.index(bk7, 0.5875618) ≈ 1.5168 atol=fit_acc
        @test GlassCat.index(bk7, 0.533) ≈ 1.519417351519283 atol=fit_acc
        @test GlassCat.index(bk7, 0.743) ≈ 1.511997032563557 atol=fit_acc
        @test GlassCat.index(bk7, 533u"nm", temperature = 35u"°C", pressure = 2.0) ≈ 1.519059871499476 atol=fit_acc

        # test other glass
        @test GlassCat.index(GlassCat.CARGILLE.OG0608, 0.578) == 1.4596475735607324

        # make sure that the other functions work
        GlassCat.plot_indices(GlassCat.SCHOTT.N_BK7; polyfit = true, fiterror = true)
    end
end # testset GlassCat
