# MIT License

# Copyright (c) Microsoft Corporation.

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE

@testset "GlassCat" begin
    # ensure there aren't any ambiguities
    @test isempty(detect_ambiguities(GlassCat))
    @test isempty(detect_unbound_args(GlassCat))

    include(joinpath("..", "..", "src", "GlassCat", "GlassTypes.jl"))
    include(joinpath(@__DIR__, "..", "..", "deps", "generate.jl"))

    CATALOG_NAME = "TEST_CAT"
    SOURCE_DIR = joinpath(@__DIR__, "..")
    SOURCE_FILE = joinpath(SOURCE_DIR, "$(CATALOG_NAME).agf")
    tmpdir = mktempdir()
    MAIN_FILE = joinpath(tmpdir, "AGF_TEST_CAT.jl")
    JL_DIR = tmpdir

    cat = Dict()

    @testset "Parsing Tests" begin
        cat = sourcefile_to_catalog(SOURCE_FILE)

        # test that all glasses were parsed correctly
        @test "MG463_" ∈ keys(cat)
        @test "MG523" ∈ keys(cat)
        @test "_5MG448" ∈ keys(cat)
        @test "MG452_STAR" ∈ keys(cat)

        # check that parsing all the variables worked for one example
        glass = cat["MG463_"]

        @test glass["raw_name"] == "MG463&"

        @test glass["dispform"] == 2
        @test glass["Nd"] == 1.0
        @test glass["Vd"] == 0.0
        @test glass["exclude_sub"] == 5
        @test glass["status"] == 6
        @test glass["meltfreq"] == -1

        @test glass["TCE"] == 19.8
        @test glass["p"] == 4.63
        @test glass["ΔPgF"] == 0.0
        @test glass["ignore_thermal_exp"] == 3

        @test glass["C1"] == 6.781409960
        @test glass["C2"] == 0.09868538020
        @test glass["C3"] == 0.003469529960
        @test glass["C4"] == 2.277823700
        @test glass["C5"] == 2.126055280
        @test glass["C6"] == 3227.490100
        @test glass["C7"] == 0.000000000
        @test glass["C8"] == 0.000000000
        @test glass["C9"] == 0.000000000
        @test glass["C10"] == 0.000000000

        @test glass["D₀"] == 2.226000000E-05
        @test glass["D₁"] == 5.150100000E-08
        @test glass["D₂"] == 1.309700000E-10
        @test glass["E₀"] == 4.382400000E-05
        @test glass["E₁"] == 4.966000000E-07
        @test glass["λₜₖ"] == 3.728000000E-01
        @test glass["temp"] == 2.002000000E+01

        @test glass["relcost"] == 1.1
        @test glass["CR"] == 1.2
        @test glass["FR"] == 1.3
        @test glass["SR"] == 1.4
        @test glass["AR"] == 1.5
        @test glass["PR"] == 1.6

        @test glass["λmin"] == 2.0
        @test glass["λmax"] == 14.0

        @test glass["transmission"] == [[2.0, 1.0, 3.0], [14.0, 4.0, 5.0]]
    end

    @testset "Module Gen Tests" begin
        generate_jls([CATALOG_NAME], MAIN_FILE, JL_DIR, SOURCE_DIR)
        include(MAIN_FILE)

        @test :MG463_ ∈ names(TEST_CAT, all = true)
        @test :MG523 ∈ names(TEST_CAT, all = true)
        @test :_5MG448 ∈ names(TEST_CAT, all = true)
        @test :MG452_STAR ∈ names(TEST_CAT, all = true)

        # test that the struct matches the dict
        glassd = cat["MG452_STAR"]
        glass = TEST_CAT.MG452_STAR

        @test glass.dispform == glassd["dispform"]

        @test glass.Nd == glassd["Nd"]
        @test glass.Vd == glassd["Vd"]
        @test glass.exclude_sub == glassd["exclude_sub"]
        @test glass.status == glassd["status"]
        @test glass.meltfreq == glassd["meltfreq"]

        @test glass.TCE == glassd["TCE"]
        @test glass.p == glassd["p"]
        @test glass.ΔPgF == glassd["ΔPgF"]
        @test glass.ignore_thermal_exp == glassd["ignore_thermal_exp"]

        @test glass.C1 == glassd["C1"]
        @test glass.C2 == glassd["C2"]
        @test glass.C3 == glassd["C3"]
        @test glass.C4 == glassd["C4"]
        @test glass.C5 == glassd["C5"]
        @test glass.C6 == glassd["C6"]
        @test glass.C7 == glassd["C7"]
        @test glass.C8 == glassd["C8"]
        @test glass.C9 == glassd["C9"]
        @test glass.C10 == glassd["C10"]

        @test glass.D₀ == glassd["D₀"]
        @test glass.D₁ == glassd["D₁"]
        @test glass.D₂ == glassd["D₂"]
        @test glass.E₀ == glassd["E₀"]
        @test glass.E₁ == glassd["E₁"]
        @test glass.λₜₖ == glassd["λₜₖ"]
        @test glass.temp == glassd["temp"]

        @test glass.relcost == glassd["relcost"]
        @test glass.CR == glassd["CR"]
        @test glass.FR == glassd["FR"]
        @test glass.SR == glassd["SR"]
        @test glass.AR == glassd["AR"]
        @test glass.PR == glassd["PR"]

        @test glass.λmin == glassd["λmin"]
        @test glass.λmax == glassd["λmax"]
    end

    @testset "Glass Tests" begin
        include(MAIN_FILE)

        g = TEST_CAT.MG523
        @test isapprox(OpticSim.GlassCat.index(g, ((g.λmin + g.λmax) / 2) * u"μm"), 3.1560980389455593, atol = 1e-14)
        @test_throws ErrorException OpticSim.GlassCat.index(g, (g.λmin - 1) * u"μm")
        @test_throws ErrorException OpticSim.GlassCat.index(g, (g.λmax + 1) * u"μm")
        @test isapprox(OpticSim.GlassCat.absairindex(500 * u"nm"), 1.0002741948670688, atol = 1e-14)
        @test isapprox(OpticSim.GlassCat.absairindex(600 * u"nm", temperature = 35 * u"°C", pressure = 2.0), 1.0005179096900811, atol = 1e-14)
        @test OpticSim.GlassCat.index(OpticSim.GlassCat.Air, 500 * u"nm") == 1.0
        @test OpticSim.GlassCat.index(OpticSim.GlassCat.Air, 600 * u"nm") == 1.0
        # @test isapprox(index(OpticSim.GlassCat.Vacuum, 500 * u"nm"), 0.9997258841958212, atol = 1e-14)

        # test against true values
        g = OpticSim.GlassCat.SCHOTT.N_BK7
        @test isapprox(OpticSim.GlassCat.index(g, 533 * u"nm"), 1.519417351519283, atol = 1e-14)
        @test isapprox(OpticSim.GlassCat.index(g, 533 * u"nm", temperature = 35 * u"°C"), 1.519462486258311, atol = 1e-14)
        @test isapprox(OpticSim.GlassCat.index(g, 533 * u"nm", pressure = 2.0), 1.518994119690216, atol = 1e-14)
        @test isapprox(OpticSim.GlassCat.index(g, 533 * u"nm", temperature = 35 * u"°C", pressure = 2.0), 1.519059871499476, atol = 1e-14)

        # test transmission
        @test OpticSim.GlassCat.absorption(OpticSim.GlassCat.Air, 500 * u"nm") == 0.0
        @test OpticSim.GlassCat.absorption(OpticSim.GlassCat.Air, 600 * u"nm") == 0.0
        # TODO these are currently taken from this package (i.e. regression tests), ideally we would get true values somehow
        @test OpticSim.GlassCat.absorption(g, 500 * u"nm") == 0.0002407228930225164
        @test OpticSim.GlassCat.absorption(g, 600 * u"nm") == 0.00022060722752440445
        @test OpticSim.GlassCat.absorption(g, 3000 * u"nm") == 0.04086604990127926
        @test OpticSim.GlassCat.absorption(g, 600 * u"nm", temperature = 35 * u"°C", pressure = 2.0) == 0.00022075540719494738

        # test that everything is alloc-less
        @test (@allocated OpticSim.GlassCat.absorption(g, 600 * u"nm")) == 0
        @test (@allocated OpticSim.GlassCat.index(g, 533 * u"nm")) == 0
        @test (@allocated OpticSim.GlassCat.SCHOTT.N_BK7) == 0

        @test OpticSim.GlassCat.glassforid(OpticSim.GlassCat.glassid(OpticSim.GlassCat.SCHOTT.N_BK7)) == OpticSim.GlassCat.SCHOTT.N_BK7

        @test OpticSim.GlassCat.isair(OpticSim.GlassCat.Air) == true
        @test OpticSim.GlassCat.isair(OpticSim.GlassCat.SCHOTT.N_BK7) == false

        # test MIL and model glasses
        fit_acc = 0.001
        bk7 = OpticSim.GlassCat.glassfromMIL(517642)
        @test OpticSim.GlassCat.glassforid(OpticSim.GlassCat.glassid(bk7)) == bk7
        @test OpticSim.GlassCat.isapprox(OpticSim.GlassCat.index(bk7, 0.5875618), 1.517, atol = fit_acc)
        @test OpticSim.GlassCat.isapprox(OpticSim.GlassCat.index(bk7, 0.533), 1.519417351519283, atol = fit_acc)
        @test OpticSim.GlassCat.isapprox(OpticSim.GlassCat.index(bk7, 0.743), 1.511997032563557, atol = fit_acc)
        @test OpticSim.GlassCat.isapprox(OpticSim.GlassCat.index(bk7, 533 * u"nm", temperature = 35 * u"°C", pressure = 2.0), 1.519059871499476, atol = fit_acc)

        t2 = OpticSim.GlassCat.glassfromMIL(1.135635)
        @test OpticSim.GlassCat.isapprox(OpticSim.GlassCat.index(t2, 0.5875618), 2.135, atol = fit_acc)

        fit_acc = 0.0001
        bk7 = OpticSim.GlassCat.modelglass(1.5168, 64.167336, -0.0009)
        @test OpticSim.GlassCat.glassforid(OpticSim.GlassCat.glassid(bk7)) == bk7
        @test isapprox(OpticSim.GlassCat.index(bk7, 0.5875618), 1.5168, atol = fit_acc)
        @test isapprox(OpticSim.GlassCat.index(bk7, 0.533), 1.519417351519283, atol = fit_acc)
        @test isapprox(OpticSim.GlassCat.index(bk7, 0.743), 1.511997032563557, atol = fit_acc)
        @test isapprox(OpticSim.GlassCat.index(bk7, 533 * u"nm", temperature = 35 * u"°C", pressure = 2.0), 1.519059871499476, atol = fit_acc)

        # test other glass
        @test OpticSim.GlassCat.index(OpticSim.GlassCat.CARGILLE.OG0608, 0.578) == 1.4596475735607324

        # make sure that the other functions work
        OpticSim.GlassCat.plot_indices(OpticSim.GlassCat.SCHOTT.N_BK7; polyfit = true, fiterror = true)
    end
end # testset GlassCat
