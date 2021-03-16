
using Test
using Unitful
using StaticArrays
import Unitful: Units

# ensure there aren't any ambiguities
@test isempty(detect_ambiguities(GlassCat))
@test isempty(detect_unbound_args(GlassCat))

@test_throws ErrorException Opticks.GlassCat.parse_glass_file!(Dict(), "NOFILE")

cat = Dict()

@testset "Parsing Tests" begin
    Opticks.GlassCat.parse_glass_file!(cat, "./TEST_CAT.AGF")
    # test that all glasses were parsed correctly
    @test "TEST_CAT" ∈ keys(cat)
    @test "MG463_" ∈ keys(cat["TEST_CAT"])
    @test "MG523" ∈ keys(cat["TEST_CAT"])
    @test "_5MG448" ∈ keys(cat["TEST_CAT"])
    @test "MG452_STAR" ∈ keys(cat["TEST_CAT"])

    # check that parsing all the variables worked for one example
    glass = cat["TEST_CAT"]["MG463_"]

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

    @test Opticks.GlassCat.load_glass_db(".") == cat
end

@testset "Module Gen Tests" begin
    Opticks.GlassCat.generate_cat_jl(cat, "./TEST_CAT.jl")
    include("./TEST_CAT.jl")

    @test :MG463_ ∈ names(TEST_CAT, all = true)
    @test :MG523 ∈ names(TEST_CAT, all = true)
    @test :_5MG448 ∈ names(TEST_CAT, all = true)
    @test :MG452_STAR ∈ names(TEST_CAT, all = true)

    # test that the struct matches the dict
    glassd = cat["TEST_CAT"]["MG452_STAR"]
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

    rm("./TEST_CAT.jl")
end

@testset "Glass Tests" begin
    g = TEST_CAT.MG523
    @test isapprox(index(g, ((g.λmin + g.λmax) / 2) * u"μm"), 3.1560980389455593, atol = 1e-14)
    @test_throws ErrorException index(g, (g.λmin - 1) * u"μm")
    @test_throws ErrorException index(g, (g.λmax + 1) * u"μm")
    @test isapprox(absairindex(500 * u"nm"), 1.0002741948670688, atol = 1e-14)
    @test isapprox(absairindex(600 * u"nm", temperature = 35 * u"°C", pressure = 2.0), 1.0005179096900811, atol = 1e-14)
    @test index(Opticks.GlassCat.Air, 500 * u"nm") == 1.0
    @test index(Opticks.GlassCat.Air, 600 * u"nm") == 1.0
    @test isapprox(index(Opticks.GlassCat.Vacuum, 500 * u"nm"), 0.9997258841958212, atol = 1e-14)

    # test against true values
    g = Opticks.GlassCat.SCHOTT.N_BK7
    @test isapprox(index(g, 533 * u"nm"), 1.519417351519283, atol = 1e-14)
    @test isapprox(index(g, 533 * u"nm", temperature = 35 * u"°C"), 1.519462486258311, atol = 1e-14)
    @test isapprox(index(g, 533 * u"nm", pressure = 2.0), 1.518994119690216, atol = 1e-14)
    @test isapprox(index(g, 533 * u"nm", temperature = 35 * u"°C", pressure = 2.0), 1.519059871499476, atol = 1e-14)

    # test transmission
    @test absorption(Opticks.GlassCat.Air, 500 * u"nm") == 0.0
    @test absorption(Opticks.GlassCat.Air, 600 * u"nm") == 0.0
    # TODO these are currently taken from this package (i.e. regression tests), ideally we would get true values somehow
    @test absorption(g, 500 * u"nm") == 0.00016032085590155288
    @test absorption(g, 600 * u"nm") == 0.00015028195510849153
    @test absorption(g, 3000 * u"nm") == 0.022484756726141652
    @test absorption(g, 600 * u"nm", temperature = 35 * u"°C", pressure = 2.0) == 0.0001502079956238452

    # test that everything is alloc-less
    @test (@allocated absorption(g, 600 * u"nm")) == 0
    @test (@allocated index(g, 533 * u"nm")) == 0
    @test (@allocated Opticks.GlassCat.SCHOTT.N_BK7) == 0

    @test glassforid(glassid(Opticks.GlassCat.SCHOTT.N_BK7)) == Opticks.GlassCat.SCHOTT.N_BK7

    @test isair(Opticks.GlassCat.Air) == true
    @test isair(Opticks.GlassCat.SCHOTT.N_BK7) == false

    # test MIL and model glasses
    fit_acc = 0.001
    bk7 = Opticks.GlassCat.glassfromMIL(517642)
    @test glassforid(glassid(bk7)) == bk7
    @test isapprox(index(bk7, 0.5875618), 1.517, atol = fit_acc)
    @test isapprox(index(bk7, 0.533), 1.519417351519283, atol = fit_acc)
    @test isapprox(index(bk7, 0.743), 1.511997032563557, atol = fit_acc)
    @test isapprox(index(bk7, 533 * u"nm", temperature = 35 * u"°C", pressure = 2.0), 1.519059871499476, atol = fit_acc)

    t2 = Opticks.GlassCat.glassfromMIL(1.135635)
    @test isapprox(index(t2, 0.5875618), 2.135, atol = fit_acc)

    fit_acc = 0.0001
    bk7 = Opticks.GlassCat.modelglass(1.5168, 64.167336, -0.0009)
    @test glassforid(glassid(bk7)) == bk7
    @test isapprox(index(bk7, 0.5875618), 1.5168, atol = fit_acc)
    @test isapprox(index(bk7, 0.533), 1.519417351519283, atol = fit_acc)
    @test isapprox(index(bk7, 0.743), 1.511997032563557, atol = fit_acc)
    @test isapprox(index(bk7, 533 * u"nm", temperature = 35 * u"°C", pressure = 2.0), 1.519059871499476, atol = fit_acc)

    # test other glass
    @test index(Opticks.GlassCat.CARGILLE.OG0608, 0.578) == 1.4596475735607324

    # make sure that the other functions work
    plot_indices(Opticks.GlassCat.SCHOTT.N_BK7; polyfit = true, fiterror = true)
end
