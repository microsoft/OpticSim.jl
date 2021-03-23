### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ d6c83680-879e-11eb-31d4-7dbda7e93e48
begin
	import PlutoUI

	import JSServe
	using Makie
	using AbstractPlotting
	using AbstractPlotting.MakieLayout
	using OpticSim
	using OpticSim.Emitters
	using OpticSim.Emitters.Geometry
	
	defs = NB.Defs("ran")

	NB.DefsClearHTML(defs)
	
	NB.DefsAddHTML(defs, 
		NB.HTMLFromObj(
			JSServe.Page()	
		)
	)
	
	# this function is needed to allow the visualization scene to be displayed inside a pluto notebook
	function Makie.display(obj)
		@info "RG: $obj"
		return obj
	end
	
	#text of this cell to appear in the notebook
	NB.DefsAddHTML(defs, NB.HTMLFromObj(md"Basic Initialization (code)"))

	PlutoUI.Show(MIME"text/html"(), NB.DefsHTML(defs))	
end

# ╔═╡ a891a210-87bc-11eb-311d-4d020986fe19
begin
	using DataFrames
	
	# defining the optical system
	sys_cooke = AxisymmetricOpticalSystem(
		DataFrame(Surface = [:Object, 1, 2, 3, :Stop, 5, 6, :Image],
				  Radius = [Inf, 26.777, 66.604, -35.571, 35.571, 35.571, -26.777, Inf],
				  Thickness = [Inf, 4.0, 2.0, 4.0, 2.0, 4.0, 44.748, missing],
				  Material = [OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_SK16, OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_SF2, OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_SK16, OpticSim.GlassCat.Air, missing],
				  SemiDiameter = [Inf, 8.580, 7.513, 7.054, 6.033, 7.003, 7.506, 15.0]))
	@show sys_cooke
end

# ╔═╡ 6ce01760-879e-11eb-2b13-4d3d07c4b4ce
md"# Various Examples"

# ╔═╡ 0d8de290-8ba1-11eb-0ae9-79fed06aae32
md"## Zoom Lenses"

# ╔═╡ e88e92c0-8ba3-11eb-11ba-d7e44bf5373e
md"## Focusing"

# ╔═╡ 576fb1c0-8ba3-11eb-1c10-fdb51172e2c9
md"## Lens Construction"

# ╔═╡ 0b6553c2-8860-11eb-0c5e-89632a8dd612
md"## Cooke Triplet Example"

# ╔═╡ 5f77cfd0-8854-11eb-377c-ef975f9abf63
md"## Initialization Stuff - Can be ignored"

# ╔═╡ 68e9b210-87ad-11eb-0f3a-5bb2dbf7d86c
begin
	NB.DefsClearHTML(defs)

	default_resolution_info = NB.GetVarInfo(
		@bind default_resolution PlutoUI.Select(["Small", "Medium", "Large"], default="Medium")
	)
	makie_backend_info = NB.GetVarInfo(
		@bind makie_backend PlutoUI.Select(["Static", "Web"], default="Static")
	)
	show_toc_info = NB.GetVarInfo(
		@bind show_toc PlutoUI.CheckBox(default=true)
	)

	fb3_items = [
		"Options",
		"@  Drawings Size $(default_resolution_info.html)",
		"Document Options",
		"@  $(show_toc_info.html) Show Table Of Content", 
		"@  Makie Backend $(makie_backend_info.html) <i><b>(Leave as Static for now)</b></i>",
		]	

	NB.DefsAddHTML(defs,
		NB.HTMLFloatingBox(
			fb3_items, 
			name="plutoui-docmenu", 
			header="Document Menu",
			width="25%",
			top="45%",
		)
	)

	# wide document layout
	NB.DefsAddHTML(defs, NB.HTMLNewDocLayout())

	# add TOC
	NB.DefsAddHTML(defs, 
		NB.HTMLFromObj(
			PlutoUI.TableOfContents(title = "Document Table of Content", depth = 4)
		)*
		NB.HTMLFixTOC()
	)
	
	
	#text of this cell to appear in the notebook
	NB.DefsAddHTML(defs, NB.HTMLFromObj(md"Define Floating Menu, TOC and New Layout (code)"))
	
	PlutoUI.Show(MIME"text/html"(), NB.DefsHTML(defs))	
end


# ╔═╡ 3a5d3ba0-87ae-11eb-1717-93be0b802cab
begin
	drawing = 1
	NB.SetBackend(defs, makie_backend)
	md"Set Makie backend"
end

# ╔═╡ 5db83180-8ba1-11eb-0364-ad9cfed99509
begin
	drawing	# notebook only - create dependency on drawing backend - comment if running in REPL

	# glass entrance lens on telescope
	topsurf = Plane(OpticSim.SVector(0.0, 0.0, 1.0), OpticSim.SVector(0.0, 0.0, 0.0), interface = FresnelInterface{Float64}(OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air), vishalfsizeu = 12.00075, vishalfsizev = 12.00075)
	botsurf = AcceleratedParametricSurface(ZernikeSurface(12.00075, radius = -1.14659768e+4, aspherics = [(4, 3.68090959e-7), (6, 2.73643352e-11), (8, 3.20036892e-14)]), 17, interface = FresnelInterface{Float64}(OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air))
	coverlens = csgintersection(leaf(Cylinder(12.00075, 1.4)), csgintersection(leaf(topsurf), leaf(botsurf, RigidBodyTransform(OpticSim.rotmatd(0, 180, 0), OpticSim.SVector(0.0, 0.0, -0.65)))))
	# big mirror with a hole in it
	bigmirror = ConicLens(OpticSim.GlassCat.SCHOTT.N_BK7, -72.65, -95.2773500000134, 0.077235, Inf, 0.0, 0.2, 12.18263, frontsurfacereflectance = 1.0)
	bigmirror = csgdifference(bigmirror, leaf(Cylinder(4.0, 0.3, interface = opaqueinterface()), OpticSim.translation(0.0, 0.0, -72.75)))
	# small mirror supported on a spider
	smallmirror = SphericalLens(OpticSim.GlassCat.SCHOTT.N_BK7, -40.65, Inf, -49.6845, 1.13365, 4.3223859, backsurfacereflectance = 1.0)
	obscuration1 = OpticSim.Circle(4.5, OpticSim.SVector(0.0, 0.0, 1.0), OpticSim.SVector(0.0, 0.0, -40.649), interface = opaqueinterface())
	obscurations2 = Spider(3, 0.5, 12.0, OpticSim.SVector(0.0, 0.0, -40.65))
	# put it together with the detector
	la = LensAssembly(coverlens(), bigmirror(), smallmirror(), obscuration1, obscurations2...)
	det = OpticSim.Circle(3.0, OpticSim.SVector(0.0, 0.0, 1.0), OpticSim.SVector(0.0, 0.0, -92.4542988), interface = opaqueinterface())
	tele = CSGOpticalSystem(la, det)

	Vis.drawtracerays(tele, raygenerator = UniformOpticalSource(CollimatedSource(GridRectOriginPoints(5, 5, 10.0, 10.0, position = OpticSim.SVector(0.0, 0.0, 20.0))), 0.55), trackallrays = true, colorbynhits = true, test = true)
end

# ╔═╡ f08a07c0-8ba3-11eb-382a-23dc079b623a
begin
	drawing	# notebook only - create dependency on drawing backend - comment if running in REPL

	local rect = Rectangle(5.0, 5.0, OpticSim.SVector(0.0, 0.0, 1.0), OpticSim.SVector(0.0, 0.0, 0.0))
	local int = HologramInterface(OpticSim.SVector(0.0, -3.0, -20.0), ConvergingBeam, OpticSim.SVector(0.0, 0.0, -1.0), CollimatedBeam, 0.55, 9.0, OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, OpticSim.GlassCat.Air, OpticSim.GlassCat.Air, 0.05, false)
	local obj = HologramSurface(rect, int)
	local sys = CSGOpticalSystem(LensAssembly(obj), Rectangle(10.0, 10.0, OpticSim.SVector(0.0, 0.0, 1.0), OpticSim.SVector(0.0, 0.0, -25.0), interface = opaqueinterface()))
	Vis.drawtracerays(sys; raygenerator = UniformOpticalSource(CollimatedSource(GridRectOriginPoints(5, 5, 3.0, 3.0, position = OpticSim.SVector(0.0, 0.0, 10.0), direction = OpticSim.SVector(0.0, 0.0, -1.0))), 0.55), trackallrays = true, rayfilter = nothing, test = true)
end

# ╔═╡ 5e2debd0-8ba3-11eb-0c4e-09d110230e2d
begin
	drawing	# notebook only - create dependency on drawing backend - comment if running in REPL

	topsurface = leaf(AcceleratedParametricSurface(QTypeSurface(9.0, radius = -25.0, conic = 0.3, αcoeffs = [(1, 0, 0.3), (1, 1, 1.0)], βcoeffs = [(1, 0, -0.1), (2, 0, 0.4), (3, 0, -0.6)], normradius = 9.5), interface = FresnelInterface{Float64}(OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air)), OpticSim.translation(0.0, 0.0, 5.0))
	botsurface = leaf(Plane(0.0, 0.0, -1.0, 0.0, 0.0, -5.0, vishalfsizeu = 9.5, vishalfsizev = 9.5, interface = FresnelInterface{Float64}(OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air)))
	barrel = leaf(Cylinder(9.0, 20.0, interface = FresnelInterface{Float64}(OpticSim.GlassCat.SCHOTT.N_BK7, OpticSim.GlassCat.Air, reflectance = zero(Float64), transmission = zero(Float64))))
	lens = csgintersection(barrel, csgintersection(topsurface, botsurface))(RigidBodyTransform{Float64}(0.0, Float64(π), 0.0, 0.0, 0.0, -5.0))
	sys = CSGOpticalSystem(LensAssembly(lens), Rectangle(15.0, 15.0, [0.0, 0.0, 1.0], [0.0, 0.0, -67.8], interface = opaqueinterface()))
	Vis.drawtracerays(sys, test = true, trackallrays = true)	
end

# ╔═╡ def239f0-87bc-11eb-2edb-2f859ac41bee
begin
	drawing	# notebook only - create dependency on drawing backend - comment if running in REPL
	
	local s1 = Sources.Source(sourcenum = 1, origins=Origins.Hexapolar(5, 8.0, 8.0), directions=Directions.Constant())
	
	# create a second source and rotate it by -5 degs
	local s2 = Sources.Source(sourcenum = 2, origins=Origins.Hexapolar(5, 8.0, 8.0), directions=Directions.Constant(), transform=Transform(origin(), rotYmat(deg2rad(-8)) * unitZ()))
	
	# create the "ray generator"
	local combined_sources = Sources.CompositeSource(Transform(Vec3D(0.0, 0.0, 10.0), unitZ() * -1), [s1 s2])
	
	# and draw the system + the generated rays
	Vis.drawtracerays(sys_cooke, raygenerator = combined_sources, test = true, trackallrays = true, colorbysourcenum = true, drawgen = false)
end

# ╔═╡ 96e423a0-885a-11eb-02a3-8704e8dbdab6
begin
	function Vis.scene(resolution = (1000, 1000))
		@info "RG: Vis.Scene Replacement"
		scene, layout = MakieLayout.layoutscene(resolution = resolution)
		Vis.set_current_main_scene(scene)
		lscene = layout[1, 1] = MakieLayout.LScene(scene, scenekw = (camera = Makie.cam3d_cad!, axis_type = Makie.axis3d!, raw = false))
		Vis.set_current_3d_scene(lscene)
		return scene, lscene
	end
	
	function resolution()
		if (default_resolution == "Small")
			return (300, 300)
		elseif (default_resolution == "Medium")
			return (500, 500)
		else
			return (1000, 1000)
		end
	end
end

# ╔═╡ Cell order:
# ╟─6ce01760-879e-11eb-2b13-4d3d07c4b4ce
# ╟─0d8de290-8ba1-11eb-0ae9-79fed06aae32
# ╠═5db83180-8ba1-11eb-0364-ad9cfed99509
# ╠═e88e92c0-8ba3-11eb-11ba-d7e44bf5373e
# ╠═f08a07c0-8ba3-11eb-382a-23dc079b623a
# ╠═576fb1c0-8ba3-11eb-1c10-fdb51172e2c9
# ╠═5e2debd0-8ba3-11eb-0c4e-09d110230e2d
# ╟─0b6553c2-8860-11eb-0c5e-89632a8dd612
# ╠═a891a210-87bc-11eb-311d-4d020986fe19
# ╠═def239f0-87bc-11eb-2edb-2f859ac41bee
# ╟─5f77cfd0-8854-11eb-377c-ef975f9abf63
# ╠═d6c83680-879e-11eb-31d4-7dbda7e93e48
# ╠═68e9b210-87ad-11eb-0f3a-5bb2dbf7d86c
# ╠═3a5d3ba0-87ae-11eb-1717-93be0b802cab
# ╠═96e423a0-885a-11eb-02a3-8704e8dbdab6
