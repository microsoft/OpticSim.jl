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

# ╔═╡ 67c8fac0-974e-11eb-2898-d14b79387954
begin
	init_notebook = true
	using OpticSim, OpticSim.Geometry, OpticSim.Emitters
	import OpticSim.NotebooksUtils as NB

	NB.InitNotebook()
end

# ╔═╡ a891a210-87bc-11eb-311d-4d020986fe19
begin
	using DataFrames
	
	# defining the optical system
	sys = AxisymmetricOpticalSystem(
		DataFrame(Surface = [:Object, 1, 2, 3, :Stop, 5, 6, :Image],
				  Radius = [Inf, 26.777, 66.604, -35.571, 35.571, 35.571, -26.777, Inf],
				  Thickness = [Inf, 4.0, 2.0, 4.0, 2.0, 4.0, 44.748, missing],
				  Material = [OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_SK16, OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_SF2, OpticSim.GlassCat.Air, OpticSim.GlassCat.SCHOTT.N_SK16, OpticSim.GlassCat.Air, missing],
				  SemiDiameter = [Inf, 8.580, 7.513, 7.054, 6.033, 7.003, 7.506, 15.0]))
	@show sys
end

# ╔═╡ d6c83680-879e-11eb-31d4-7dbda7e93e48
begin
	init_notebook
	
	import PlutoUI

	using Makie
	using AbstractPlotting
	using AbstractPlotting.MakieLayout
	
	defs = OpticSim.NotebooksUtils.Defs("ran")

	NB.DefsClearHTML(defs)
	
	# this function is needed to allow the visualization scene to be displayed inside a pluto notebook
	function Makie.display(obj)
		# @info "RG: $obj"
		return obj
	end
	
	#text of this cell to appear in the notebook
	NB.DefsAddHTML(defs, NB.HTMLFromObj(md"Basic Initialization (code)"))

	PlutoUI.Show(MIME"text/html"(), NB.DefsHTML(defs))	
end

# ╔═╡ 6ce01760-879e-11eb-2b13-4d3d07c4b4ce
md"# Introductions To Emitters"

# ╔═╡ 87954830-8854-11eb-2baf-d3b79b804e07
md"## Basic Emitters"

# ╔═╡ 9608ed90-8854-11eb-080d-4d1f60d596f5
md"""
Emitters are defined by Pixels and Spatial Layouts. An emitter has a spectrum, and an optical power distribution over the hemisphere. These are intrinsic physical properties of the emitter.

The **basic emitter** is constructed as a combination of 4 basic elements and a 3D transform. The basic elements include:
- Spectrum
- Angular Power Distribution
- Rays Origins Distribution
- Rays Directions Distribution

The OpticSim package comes with various implementations of each of these basic elements:
- Spectrum - the **generate** interface returns a tuple (power, wavelength)
  * **Uniform** - A flat spectrum bounded (default: from 450nm to 680nm). the range is sampled uniformly.
  * **DeltaFunction** - Constant wave length.
  * **Measured** - measured spectrum to compute emitter power and wavelength (created by reading CSV files – more details will follow).
- Angular Power Distribution - the interface **apply** returns an OpticalRay with modified power
  * **Lambertian**
  * **Cosine**
  * **Gaussian**
- Rays Origins Distribution - the interface **length** returns the number of samples, and **generate** returns the n'th sample.
  * **Point** - a single point
  * **RectUniform** - a uniformly sampled rectangle with user defined number of samples
  * **RectGrid** - a rectangle sampled in a grid fashion
  * **Hexapolar** - a circle (or an ellipse) sampled in an hexapolar fasion (rings)
- Rays Directions Distribution - the interface **length** returns the number of samples, and **generate** returns the n'th sample.
  * **Constant**
  * **RectGrid**
  * **UniformCone**
  * **HexapolarCone**


"""

# ╔═╡ d970a5b0-8858-11eb-1425-cdfb68e222ba
md"### Examples of basic emitters"

# ╔═╡ f21787a0-8858-11eb-19b1-914e835922eb
md"#### Point Origin Samples"

# ╔═╡ 5efca830-885b-11eb-3f76-3d9ac33369d8
md"#### Rectangle and Ellipse Origins samples"

# ╔═╡ 4d265a60-885c-11eb-2b6d-7b244ccfd172
md"#### Angular Power Distribution"

# ╔═╡ 4553fb50-88d7-11eb-3fc3-6ddfeb3f1910
md"In the following drawings, the length of the sampled rays is relative to its power."

# ╔═╡ fd3287f0-885f-11eb-279d-990c4222ba83
md"## Composite Emitters"

# ╔═╡ 03fc33b0-8860-11eb-2a20-0f8c2b1f79a3
md"""
**Composite Emitters** are emitters that can hold a list of basic emitters or composite ones. Composite emitters allow you to combine multiple basic emitters, such as an R, G and B emitters in to a Pixel emitter, which then can be combined in to a "Display" composite emitters.
The Transform data structure allow the user to define the geometric arrangment of the composited elements.

"""

# ╔═╡ 94bffc60-8860-11eb-37ca-b3203b46693a
md"### Display example"

# ╔═╡ b39c3952-8860-11eb-13a4-0b5004282a09
begin
	# construct the emitter's basic components
	local S = Spectrum.Uniform()
	local P = AngularPower.Lambertian()
	local O = Origins.RectGrid(1.0, 1.0, 3, 3)
	local D = Directions.HexapolarCone(deg2rad(5.0), 3)	
	
	# construct the srouce. in this example a "pixel" source will contain only one source as we are simulating a "b/w" display. 
	# for RGB displays we can combine 3 sources to simulate "a pixel".
	local Tr = Transform(Geometry.Vec3(0.5, 0.5, 0.0))
	local source1 = Sources.Source(Tr, S, O, D, P)
	
	# create a list of pixels - each one is a composite source
	local pixels = Vector{Sources.CompositeSource{Float64}}(undef, 0)
	for y in 1:10 # image_height
		for x in 1:10 # image_width
			# pixel position relative to the display's origin
			local pixel_position = Geometry.Vec3((x-1) * 1.1, (y-1) * 1.5, 0.0)
			local Tr = Transform(pixel_position)

			# constructing the "pixel"
			pixel = Sources.CompositeSource(Tr, [source1])

			push!(pixels, pixel)
		end
	end
	
	local Tr = Transform(Geometry.Vec3(0.0, 0.0, 0.0))
	my_display = Sources.CompositeSource(Tr, pixels)
	md"###### CONSTRUCTION CODE: Display is composed $(length(pixels)) pixels, esch one generating $(length(pixels[1])) rays, resulting in $(length(my_display)) total."	
end

# ╔═╡ 04d97450-8860-11eb-2bda-036792395ee4
md"----"

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
	OpticSim.NotebooksUtils.SetBackend(defs, makie_backend)
	md"Set Makie backend"
end

# ╔═╡ 96e423a0-885a-11eb-02a3-8704e8dbdab6
begin
	function Vis.scene(resolution = (1000, 1000))
		# @info "RG: Vis.Scene Replacement"
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

# ╔═╡ e92e18de-8847-11eb-2fdf-b90425c7327b
begin
	drawing	# notebook only - create dependency on drawing backend - comment if running in REPL
	
	local s = Sources.Source(origins=Origins.Point(), directions=Directions.RectGrid(π/4, π/4, 15, 15))
	Vis.draw(s, resolution=resolution(), debug=true)
end

# ╔═╡ 544d5870-8848-11eb-385e-afcce49d099d
begin
	drawing	# notebook only - create dependency on drawing backend - comment if running in REPL
	
	local s = Sources.Source(origins=Origins.Point(), directions=Directions.UniformCone(π/6, 1000))
	Vis.draw(s, resolution=resolution(), debug=true)
end

# ╔═╡ 78875510-8848-11eb-18d6-eba32897a99e
begin
	drawing	# notebook only - create dependency on drawing backend - comment if running in REPL
	
	local s = Sources.Source(origins=Origins.Point(), directions=Directions.HexapolarCone(π/6, 10))
	Vis.draw(s, resolution=resolution(), debug=true)
end

# ╔═╡ 69f94f40-885b-11eb-08b0-e9faf7969dda
begin
	drawing	# notebook only - create dependency on drawing backend - comment if running in REPL
	
	local s = Sources.Source(origins=Origins.RectGrid(1.0, 1.0, 10, 10), directions=Directions.Constant())
	Vis.draw(s, resolution=resolution(), debug=true)
end

# ╔═╡ ac0ee1b0-885b-11eb-0395-a19b37c5a472
begin
	drawing	# notebook only - create dependency on drawing backend - comment if running in REPL
	
	local s = Sources.Source(origins=Origins.RectGrid(1.0, 1.0, 3, 3), directions=Directions.HexapolarCone(π/6, 10))
	Vis.draw(s, resolution=resolution(), debug=true)
end

# ╔═╡ e14e85c0-8846-11eb-1f4a-2768e655321c
begin
	drawing	# notebook only - create dependency on drawing backend - comment if running in REPL
	
	local s = Sources.Source(origins=Origins.Hexapolar(2, 8.0, 8.0), directions=Directions.RectGrid(π/8, π/8, 5, 5))
	Vis.draw(s, resolution=resolution(), debug=true)
end

# ╔═╡ 69d61a60-885c-11eb-3f97-75e98121ec8c
begin
	drawing	# notebook only - create dependency on drawing backend - comment if running in REPL
	
	local s = Sources.Source(
		origins=Origins.Hexapolar(1, 8.0, 8.0), 
		directions=Directions.RectGrid(π/6, π/6, 15, 15), 
		# try to comment and un-comment the following lines to see the diffrence
		power=AngularPower.Cosine(10.0)
	)
	Vis.draw(s, resolution=resolution(), debug=true)
end

# ╔═╡ 0ae26f80-885d-11eb-2c2a-ad23b34daac1
begin
	drawing	# notebook only - create dependency on drawing backend - comment if running in REPL
	
	local s = Sources.Source(
		origins=Origins.RectGrid(1.0, 1.0, 3, 3), 
		directions=Directions.HexapolarCone(π/6, 10), 
		# try to comment and un-comment the following lines to see the diffrence
		power=AngularPower.Gaussian(2.0, 2.0)
	)
	Vis.draw(s, resolution=resolution(), debug=true)
end

# ╔═╡ 941bb0a0-8861-11eb-09fd-87e2022e068e
begin
	drawing	# notebook only - create dependency on drawing backend - comment if running in REPL
	
	Vis.draw(my_display;  resolution=resolution(), color=:red, debug=false)
	local rays = AbstractArray{OpticalRay{Float64, 3}}(collect(my_display)) 
	Vis.draw!(rays; color=:red, debug=false)
end

# ╔═╡ def239f0-87bc-11eb-2edb-2f859ac41bee
begin
	drawing	# notebook only - create dependency on drawing backend - comment if running in REPL
	
	local s1 = Sources.Source(sourcenum = 1, origins=Origins.Hexapolar(5, 8.0, 8.0), directions=Directions.Constant())
	
	# create a second source and rotate it by -5 degs
	local s2 = Sources.Source(sourcenum = 2, origins=Origins.Hexapolar(5, 8.0, 8.0), directions=Directions.Constant(), transform=Transform(zeros(Geometry.Vec3), rotationY(deg2rad(-8)) * unitZ3()))
	
	# create the "ray generator"
	local combined_sources = Sources.CompositeSource(Transform(Geometry.Vec3(0.0, 0.0, 10.0), unitZ3() * -1), [s1 s2])
	
	# and draw the system + the generated rays
	Vis.drawtracerays(sys, raygenerator = combined_sources,  resolution=resolution(), test = true, trackallrays = true, colorbysourcenum = true, drawgen = false)
end

# ╔═╡ Cell order:
# ╟─6ce01760-879e-11eb-2b13-4d3d07c4b4ce
# ╟─87954830-8854-11eb-2baf-d3b79b804e07
# ╟─9608ed90-8854-11eb-080d-4d1f60d596f5
# ╟─d970a5b0-8858-11eb-1425-cdfb68e222ba
# ╟─f21787a0-8858-11eb-19b1-914e835922eb
# ╠═e92e18de-8847-11eb-2fdf-b90425c7327b
# ╠═544d5870-8848-11eb-385e-afcce49d099d
# ╠═78875510-8848-11eb-18d6-eba32897a99e
# ╠═5efca830-885b-11eb-3f76-3d9ac33369d8
# ╠═69f94f40-885b-11eb-08b0-e9faf7969dda
# ╠═ac0ee1b0-885b-11eb-0395-a19b37c5a472
# ╠═e14e85c0-8846-11eb-1f4a-2768e655321c
# ╟─4d265a60-885c-11eb-2b6d-7b244ccfd172
# ╟─4553fb50-88d7-11eb-3fc3-6ddfeb3f1910
# ╠═69d61a60-885c-11eb-3f97-75e98121ec8c
# ╠═0ae26f80-885d-11eb-2c2a-ad23b34daac1
# ╟─fd3287f0-885f-11eb-279d-990c4222ba83
# ╟─03fc33b0-8860-11eb-2a20-0f8c2b1f79a3
# ╟─94bffc60-8860-11eb-37ca-b3203b46693a
# ╠═b39c3952-8860-11eb-13a4-0b5004282a09
# ╠═941bb0a0-8861-11eb-09fd-87e2022e068e
# ╟─04d97450-8860-11eb-2bda-036792395ee4
# ╟─0b6553c2-8860-11eb-0c5e-89632a8dd612
# ╟─a891a210-87bc-11eb-311d-4d020986fe19
# ╠═def239f0-87bc-11eb-2edb-2f859ac41bee
# ╟─5f77cfd0-8854-11eb-377c-ef975f9abf63
# ╟─67c8fac0-974e-11eb-2898-d14b79387954
# ╟─d6c83680-879e-11eb-31d4-7dbda7e93e48
# ╟─68e9b210-87ad-11eb-0f3a-5bb2dbf7d86c
# ╟─3a5d3ba0-87ae-11eb-1717-93be0b802cab
# ╟─96e423a0-885a-11eb-02a3-8704e8dbdab6
