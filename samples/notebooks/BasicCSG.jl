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

# ╔═╡ 576fb1c0-8ba3-11eb-1c10-fdb51172e2c9
md"# Basic CSG (more to come soon)"

# ╔═╡ aa47f750-8c1f-11eb-21d7-d3a969fb1f6d
begin
	csg1_top_surface_caption = "Top Surface"
	csg1_bottom_surface_caption = "Bottom Surface"
	csg1_cylinder_surface_caption = "Cylinder"
	csg1_result_surface_caption = "Intersection Result"
	csg1_surfaces = [
		csg1_top_surface_caption, 
		csg1_bottom_surface_caption, 
		csg1_cylinder_surface_caption, 
		csg1_result_surface_caption,
	]

	NB.DefsClearHTML(defs)
	
	csg1_selected_surfaces_info = NB.GetVarInfo(@bind csg1_selected_surfaces PlutoUI.MultiSelect(csg1_surfaces, default=csg1_surfaces))
	
	cyl_rot_x_info = NB.GetVarInfo(@bind cyl_rot_x NB.UISlider(0:30, 0))
	cyl_rot_y_info = NB.GetVarInfo(@bind cyl_rot_y NB.UISlider(0:30, 0))
	cyl_rot_z_info = NB.GetVarInfo(@bind cyl_rot_z NB.UISlider(0:30, 0))

	NB.DefsAddHTML(defs, 
		"<div>" * csg1_selected_surfaces_info.html *"</div>" *
		"Rotation: (x, y, z):" * cyl_rot_x_info.html * 
		cyl_rot_y_info.html *
		cyl_rot_z_info.html 
	)
		
	PlutoUI.Show(MIME"text/html"(), NB.DefsHTML(defs))	
end

# ╔═╡ 5cc45e60-8c1f-11eb-3e4a-17c26a3214f6
begin
	function beziersurface()
		points = map(
			x -> collect(x),
			[
				(0.0, 0.0, 0.0) (0.0, 0.33, 0.0) (0.0, 0.66, 0.0) (0.0, 1.0, 0.0)
				(0.33, 0.0, 0.0) (0.33, 0.33, 1.0) (0.33, 0.66, 1.0) (0.33, 1.0, 0.0)
				(0.66, 0.0, 0.0) (0.66, 0.33, 1.0) (0.66, 0.66, 1.0) (0.66, 1.0, 0.0)
				(1.0, 0.0, 0.0) (1.0, 0.33, 0.0) (1.0, 0.66, 0.0) (1.0, 1.0, 0.0)
			],
		)
		return BezierSurface{OpticSim.Euclidean,Float64,3,3}(points)
	end	
	
	# canonic bezier surface
	csg1_surf1 = AcceleratedParametricSurface(beziersurface(), 25);
	
	# two transformed copies of the canonic bezier surface
	csg1_surf2 = leaf(csg1_surf1, OpticSim.translation(-0.5, -0.5, 0.0))
	csg1_surf3 = leaf(csg1_surf1, RigidBodyTransform{Float64}(0.0, Float64(π), 0.0, 0.5, -0.5, 0.0))
	
	# transformed cilinder
	csg1_surf4_1 = leaf(Cylinder(0.3, 1.0), OpticSim.translation(0.0, 0.0, 0.0))
	csg1_surf4 = leaf(csg1_surf4_1, OpticSim.rotation(deg2rad(cyl_rot_x), deg2rad(cyl_rot_y), deg2rad(cyl_rot_z)))
	
	# intersection result	
	csg1_surf5 = OpticSim.csgintersection(csg1_surf2,csgintersection(csg1_surf4, csg1_surf3))()
	
	md"## Define Surfaces and Perform CSG Optrations (code)"
end


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

# ╔═╡ 8ba2f700-8c1f-11eb-209f-b76b9713b576
begin
	drawing	# notebook only - create dependency on drawing backend - comment if running in REPL
	
	Vis.draw(OpticSim.SVector(0.0, 0.0, 0.0), markersize=0)

	if (csg1_top_surface_caption in csg1_selected_surfaces)
		Vis.draw!(
			csg1_surf2; 
			wireframe=true, 
			linewidth=1, 
			color=:orange, 
			shaded=false, 
			normals=false, 
			numdivisions=50
		)
	end
	if (csg1_bottom_surface_caption in csg1_selected_surfaces)
		Vis.draw!(
			csg1_surf3; 
			wireframe=true, 
			linewidth=1, 
			color=:blue, 
			shaded=false, 
			normals=false, 
			numdivisions=50
		)
	end
	if (csg1_cylinder_surface_caption in csg1_selected_surfaces)
		Vis.draw!(
			csg1_surf4; 
			wireframe=true, 
			linewidth=1, 
			color=:red, 
			shaded=false, 
			normals=false, 
			numdivisions=50
		)
	end
	if (csg1_result_surface_caption in csg1_selected_surfaces)
		Vis.draw!(
			csg1_surf5; 
			wireframe=true, 
			linewidth=1, 
			color=:green, 
			shaded=true, 
			normals=false, 
			numdivisions=50
		)
	end
	Vis.current_main_scene
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
# ╟─576fb1c0-8ba3-11eb-1c10-fdb51172e2c9
# ╟─5cc45e60-8c1f-11eb-3e4a-17c26a3214f6
# ╟─aa47f750-8c1f-11eb-21d7-d3a969fb1f6d
# ╟─8ba2f700-8c1f-11eb-209f-b76b9713b576
# ╟─5f77cfd0-8854-11eb-377c-ef975f9abf63
# ╟─d6c83680-879e-11eb-31d4-7dbda7e93e48
# ╟─68e9b210-87ad-11eb-0f3a-5bb2dbf7d86c
# ╟─3a5d3ba0-87ae-11eb-1717-93be0b802cab
# ╟─96e423a0-885a-11eb-02a3-8704e8dbdab6
