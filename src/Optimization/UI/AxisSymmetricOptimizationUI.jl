# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

module AxisSymmetricOptimizationUI

const INTERACTIVE_HTML = false

using OpticSim
using Blink, WebIO, Interact, Interact.CSSUtil, DataFrames, CSV, TableView
using ..BlinkUtils
using ..AxisSymmetricOptimization

#
import Makie
import Makie.AbstractPlotting
import Makie.AbstractPlotting.MakieLayout
if (INTERACTIVE_HTML)
    import WGLMakie
else
    import GLMakie
end


function rg_showtable_sync!(w, schema, names, types, rows, coldefs, tablelength, id, options)
    # RG - disable sorting for all columns
    for c in coldefs
        c[:sortable] = false
    end
      
    options[:rowData] = TableView.JSONText(TableView.table2json(schema, rows, types))
    license = get(ENV, "AG_GRID_LICENSE_KEY", nothing)
    handler = @js function (RowNumberRenderer, agGrid)
        @var gridOptions = $options
        @var el = document.getElementById($id)
        gridOptions.components = Dict(
            "rowNumberRenderer" => RowNumberRenderer
        )
        if $(license !== nothing)
            agGrid.LicenseManager.setLicenseKey($license)
        end
        this.table = @new agGrid.Grid(el, gridOptions)
        gridOptions.columnApi.autoSizeAllColumns()
    end
    onimport(w, handler)
end

function run()

    @eval begin
        # this function is needed to allow the visualization scene to be displayed inside a blink window
        function Makie.display(obj)
            # @info "RG: $obj"
            return obj
        end

        function Vis.scene(resolution = (1000, 1000))
            # @info "RG: Vis.Scene Replacement"
            scene, layout = MakieLayout.layoutscene(resolution = resolution)
            Vis.set_current_main_scene(scene)
            lscene = layout[1, 1] = MakieLayout.LScene(scene, scenekw = (camera = Makie.cam3d_cad!, axis_type = Makie.axis3d!, raw = false))
            Vis.set_current_3d_scene(lscene)
            return scene, lscene
        end

        function TableView._showtable_sync!(w, schema, names, types, rows, coldefs, tablelength, id, options)
            # RG - disable sorting for all columns
            for c in coldefs
                c[:sortable] = false
            end
            rg_showtable_sync!(w, schema, names, types, rows, coldefs, tablelength, id, options)  
        end
    end


    if (INTERACTIVE_HTML)
        WGLMakie.activate!()
    else
        GLMakie.activate!()
    end
    AbstractPlotting.__init__()
    AbstractPlotting.inline!(true)

    # a dummy call to a Makie rendering function just to make sure makie is compiled by Julia
    Makie.scatter(rand(10))

    # "global" variables for this function
    ax = nothing
    model = nothing
    statistics = nothing
    slider_iterations = nothing

    #
    # load the defualt file (name of the script with an ".html" extention)
    #

    # this is an optional list of stings to replace inside the html - somthing like an internal macro 
    # here we are trying to provide "absolute" path to local files relative to the files location
    a = Base.replace(@__FILE__, "\\" => "/")
    @show a
    replace = [                         
        """src="./""" => """src="$(Base.replace(splitdir(@__FILE__)[1], "\\" => "/"))/""",
    ]

    html = BlinkUtils.HTMLParserDefault(@__FILE__, replace_in_html = replace)
    
    # not nececery to change folder for now
    # BlinkUtils.ChangeDirectory(@__FILE__)         

    # function run_electron_js_script(script)
    #     script = replace(script, "@@@" => """windows["$(win.id)"]""")
    #     @info script
    #     js(win.shell, Blink.JSString(script), callback=false)
    # end


    function set_zoom_factor(val)
        @info "Zoom $val"
        # js(win.shell, Blink.JSString("""
        # windows["$(win.id)"].webContents.setZoomFactor($val);
        # """), callback=false)
        BlinkUtils.run_electron_js_script(win, "@@@.webContents.setZoomFactor($val)")
    end

    # Utils

    # i removed the zoom slide for now - keyboard shortcuts are much better for this type of activity
    # slider_zoom = slider([0.2, 0.5, 1.0, 1.5, 2, 3, 4, 5], value=1, label="Zoom")
    # on(slider_zoom) do val
    #     set_zoom_factor(val)
    # end
    # BlinkUtils.observable(html, "ZoomSlider")[] = dom"div"(slider_zoom)


    #
    # Prescription Filename
    #

    prescription_filename = Observable{Any}(String)
    on(prescription_filename) do val
        @info "Got a new prescription: [$val]"
   
        # html_val = replace(val, "\\" => "\\\\")
        # js(win, Blink.JSString("""document.getElementById("myspan").textContent="$html_val";"""), callback=false)
        BlinkUtils.observable(html, "Filename")[] = dom"div"(Node(:p, val))

        load_dataframe()
    end


    function data_cell_change(arr, msg)
        try
            row = msg["row"] + 1 # zero-indexed in JS
            col = msg["col"]
            @info "Cell Change: $row, $col"
            # col = parse(Int, match(r"Column(\d+)", msg["col"])[1])
            # arr[row, col] = parse(eltype(arr), msg["new"])
            arr[row, col] = msg["new"]
            # @info msg
        catch e
            show
            set_message("Error updating value\n$e")
            show_table()
        end
    end

    function show_table()
        table = TableView.showtable(data_table[], cell_changed = msg -> data_cell_change(data_table[], msg))

        set_visibility("div.rg-buttons", true)
        set_visibility("div.rg-parameters", true)
        set_visibility("div.drawings", true)
        set_visibility("div.history", true)

        BlinkUtils.observable(html, "DataFrame")[] = dom"div"(table)
    end

    #
    # data table
    #
    data_table = Observable{Any}(DataFrame)
    data_table_clean = Observable{Any}(DataFrame)
    on(data_table) do val
        @info "DATA TABLE Updated"
        if (isempty(data_table[]))
            BlinkUtils.observable(html, "DataFrame")[] = dom"span"(HTML("""<p style="text-align: center">No Data Loaded</p>"""))
            return 
        end

        show_table()

        validate_data()
        show_original_drawing()
    end

    #
    # Optimization parameters
    #
    text_optimization_time = textbox(value="10", hint="Optimization time")
    BlinkUtils.observable(html, "OptimizationTime")[] = dom"span"(text_optimization_time)

    text_optimization_snapshot_rate = textbox(value="50", hint="Optimization Snapshot Rate")
    BlinkUtils.observable(html, "OptimizationSnapshotRate")[] = dom"span"(text_optimization_snapshot_rate)


    #
    # Button Load
    #
    button_load = filepicker(label = "Load prescription to optimize", multiple=false, accept=".csv")
    on(button_load) do val
        @info 1
        prescription_filename[] = val
    end
    BlinkUtils.observable(html, "FilePicker")[] = dom"span"(button_load)

    #
    # Button Test
    #
    button_test = button("Load Double Convex Example")
    on(button_test) do val
        filename = abspath(joinpath(splitdir(@__FILE__)[1], "..", "Data", "DoubleConvex.csv"))
        prescription_filename[] = filename
    end
    BlinkUtils.observable(html, "TestButton")[] = dom"span"(button_test)

    #
    # Button Optimize
    #
    button_optimize = button("Optimize")
    on(button_optimize) do val
        optimize()
    end
    BlinkUtils.observable(html, "Optimize")[] = dom"span"(button_optimize)

    #
    # Button Validate
    #
    button_validate = button("Validate")
    on(button_validate) do val
        @info "Validate"
        validate_data()
    end
    BlinkUtils.observable(html, "Validate")[] = dom"span"(button_validate)

    #
    # Button External Application
    #
    button_external_app = button("External App Try")
    on(button_external_app) do val
        @info "ExternalApp"
        BlinkUtils.run_shell_default(prescription_filename[])
    end
    BlinkUtils.observable(html, "ExternalApp")[] = dom"span"(button_external_app)


    #
    # loading the data frame
    #
    function load_dataframe()
        # load the dataframe
        fn = prescription_filename[]
        if (startswith(fn, "?"))
            return
        end
        @info "loading csv file [$fn]"
        csv = CSV.File(fn)
        data_table[] = DataFrame(csv)
        @info "done loading csv file"
    end

    function clear_data()
        prescription_filename[] = "? No File Loaded ?"
        data_table[] = DataFrame()
        set_table_stat_color("inherit")
        set_visibility("div.rg-buttons", false)
        set_visibility("div.rg-parameters", false)
        set_visibility("div.drawings", false)
        set_visibility("div.history", false)

        BlinkUtils.observable(html, "OriginalStatistics")[] = dom"div"("")
        BlinkUtils.observable(html, "OriginalDrawing")[] = dom"div"("")
        BlinkUtils.observable(html, "OptimizedStatistics")[] = dom"div"("")
        BlinkUtils.observable(html, "OptimizedDrawing")[] = dom"div"("")

        BlinkUtils.observable(html, "OptimizationGraph")[] = dom"div"("")
        BlinkUtils.observable(html, "HistorySlider")[] = dom"div"("")
        BlinkUtils.observable(html, "HistoryDataFrame")[] = dom"div"("")

    end

    function set_visibility(id, val)
        display_style = (val) ? "block" : "none"
        js(win, Blink.JSString("""
        console.log("Setting Visibility of [$id] to [$display_style]")
        var items = document.querySelectorAll("$id")
        console.log("Found " + items.length.toString() + " items")
        for (var i = 0; i < items.length; i++) {
            console.log("Setting Visibility of " + items[i].className + " to [$display_style]")
            items[i].style.display = '$(display_style)'
        }
        """), callback=false)
    end

    function set_table_stat_color(val)
        js(win, Blink.JSString("""
        console.log("Setting Table Frame Color to [$val]")
        var items = document.querySelectorAll("div.data-table-frame")
        for (var i = 0; i < items.length; i++) {
            items[i].style.border = "medium solid $val"
        }
        """), callback=false)
    end

    function set_message(msg)
        BlinkUtils.observable(html, "Message")[] = dom"div"(Node(:p, msg))
    end

    function validate_data()
        @info "validating data"

        res = true;
        try
            ax = AxisSymmetricOptimization.AxisSymmetric(data_table[])
            set_table_stat_color("MediumSeaGreen")
            set_message("Validation Passed")
        catch e
            ax = nothing
            set_table_stat_color("red")
            set_message("Validation Error\n$e")
            res = false
        end
        AxisSymmetricOptimization.set_global_ax(ax)
        return res
    end

    function show_original_drawing()
        @info "showing original drawing"
        validate_data()

        df = AxisSymmetricOptimization.original(ax)
        sys_original = AxisymmetricOpticalSystem(df)

        if (INTERACTIVE_HTML)
            resolution = BlinkUtils.resolution(html, "OriginalDrawing")
            drawing_original  = Vis.drawtracerays(sys_original, raygenerator=AxisSymmetricOptimization.emitter, trackallrays = true, test = true, linewidth=100, resolution = resolution)
        else
            drawing_original  = Vis.drawtracerays(sys_original, raygenerator=AxisSymmetricOptimization.emitter, trackallrays = true, test = true, linewidth=1)
        end
        cost_original = AxisSymmetricOptimization.cost_func(AxisSymmetricOptimization.original_values(ax)...)

        stats = vbox(
            Node(:div, "Original State"), 
            Node(:div, "Cost: $cost_original"), 
        )
        BlinkUtils.observable(html, "OriginalStatistics")[] = dom"div"(stats)
        BlinkUtils.observable(html, "OriginalDrawing")[] = dom"div"(drawing_original)
    end

    # function show_optimized_drawing()
    #     @info "showing optimized drawing"

    #     sys_optimized = AxisymmetricOpticalSystem(AxisSymmetricOptimization.optimized(ax, model))

    #     if (INTERACTIVE_HTML)
    #         resolution = BlinkUtils.resolution(html, "OriginalDrawing")
    #         @info "@@@@@@@@@@@@@@@@@@@@@@ $resolution"
    #         drawing_optimized  = Vis.drawtracerays(sys_optimized, raygenerator=AxisSymmetricOptimization.emitter, trackallrays = true, test = true, linewidth=100, resolution=resolution)
    #     else
    #         drawing_optimized  = Vis.drawtracerays(sys_optimized, raygenerator=AxisSymmetricOptimization.emitter, trackallrays = true, test = true, linewidth=1)
    #     end
    #     cost_optimized = AxisSymmetricOptimization.cost_func(AxisSymmetricOptimization.optimized_values(ax, model)...)

    #     stats = vbox(
    #         Node(:div, "Optimized State"), 
    #         Node(:div, "Cost: $cost_optimized"), 
    #     )
    #     BlinkUtils.observable(html, "OptimizedStatistics")[] = dom"div"(stats)
    #     BlinkUtils.observable(html, "OptimizedDrawing")[] = dom"div"(drawing_optimized)
    # end

    # show specific history entry
    function show_stats(val)
        @info "showing stat [$val]  [$(slider_iterations[])]"

        # in the case the slider value changed already, no reason to render this iteration
        if (slider_iterations[] != val)
            return
        end

        iterations = statistics["iterations"]
        key = (val==0) ? "original" : ((val==999999999) ? "final" : val)

        # @show key
        iteration = iterations[key]

        params = iteration["parameters"]
        cost = iteration["cost"]

        if (!isa(params[1], Float64))
            params = Tuple(x.value for x in params)
        end
        # @show params, cost
        
        df = AxisSymmetricOptimization.optimized(ax, params)

        sys_optimized = AxisymmetricOpticalSystem(df)

        if (haskey(iteration, "drawing"))
            drawing = iteration["drawing"]    
            @info "used chache"
        else
            if (INTERACTIVE_HTML)
                resolution = BlinkUtils.resolution(html, "OptimizedDrawing")
                drawing  = Vis.drawtracerays(sys_optimized, raygenerator=AxisSymmetricOptimization.emitter, trackallrays = true, test = true, linewidth=100, resolution=resolution)
            else
                drawing  = Vis.drawtracerays(sys_optimized, raygenerator=AxisSymmetricOptimization.emitter, trackallrays = true, test = true, linewidth=1)
            end
            iteration["drawing"] = drawing
        end

        stats = vbox(
            Node(:div, "Optimized State [Iteration: $key]"), 
            Node(:div, "Cost: $cost"), 
        )
        BlinkUtils.observable(html, "OptimizedStatistics")[] = dom"div"(stats)
        BlinkUtils.observable(html, "OptimizedDrawing")[] = dom"div"(drawing)

        table = TableView.showtable(df)
        BlinkUtils.observable(html, "HistoryDataFrame")[] = dom"div"(table)
    end

    function show_stats()
        @info "showing statistics"

        iterations = statistics["iterations"]
        slider_range = [ (k=="original") ? 0 : ((k=="final") ? 999999999 : (k))  for (k, v) in iterations]

        # prepare the optimization graph
        x = [Float64(v["index"]) for (k, v) in iterations]
        y = [v["cost"] for (k, v) in iterations]
        # @show x
        # @show y

        resolution = BlinkUtils.resolution(html, "OptimizationGraph")
        s = nothing
        Makie.with_theme(Makie.Theme(resolution = resolution)) do
            s = Makie.scatter(x, y, color = :green, markersize = 5)
            Makie.lines!(x, y, color = :blue)
        end            

        BlinkUtils.observable(html, "OptimizationGraph")[] = dom"div"(Makie.display(s))


        # history slider
        slider_iterations = slider(slider_range, value=999999999, label="Iterations History")
        on(slider_iterations) do val
            show_stats(val)
        end
        BlinkUtils.observable(html, "HistorySlider")[] = dom"div"(slider_iterations)

        show_stats(999999999)
    end

    function optimize()
        @info "Optimize"

        if (!validate_data())
            return
        end

        max_optimization_time = 0.0
        try
            max_optimization_time = parse(Float64, text_optimization_time[])
        catch e
            set_message("Error parsing optimization time\n$e")
            return
        end

        optimization_snapshot_rate = 10
        try
            optimization_snapshot_rate = parse(Float64, text_optimization_snapshot_rate[])
        catch e
            set_message("Error parsing optimization snapshot rate\n$e")
            return
        end

        optimization_options = Dict(
            :MaxTime => max_optimization_time,
            :SnapshotRate => optimization_snapshot_rate,
        )

        if (true)
            model = AxisSymmetricOptimization.optimize(ax, optimization_options)

            # show_optimized_drawing()

            statistics = AxisSymmetricOptimization.get_statistics()
            # @show statistics

            show_stats()
        end

    end

    # show the html in a blink window - note that you can use the title in the html file to overrigde the title supplied by this call
    win = BlinkUtils.window(html, title="OpticSim: Axis Simmetric Optimization", size=(1400, 800))

    # toggle the dev tools - can be done using F12 now
    # Blink.AtomShell.tools(win)    

    set_zoom_factor(1.0)
    BlinkUtils.setup_default_keybinds(win)

    #
    # message handlers form Julia
    #

    handle(win, "OpenLink") do args  
        BlinkUtils.run_shell_default(args)
    end


    # if (INTERACTIVE_HTML)
    #     @eval begin
    #         using JSServe

    #         function JSServe.iframe_html(server::JSServe.Server, session::JSServe.Session, route::String)
    #             @info "My JS Serve"
    #             # Display the route we just added in an iframe inline:
    #             url = JSServe.online_url(server, route)
    #             remote_origin = JSServe.online_url(server, "")
    #             style = "position: absolute; width: 100%; height: 100%; padding: 0; overflow: hidden; border: none"
    #             return DOM.iframe(src=url, id=session.id, style=style, scrolling="no")
    #         end
    #     end
    # end

    # throttle()
    Makie.scatter

    clear_data()
    nothing;
end

end # module

