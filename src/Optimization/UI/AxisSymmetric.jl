
#region AxisSymmetricOptimization 
module AxisSymmetricOptimization

using OpticSim, OpticSim.Geometry, OpticSim.Emitters
using CSV, DataFrames
using Ipopt
using JuMP

global global_model = nothing
global global_ax = nothing
global global_df = nothing
global global_counter = 0
global solution_history = []

function set_global_ax(ax) 
    global global_ax = ax
end 

global emitter = Sources.Source(origins=Emitters.Origins.Hexapolar(5, 8.0, 8.0), directions=Directions.Constant(), transform=Transform(Vec3(0.0,0.0,5.0), Geometry.unitZ3()*-1))

struct VarInfo
    _row::Int
    _col::Int
    _name::String
    _start_value::Float64
    _lower_bound::Union{Float64,Nothing}
    _upper_bound::Union{Float64,Nothing}
    _info::Dict
end

function get_number(dic, key)
    v = get(dic, key, Nothing)
    if (v != Nothing)
        v = parse(Float64, v)
    end
    return v
end

struct AxisSymmetric
    _original_definition
    _original_df::DataFrame
    _variables::AbstractArray
    
    function AxisSymmetric(def::String)
        
        # load the dataframe from the input string
        df = CSV.read(IOBuffer(def), DataFrame, delim="|", normalizenames=true, missingstrings=["missing", "NA"])
        
        # remove spaces from column names
        rename!(x -> strip(x), df)

        # parse the objects names - if name starts with : make it a symbol, otherwise treat it as an integer
        # df.Surface = [(startswith(strip(x), ":")) ? Symbol(strip(x)[2:end]) : parse(Int64, strip(x))  for x in df.Surface]

        # parse the materials names - add the OpticSim.GlassCat. prefix if needed and try to evaluate it
        # df.Material = [(startswith(strip(x), "OpticSim")) ? eval(Meta.Parse(strip(x))) : eval(Meta.parse("OpticSim.GlassCat.$(strip(x))")) for x in df.Material]

        auto_variable_count = 0
        variables = []

        # first loop - fixing CSV issues and converting items 
        for col in axes(df,2)
            col_name = names(df)[col]
            col_vector = df[!, col]
            new_col_vector = []
            for row in 1:length(col_vector)
                item = col_vector[row]
                if (col_name == "Surface")
                    if (isa(item, String))
                        item = strip(item)
                        item = startswith(item, ":") ? Symbol(item[2:end]) : parse(Int64, item)
                    end
                end
                if (col_name == "Material")
                    if (isa(item, String))
                        item = strip(item)
                        item = startswith(item, "OpticSim") ? eval(Meta.Parse(item)) : eval(Meta.parse("OpticSim.GlassCat.$(item)"))
                    end
                end
                if (col_name == "Radius")
                    if (isa(item, String))
                        item = strip(item)
                        item = (occursin(",", item)||occursin("=", item)) ? String(item) : parse(Float64,item)
                    end
                end
                if (col_name == "Thickness")
                    if (isa(item, String))
                        item = strip(item)
                        item = (occursin(",", item)||occursin("=", item)) ? String(item) : ((item=="missing") ? missing : parse(Float64,item))
                    end
                end
                push!(new_col_vector, item)
            end
            df[!, col] = new_col_vector
        end

        # loop over dataframe
        for row in axes(df,1)
            for col in axes(df,2)
                col_name = names(df)[col]
                item = df[row, col]
                if (isa(item, String))
                    parts = split(item, ",")
                    dic = Dict()
                    for part in parts
                        kv = split(part, "=")
                        dic[strip(kv[1])] = strip(kv[2])
                    end

                    auto_variable_count += 1
                    name = get(dic, "name", "autoname$(auto_variable_count)")
                    start= get_number(dic, "start")
                    lower= get_number(dic, "low")
                    upper= get_number(dic, "high")

                    var_info = VarInfo(row, col, name, start, lower, upper, Dict())
                    push!(variables, var_info)
                end
            end
        end

        return new(def, df, variables)
    
    
    

    end

    function AxisSymmetric(def::DataFrame)

        df = copy(def)

        auto_variable_count = 0
        variables = []

        # first loop - fixing CSV issues and converting items 
        for col in axes(df,2)
            col_name = names(df)[col]
            col_vector = df[!, col]
            new_col_vector = []
            for row in 1:length(col_vector)
                item = col_vector[row]
                if (col_name == "Surface")
                    if (isa(item, String))
                        item = strip(item)
                        item = startswith(item, ":") ? Symbol(item[2:end]) : parse(Int64, item)
                    end
                end
                if (col_name == "Material")
                    if (isa(item, String))
                        item = strip(item)
                        item = startswith(item, "OpticSim") ? eval(Meta.Parse(item)) : eval(Meta.parse("OpticSim.GlassCat.$(item)"))
                    end
                end
                if (col_name == "Radius")
                    if (isa(item, String))
                        item = strip(item)
                        item = (occursin(",", item)||occursin("=", item)) ? String(item) : parse(Float64,item)
                    end
                end
                if (col_name == "Thickness")
                    if (isa(item, String))
                        item = strip(item)
                        item = (occursin(",", item)||occursin("=", item)) ? String(item) : ((item=="missing") ? missing : parse(Float64,item))
                    end
                end
                push!(new_col_vector, item)
            end
            df[!, col] = new_col_vector
        end

        # loop over dataframe
        for row in axes(df,1)
            for col in axes(df,2)
                col_name = names(df)[col]
                item = df[row, col]
                if (isa(item, String))
                    parts = split(item, ",")
                    dic = Dict()
                    for part in parts
                        kv = split(part, "=")
                        dic[strip(kv[1])] = strip(kv[2])
                    end

                    auto_variable_count += 1
                    name = get(dic, "name", "autoname$(auto_variable_count)")
                    start= get_number(dic, "start")
                    lower= get_number(dic, "low")
                    upper= get_number(dic, "high")

                    var_info = VarInfo(row, col, name, start, lower, upper, Dict())
                    push!(variables, var_info)
                end
            end
        end

        return new(def, df, variables)
    end
end

function original(ax::AxisSymmetric)
    res = copy(ax._original_df)
    for var in ax._variables
        res[var._row, var._col] = var._start_value
    end    
    res[!,:Radius] = passmissing(convert).(Float64,res[!,:Radius])
    res[!,:Thickness] = passmissing(convert).(Float64, res[!,:Thickness])
    res[!,:SemiDiameter] = passmissing(convert).(Float64,res[!,:SemiDiameter])
    return res
end

function original_values(ax::AxisSymmetric)
    res = []
    for var in ax._variables
        push!(res, var._start_value)
    end    
    return tuple(res...)
end

function optimized(ax::AxisSymmetric, params)
    T = typeof(params[1])
    res = copy(ax._original_df)
    for (index, var) in enumerate(ax._variables)
        res[var._row, var._col] = params[index]
    end    
    res[!,:Radius] = passmissing(convert).(T,res[!,:Radius])
    res[!,:Thickness] = passmissing(convert).(T, res[!,:Thickness])
    res[!,:SemiDiameter] = passmissing(convert).(T,res[!,:SemiDiameter])
    return res
end

function optimized(ax::AxisSymmetric, model::AbstractModel)
    T = Float64
    res = copy(ax._original_df)
    model_vars = JuMP.all_variables(model)
    for (index, var) in enumerate(ax._variables)
        res[var._row, var._col] = value(model_vars[index])
    end    
    res[!,:Radius] = passmissing(convert).(T,res[!,:Radius])
    res[!,:Thickness] = passmissing(convert).(T, res[!,:Thickness])
    res[!,:SemiDiameter] = passmissing(convert).(T,res[!,:SemiDiameter])
    return res
end

function optimized_values(ax::AxisSymmetric, model::AbstractModel)
    res = []
    model_vars = JuMP.all_variables(model)
    for (index, var) in enumerate(ax._variables)
        push!(res, value(model_vars[index]))
    end    
    return tuple(res...)
end

function cost_func(x...) 
    global global_counter+=1
    if (mod(global_counter, 50) == 0)
        @info "INSIDE $(global_counter)"
        push!(solution_history, x)
        # show_optimized_system(x...)
    end
    try
        df = optimized(global_ax, x)  

        # global global_df = df
    
        T = typeof(x[1])
        lens = AxisymmetricOpticalSystem{T}(df)
    
        field = HexapolarField(lens, collimated = true, samples = 10)
        error = zero(T)
        hits = 0
        for r in field
            traceres = OpticSim.trace(lens, r, test = true)
    
            if !(nothing === traceres)
                hitpoint = point(traceres)
                if abs(hitpoint[1]) > eps(T) && abs(hitpoint[2]) > eps(T)
                    dist_to_axis = sqrt(hitpoint[1]^2 + hitpoint[2]^2)
                    error += dist_to_axis
                end
                hits += 1
            end
        end
    
        # if isa(radius1,ForwardDiff.Dual)
        #     println("radius1 $(realpart(radius1)) radius2 $(realpart(radius2))")
        # end
    
        error /= hits
        # println("error in my code $error")
    
        # Vis.drawtracerays(doubleconvex(radiu1,radius2),trackallrays = true, test = true)
    
        return error
    catch e
        @error "Somthing went wrong\n$(e)"
        return 1000 # rqual inf
    end
end

function optimize(ax::AxisSymmetric)
    @info "Optimizing"
    model = Model(Ipopt.Optimizer)
    set_time_limit_sec(model, 10.0)
    # @variable(model, RG[1:1], container=Array)

    params = Vector{VariableRef}()
    for var in ax._variables
        has_lb = (var._lower_bound !== nothing)
        has_ub = (var._upper_bound !== nothing)
    
        var_info = VariableInfo(has_lb, var._lower_bound, has_ub, var._upper_bound, false, NaN, true, var._start_value, false, false)
        new_var = JuMP.add_variable(model, JuMP.build_variable(error, var_info), var._name)
        push!(params, new_var)
    end

    register(model, :cost_func, length(params), cost_func, autodiff = true)
    # @variable(model, 10.0 <= rad1 <= 85.0, start = 60.0)
    # @variable(model, -150 <= rad2 <= 20.0, start = -10.0)
    if (length(params) == 1)
        @NLobjective(model, Min, cost_func(params[1]))
    else
        @NLobjective(model, Min, cost_func(params...))
    end

    global global_model = model
    global solution_history = []
    JuMP.optimize!(model)

    vars = JuMP.all_variables(model)
    println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
    for v in vars
        println("[$(name(v))]=[$(value(v))]")        
    end
    println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")

    return model
end

end # module AxisSymmetricOptimization

#endregion AxisSymmetricOptimization 

#######################################################################################
module AxisSymmetricOptimizationUI

using OpticSim
using Blink, WebIO, Interact, Interact.CSSUtil, DataFrames, CSV, TableView
using ..AxisSymmetricOptimization

#
import Makie
import Makie.AbstractPlotting
import Makie.AbstractPlotting.MakieLayout
import GLMakie


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
    end

    GLMakie.activate!()
    AbstractPlotting.__init__()
    AbstractPlotting.inline!(true)

    #
    window_defaults = Blink.@d(
        :title => "Ran Gal", 
        :width=>1200, 
        :height=>800,
    )

    win = Window(window_defaults)

    header = HTML("""
    <style type="text/css">
    .auto-style2 {
        border-style: solid;
        border-width: 3px;
    }
    .auto-style3 {
        font-family: Arial, Helvetica, sans-serif;
        text-align: center;
        font-size: xx-large;
    }
    </style>

    <p class="auto-style3"><strong>Optimization UI</strong></p>
    """)

    perscription_filename = Observable{Any}(String)

    button_load = filepicker(label = "Choose a perscription to optimize", multiple=false, accept=".csv")
    button_test = button("Load Double Convex Example")
    span_filename = HTML("""<span id="myspan"> ?filename? </span>""")
    buttons = vbox(hbox(button_load, button_test), span_filename)
    # data = Observable{Any}(DataFrame)


    on(button_test) do val
        # @info 1
        filename = abspath(joinpath(splitdir(@__FILE__)[1], "..", "Data", "DoubleConvex.csv"))
        perscription_filename[] = filename
    end

    on(button_load) do val
        @info 3
        perscription_filename[] = val
    end

    # perscription filename changes
    on(perscription_filename) do val
        @info "Got an update: ", val
        html_val = replace(val, "\\" => "\\\\")
        js(win, Blink.JSString("""document.getElementById("myspan").textContent="$html_val";"""), callback=false)

        # load the dataframe
        @info "loading csv"
        csv = CSV.File(val)
        data_table[] = DataFrame(csv)
        @info "done"
    end

    div_data_table = Observable{Any}(dom"div"())
    div_drawing_1 = Observable{Any}(dom"div"())
    div_drawing_2 = Observable{Any}(dom"div"())
    data_table = Observable{Any}(DataFrame)

    function data_cell_change(arr, msg)
        row = msg["row"] + 1 # zero-indexed in JS
        col = msg["col"]
        # col = parse(Int, match(r"Column(\d+)", msg["col"])[1])
        # arr[row, col] = parse(eltype(arr), msg["new"])
        arr[row, col] = msg["new"]
        @info "Cell Change: $row, $col"
        @info msg
    end

    # data was loaded
    on(data_table) do val
        # @info "Got an update to data: ", val

        table = TableView.showtable(data_table[], cell_changed = msg -> data_cell_change(data_table[], msg))


        button_optimize = button("Optimize")
        button_validate = button("Validate (test)")

        buttons = hbox(button_optimize, button_validate)

        all = vbox(table, buttons)

        on(button_optimize) do val
            optimize()        
        end

        div_data_table[] = dom"div"(all)
    end


    function optimize()
        @info "Optimize"
        ax = AxisSymmetricOptimization.AxisSymmetric(data_table[])
        AxisSymmetricOptimization.set_global_ax(ax)

        df = AxisSymmetricOptimization.original(ax)
        sys_original = AxisymmetricOpticalSystem(df)
        cost_original = AxisSymmetricOptimization.cost_func(AxisSymmetricOptimization.original_values(ax)...)
        
        # field = HexapolarField(sys_original, collimated = true, samples = 4)
        
        if (true)

            m = AxisSymmetricOptimization.optimize(ax)
        
            sys_optimized = AxisymmetricOpticalSystem(AxisSymmetricOptimization.optimized(ax, m))
            cost_optimized = AxisSymmetricOptimization.cost_func(AxisSymmetricOptimization.optimized_values(ax, m)...)
        
            # println("Before Optimization: [$(AxisSymmetricOptimization.original_values(ax))]  Error: $(AxisSymmetricOptimization.cost_func(AxisSymmetricOptimization.original_values(ax)...))")
            # println("After Optimization:  [$(AxisSymmetricOptimization.optimized_values(ax, m))]  Error: $(AxisSymmetricOptimization.cost_func(AxisSymmetricOptimization.optimized_values(ax, m)...))")
        
            #Vis.drawtracerays(sys_optimized, raygenerator=emitter, trackallrays = true, test = true)
        end
        
        drawing_original  = Vis.drawtracerays(sys_original, raygenerator=AxisSymmetricOptimization.emitter, trackallrays = true, test = true, linewidth=1)
        drawing_optimized = Vis.drawtracerays(sys_optimized, raygenerator=AxisSymmetricOptimization.emitter, trackallrays = true, test = true)

        original = vbox(
            Node(:div, "Original State"), 
            Node(:div, "Cost: $cost_original"), 
            drawing_original
        )
        optimized = vbox(
            Node(:div, "Optimized State"), 
            Node(:div, "Cost: $cost_optimized"), 
            drawing_optimized
        )

        drawing_compared = hbox(original, optimized)

        div_drawing_1[] = Node(:div, drawing_compared)

        # @info typeof(drawing1)
        # div_drawings[] = dom"div"("Drawing 1", drawing1)
        # div_drawing_1[] = Node(:div, "Original State", drawing1)

        @info "Optimize Done"

    end


    ui = dom"div"(header, buttons, div_data_table, div_drawing_1, div_drawing_2)

    body!(win, ui)
    # Blink.AtomShell.opentools(win)

    return win
end    


end # module


