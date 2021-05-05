
module AxisSymmetricOptimization

using OpticSim, OpticSim.Geometry, OpticSim.Emitters
using CSV, DataFrames
using Ipopt
using JuMP
import DataStructures

global global_model = nothing
global global_ax = nothing
global global_df = nothing
global global_counter = 0
global solution_history = []
global statistics = Dict()
global options = Dict()

function set_global_ax(ax) 
    global global_ax = ax
end 

global emitter = Sources.Source(origins=Emitters.Origins.Hexapolar(5, 8.0, 8.0), directions=Directions.Constant(), transform=Transform(Vec3(0.0,0.0,5.0), Geometry.unitZ3()*-1))

mutable struct VarInfo
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

function try_parse_variable_1(row, col, str)

    float_re = "[+-]?(?:[0-9]+(?:[.][0-9]*)?|[.][0-9]+)"

    re = Regex(
        "^\\s*"*                                    # white spaces at start
        "(?:"*
            "(?<lower>(?:$float_re))" *             # lower bound
            "\\s*"*                                 # white spaces 
            "(?<lower_cond><=?)" *                  # lower condition
        ")?"*
        "\\s*"*                                     # white spaces 
        "(?<name>(?:[@a-zA-Z][a-zA-Z0-9_@]*))" *    # variable name
        "\\s*"*                                     # white spaces 
        "(?<assignment>=)" *                        # assignment
        "\\s*"*                                     # white spaces 
        "(?<middle>(?:$float_re))" *                # starting value
        "\\s*"*                                     # white spaces 
        "(?:"*
            "(?<upper_cond><=?)" *                  # uppwe condition
            "\\s*"*                                 # white spaces 
            "(?<upper>(?:$float_re))" *             # upper bound
        ")?"*
        "\\s*\$"                                    # white spaces at end
    )
    
    m=match(re,str)
    if (m === nothing)
        return nothing
    end

    name = m["name"]
    start = parse(Float64, m["middle"])
    lower = m["lower"]===nothing ? nothing : parse(Float64, m["lower"])
    upper = m["upper"]===nothing ? nothing : parse(Float64, m["upper"])

    var_info = VarInfo(row, col, name, start, lower, upper, Dict())

    return var_info;
end


struct AxisSymmetric
    _original_definition
    _original_df::DataFrame
    _variables::AbstractArray
    
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

                    # try the new version of variable decleration
                    var_info = try_parse_variable_1(row, col, item)
                    if (var_info !== nothing)
                        if (var_info._name == "@")
                            auto_variable_count += 1
                            var_info._name = "autoname$(auto_variable_count)"
                        end
                    # try the old version of variable decleration
                    else
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
                    end

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

        if (mod(global_counter, get(options, :SnapshotRate, 50)) == 0)
            @info "INSIDE $(global_counter)"
            push!(solution_history, x)
            info = Dict()
            statistics["iterations"][global_counter] = Dict("index" => global_counter, "parameters" => x, "cost" => error)
            # show_optimized_system(x...)
        end
            
        return error
    catch e
        @error "Somthing went wrong in calculating the cost function\n$(e)"
        return 1000 # rqual inf
    end
end

function optimize(ax::AxisSymmetric, optimization_options::Dict = Dict())
    @info "Optimizing"

    global options = copy(optimization_options)

    model = Model(Ipopt.Optimizer)
    set_time_limit_sec(model, get(options, :MaxTime, 10.0))
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

    global global_counter = 0
    global global_model = model
    global solution_history = []
    global statistics = Dict()

    statistics["max_time"] = get(options, :MaxTime, 10.0)
    statistics["iterations"] = DataStructures.OrderedDict()

    original_params = AxisSymmetricOptimization.original_values(ax)
    cost_original = AxisSymmetricOptimization.cost_func(original_params...)
    statistics["iterations"]["original"] = Dict("index" => global_counter, "parameters" => original_params, "cost" => cost_original)

    JuMP.optimize!(model)

    # save final statistics
    optimized_params = AxisSymmetricOptimization.optimized_values(ax, model)
    cost_optimized = AxisSymmetricOptimization.cost_func(optimized_params...)
    statistics["iterations"]["final"] = Dict("index" => global_counter, "parameters" => optimized_params, "cost" => cost_optimized)

    # debug dump of variables
    if (true)
        vars = JuMP.all_variables(model)
        println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
        for v in vars
            println("[$(name(v))]=[$(value(v))]")        
        end
        println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
    end

    return model
end


function get_statistics()
    return copy(statistics)
end


end # module AxisSymmetricOptimization




