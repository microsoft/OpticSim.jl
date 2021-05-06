# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

module BlinkUtils

# using Blink, WebIO, Interact, Interact.CSSUtil, DataFrames, CSV, TableView, Widgets
using Blink, WebIO



# function WebIO.render(x::HTML{String}) 
#     # @info "--------------------------- $x"
#     # return dom"div"(; setInnerHtml=WebIO.richest_html(x))
#     return Node(:span ; setInnerHtml=WebIO.richest_html(x))
# end

# function WebIO.richest_html(val)
#     # @info "RG $val"
#     topmime = WebIO.richest_mime(val)
#     str_repr = WebIO.stringmime(topmime, val)
#     if topmime == MIME("text/html")
#         # @info "1 $str_repr"
#         return str_repr # |> WebIO.encode_scripts
#     elseif topmime in map(MIME, ["image/png", "image/jpeg"])
#         @info "2"
#         return "<img src='data:image/png;base64,$str_repr'></img>"
#     elseif topmime == MIME("image/svg+xml")
#         @info "3"
#         return str_repr
#     else
#         @info "4 $str_repr" 
#         return str_repr # |> WebIO.encode_scripts
#         # return "<pre>$str_repr</pre>"
#     end
# end


struct HTMLParser
    _original_text
    _text
    _variables
end

function HTMLParser(filename::String; kwargs...)
    sl = readlines(filename)
    return BlinkUtils.HTMLParser(sl; kwargs...)
end

function HTMLParser(text::AbstractArray{String}; replace_in_html = nothing, kwargs...)
    original_text = copy(text)
    text = []
    variables = Dict()
    white_spaces = [' ', '\t', '\r', '\n']
    index_offset = 0
    for (i, s) in enumerate(original_text)
        # skip a weird ID for html files that exists in windows - no clue about linux
        if (i == 1 && s[1] == '\ufeff')
            s = s[4:end]
        end

        ss = strip(s, white_spaces )
        
        # replace strings in html
        if (replace_in_html !== nothing)
            for p in replace_in_html 
                # @info p
                s = replace(s, p)
            end
        end

        if (startswith(ss, "<!--Julia("))
            push!(text, s)
            # @info "Found: [$s]"

            # re = Regex("\<\!--Julia\((?<dic>*.)\)--\>")
            re = Regex(raw"<!--Julia\((?<dic>.*)\)-->")
            m=match(re,ss)
            # @info m

            if (m === nothing)
                @error "Can't parse the line [$ss]"
                continue
            end
            dic_text = m["dic"]
            parts = split(dic_text, ",")
            dic = Dict()
            for part in parts
                kv = split(strip(part, white_spaces), "=")
                dic[String(strip(kv[1], white_spaces))] = String(strip(kv[2], white_spaces))
            end
            dic["original_index"] = i
            index_offset += 1
            dic["index"] = i + index_offset

            if (get(dic, "type", "?") == "div" )
                obs = Observable{Any}(Node(:div, dic["name"]))
                scope = Scope(dom=obs)
                dic["obs"] = obs
                dic["scope"] = scope
                # push!(text, scope)
                id = WebIO.scopeid(scope)
                push!(text, """<span class="rg-hook" id="$id">hook</span>""")

                variables[get(dic, "name", "?")] = dic
                @info "Variable [$(get(dic, "name", "?"))]"
            elseif (get(dic, "type", "?") == "span" )
                obs = Observable{Any}(Node(:span, dic["name"]))
                scope = Scope(dom=obs)
                dic["obs"] = obs
                dic["scope"] = scope
                # push!(text, scope)
                id = WebIO.scopeid(scope)
                push!(text, """<span class="rg-hook" id="$id">hook</span>""")

                variables[get(dic, "name", "?")] = dic
                @info "Variable [$(get(dic, "name", "?"))]"
            else
                @error "Can't parse the type of [$ss]"
                continue
            end

            # @info dic
        else
            push!(text, s)
        end
    end
    return HTMLParser(original_text, text, variables)
end

function HTMLParserDefault(script_filename::String; kwargs...)
    fn = splitext(script_filename)[1] * ".html"
    return HTMLParser(fn; kwargs...)
end

# don't think we need this function as electron defines its own rool location
# i'm planning to handle local files with macro replacements in the html code instead of trying to understand the electron "relative" location
function ChangeDirectory(script_filename::String)
    folder = splitdir(script_filename)[1]
    cd(folder)
    @info "Working directory is now [$(pwd())]"
end

variables(obj::HTMLParser) = keys(p._variables)
observable(obj::HTMLParser, name::String) = obj[name]["obs"]
function node(obj::HTMLParser) 
    res = []

    # insert all observables
    for (k,v) in obj._variables
        push!(res, v["scope"])
    end

    ss = ""
    for (i, s) in enumerate(obj._text)
        if (isa(s, String))
            if (ss == "")
                ss = s
            else
                ss = ss * "\n" * s
            end
        else
            if (ss != "")
                push!(res, HTML(ss))
                ss = ""
            end
            push!(res, s)
            # @info s
        end
        # if (i == 41 || i == 64)
        #     push!(res, HTML(ss))
        #     ss = ""
        # end
    
    end
    if (ss != "")
        push!(res, HTML(ss))
        ss = ""
    end

    # for i in res
    #     @info typeof(i)
    #     @info i
    # end
    return Node(:html, res...)
end

function Base.getindex(obj::HTMLParser, key::String)
    return obj._variables[key]
end


function window(node::HTMLParser;
                title = "Blink Window",
                size = (800, 600),
                )

    window_defaults = Blink.@d(
        :title => title, 
        :width=>size[1], 
        :height=>size[2],
        # this will allow us to load local file which is a security risk
        :webPreferences => Blink.@d(:webSecurity => false)
    )
    
    win = Window(window_defaults)
            
    ui = BlinkUtils.node(node)

    myscope = Scope(
        dom=ui,
        # imports=[Asset("foo.js"), "bar.css", "spam" => "spam.js"],
        mount_callbacks = Blink.JSString[
            Blink.JSString(""" 
            function ()
            {  
              var  hooks = document.querySelectorAll("span.rg-hook");
              console.log("Found " + hooks.length.toString() + " hooks")
              var obs = document.querySelectorAll("span.webio-scope");
              console.log("Found " + obs.length.toString() + " observables")
    
              /* move all hooks */
              for (var i = 0; i < hooks.length; i++) {
                var hook_id = hooks[i].getAttribute("id");
                console.log("hook ID " + hook_id)
                for (var j = 0; j < obs.length; j++) {
                  var ob_id = obs[j].getAttribute("data-webio-scope-id");
                  console.log("OB ID " + ob_id)
                  if (ob_id == hook_id)
                  {
                    hooks[i].innerHTML = "";
                    hooks[i].appendChild(obs[j]);
                    console.log("RG HOOk MOVED " + ob_id )
                    break;
                  }
                }
              }
    
              /* remove all place holders */
              var ph = document.querySelectorAll("span.rg-hook-placeholder")
              for (var i = 0; i < ph.length; i++) {
                console.log("removing place-holder")
                var n = ph[i]
                n.parentElement.removeChild(n)
              }
    
            } 
            """),
        ]
    )
    
    
    body!(win, myscope, async=false)
    # Blink.AtomShell.opentools(win)
    
    return win
end


function run_electron_js_script(win, script)
    script = replace(script, "@@@" => """windows["$(win.id)"]""")
    js(win.shell, Blink.JSString(script), callback=false)
end

function setup_default_keybinds(win)
    BlinkUtils.run_electron_js_script(win, """
    @@@.webContents.on('before-input-event', (event, input) => {
        if (input.type == 'keyDown') {

            //
            // ctrl-0 - reset zoom factor
            //
            if (input.control && input.key.toLowerCase() === '0') {
                @@@.webContents.setZoomFactor(1.0)
                event.preventDefault()

            //
            // ctrl-+ - zoom in
            //
            } else if (input.control && input.key.toLowerCase() === '=') {
                @@@.webContents.getZoomFactor(function(factor){
                    if (factor < 5) {
                        @@@.webContents.setZoomFactor(factor + 0.2)
                    }
                });
                event.preventDefault()

            //
            // ctrl-- - zoom out
            //
            } else if (input.control && input.key.toLowerCase() === '-') {
                @@@.webContents.getZoomFactor(function(factor){
                    if (factor > 0.4) {
                        @@@.webContents.setZoomFactor(factor - 0.2)
                    }
                });
                event.preventDefault()

            //
            // F12 - toggle dev tools    
            //
            } else if (input.key == 'F12') {
                @@@.webContents.toggleDevTools()
            }
        }
      })    
    """)
end

# shource code is based on https://github.com/tpapp/DefaultApplication.jl/blob/master/src/DefaultApplication.jl
function run_shell_default(command; wait = false)
    @static if Sys.isapple()
        run(`open $(command)`; wait = wait)
    elseif Sys.islinux() || Sys.isbsd()
        run(`xdg-open $(command)`; wait = wait)
    elseif Sys.iswindows()
        cmd = get(ENV, "COMSPEC", "cmd.exe")
        run(`$(cmd) /c start $(command)`; wait = wait)
    else
        @warn("Opening files the default application is not supported on this OS.", KERNEL = Sys.KERNEL)
    end
end

function resolution(html::HTMLParser, var_name::String)
    var = html[var_name]
    resolution = (parse(Int, get(var, "width", 500))-10, parse(Int, get(var, "height", 500))-10)
    return resolution
end


end # module BlinkUtils

# # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# module test

# using Blink, WebIO, Interact, Interact.CSSUtil, DataFrames, CSV, TableView, Widgets
# using ..BlinkUtils


# filename = raw"D:\Projects\Rays\Github\galran\JUMPTest1\AxisymmetricOptimizationTemplate.html"


# p = BlinkUtils.HTMLParser(filename)

# win = BlinkUtils.window(p, title="OpticSim: Axis Simmetric Optimization", size=(1400, 800))


# BlinkUtils.observable(p, "Dummy1")[] = dom"span"("NOT DUMMY")
# BlinkUtils.observable(p, "FilePicker")[] = dom"span"("NOT DUMMY")

# end # module test
