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

module NotebooksUtils   

# import Pluto
import PlutoUI
# import Markdown
import Format

import Makie
import Makie.AbstractPlotting
import Makie.AbstractPlotting.MakieLayout
import WGLMakie
import GLMakie

mutable struct Defs
    authors::String     # authors
    desc::String        # description of the notebook
    width::Int64
    html::String        # work area to accumulate html syntaxt in order to save cell space (ouput at the end of a single cell)
end
Defs() = Defs("?Authors?")
Defs(authors) = Defs(authors, "")
Defs(authors, desc) = Defs(authors, desc, 1000, "")

function fix_path(path)
    return replace(abspath(path), "\\"=>"/")
end
    

"""
    function run(notebook_filename)

Launch Pluto and allow teh user to open a specific notebook.
"""
function run(notebook_filename)
    run(; path=notebook_filename)
end

"""
    function run(; port=nothing, path=nothing, sysimage_file=nothing )

Launch Pluto and allow teh user to open a specific notebook.
Also allow the usage of a sysimage file for faster loading.        
"""
function run(; port=nothing, path=nothing, sysimage_file=nothing, auto_detect_sysimage = false )

    @eval begin
        import Pluto
        
        local file_path = $path
        local sysimage_file_path = $sysimage_file
        
        if (file_path === nothing)
            @info "Launching Pluto"
        else
            file_path = fix_path(file_path)
            if (!isfile(file_path))
                @warn "Can't find the notebook file [$(file_path)] - trying to locate it in samples"
                base_filename = splitdir(file_path)[2]
                file_path = fix_path(joinpath(splitdir(splitdir(@__DIR__)[1])[1], "samples", "notebooks", base_filename))
                if (!isfile(file_path))
                    @warn "Can't find the notebook file [$(file_path)] in the samples folder - launching Pluto without a file"
                else
                    @info "Launching Notebook [$(file_path)]"
                end
            else
                @info "Launching Notebook [$(file_path)]"
            end
        end

        # handling the sysimage file
        if (sysimage_file_path !== nothing)
            if (!isfile(sysimage_file_path))
                @warn "Can't find the sysimage file [$(sysimage_file_path)] - launching Pluto without a sysimage"
                sysimage_file_path = nothing
            end
        else
            if ($auto_detect_sysimage)
                if Sys.iswindows()
                    if (!isfile("JuliaSysimage.dll"))
                        sysimage_file_path = "JuliaSysimage.dll"
                    end
                elseif Sys.islinux()
                    if (!isfile("JuliaSysimage.so"))
                        sysimage_file_path = "JuliaSysimage.so"
                    end
                end
            end
        end

        local options = Pluto.Configuration.from_flat_kwargs(
            port=$port, 
            notebook=file_path, 
            launch_browser=true,
            host="127.0.0.1",
            require_secret_for_open_links = false,
            require_secret_for_access = false,
            run_notebook_on_load = true,
            sysimage=sysimage_file_path,
            # banner = "yes",
            )
        
        session = Pluto.ServerSession(options=options)
        
        Pluto.run(session)
    end
end

"""
    run_sample(sample_name::String)

Launch Pluto and allow the user to open a specific sample notebook. If a notebook of the same name exists in the current working folder,
it will be opened in Pluto, otherwise, the original sample notebook will be copied to the current folder and be used.
This beheviour will prevent users from updating the original sample notebook. 
"""
function run_sample(sample_name::String, edit_original::Bool = false)
    folder, basename = splitdir(sample_name)
    filename = "DummyFile,jl"
    if isempty(folder)
        if (isfile(basename))
            # file already exists in current folder - launch pluto and point to it
            filename = fix_path(abspath(basename))
        else
            # file does not exists in current directory
            file_in_samples = fix_path(joinpath(splitdir(splitdir(@__DIR__)[1])[1], "samples", "notebooks", basename))
            if (isfile(file_in_samples))
                if (edit_original)
                    @warn "You are about to edit the original sample notebook"
                    filename = file_in_samples
                else
                    # file exists - need to be copied into the current folder
                    src_fn = file_in_samples
                    dst_fn = fix_path(abspath(basename))
                    @info "Copying the sample notebook [$basename] from the samples directory to the current one\n$src_fn ->\n$dst_fn"
                    cp(src_fn, dst_fn)
                    filename = dst_fn
                end
            else
                # file does not exist in samples - need to issue an error
                # @error "Can't find a sample notebook named [file_in_samples] - please check again" 
                throw(ArgumentError("Can't find a sample notebook named [$file_in_samples] - please check again" ))
            end
        end
    end

    @info "Running sample notebook from [$filename]"
    run(; path=filename)
end


mutable struct VariableInfo
    bond::Any
    html::String
end
function GetVarInfo(bond)
    res = VariableInfo(bond, HTMLFromObj(bond))
    return res
end


struct UISlider 
    range::AbstractRange
    default::Number
    dec::Int16
end
UISlider() = UISlider(1:10, 1, -1)
UISlider(r) = UISlider(r, r.start, -1)
UISlider(r, def) = UISlider(r, def, -1)
# UISlider(r, def, format::String) = UISlider(r, def, format)

function Base.show(io::IO, ::MIME"text/html", slider::UISlider)
    if (slider.dec == -1)
        print(io, """
            <input type="range" 
            min="$(first(slider.range))" 
            step="$(step(slider.range))"
            max="$(last(slider.range))" 
            value="$(slider.default)"
            oninput="this.nextElementSibling.value=this.value">
            <output>$(slider.default)</output>""")
    else
        fmt = "%.$(slider.dec)f"
        print(io, """
            <input type="range" 
            min="$(first(slider.range))" 
            step="$(step(slider.range))"
            max="$(last(slider.range))" 
            value="$(slider.default)"
            oninput="this.nextElementSibling.value=parseFloat(this.value).toFixed($(slider.dec))">
            <output>$(Format.cfmt( fmt, slider.default ))</output>""")
    end
end




function SetDefs(defs::Defs)

    @info "I'm in SetDefs"
    ret = PlutoUI.Show(MIME"text/html"(), """
            <style>
            main {
                max-width: $(defs.width)px;
            }
            </style>
            """)

    #@info "The return type $(typeof(ret))"
                
    return ret
end

function DefsClearHTML(defs::Defs)
    defs.html = ""
end

function DefsAddHTML(defs::Defs, html::String)
    defs.html = defs.html * html
end

function DefsHTML(defs::Defs)
    return defs.html
end


# we don't need to use this function anymore - we can extract the html of items and use it using the PlutoUI.show command
# function SetHTMLMarkdown(val::Bool)
#     if (val)
#         @info "SetHTMLMarkdown: Allow Raw HTML Tags in Markdown"
#         @eval Markdown.htmlesc(io::IO, s::AbstractString) = print(io,s)
#     else
#         @info "SetHTMLMarkdown: Original process - Raw HTML is not allowed"
#         @eval Markdown.htmlesc(io::IO, s::AbstractString) = 
#                 for ch in s
#                 print(io, get(Markdown._htmlescape_chars, ch, ch))
#             end	        
#     end
# end

"""
    function SetBackend(defs::Defs, be::String)

    this is my first comment try        
"""
function SetBackend(defs::Defs, be::String)

    if (be == "Web")
        @info "Makie backend set to WEB (WGLMakie)"
        WGLMakie.activate!()
        AbstractPlotting.__init__()
        AbstractPlotting.inline!(true)
    else 
        @info "Makie backend set to STATIC (GLMakie)"
        GLMakie.activate!()
        AbstractPlotting.__init__()
        AbstractPlotting.inline!(true)
    end

end

"""
    InitNotebook(; port=8449)

initialize the JSServe package.
"""
function InitNotebook(; port=8449)
    @eval begin
        try
            import JSServe
            local port = 8449 # the port you want
            JSServe.JSSERVE_CONFIGURATION.listen_port[] = port 
            JSServe.JSSERVE_CONFIGURATION.external_url[] = "http://localhost:$(port)"
            JSServe.JSSERVE_CONFIGURATION.content_delivery_url[] = "http://localhost:$(port)"
            return JSServe.Page() # needs to get displayed by Pluto
        catch e
            @warn "Can't initialize the JSServe package\n$e"
        end
    end
end

function HTMLFromObj(obj)
    io = IOBuffer()
	Base.show(io, MIME"text/html"(), obj)
    res = String(take!(io))
 
    return res
end

# function MDFromString(str)
#     # SetHTMLMarkdown(true)

# 	io = IOBuffer(str)
#     res_md = Markdown.parse(io, )

#     # SetHTMLMarkdown(false)

#     return res_md
# end


function HTMLFix(html::String)
    html = replace(html, "\r" => "")
    html = replace(html, "\n" => "")
    return html
end

function HTMLFloatingBox(items; name="plutoui-genericfloatingbox", header="?? header ??", kwargs...)
    
    res = ""

    res = res * 
          """<nav class="$name aside indent">""" *
          """<header>$header</header>""" *
          """<section>""" * 
          """<span>"""


    for item in items
        item_level = 1
        if     (startswith(item, "@ "))
            item_level = 2
        elseif (startswith(item, "@@ "))
            item_level = 3
        elseif (startswith(item, "@@@ "))
            item_level = 4
        elseif (startswith(item, "@@@@ "))
            item_level = 5
        elseif (startswith(item, "@@@@@ "))
            item_level = 6
        elseif (startswith(item, "@@@@@@ "))
            item_level = 7
        end
        
        if (startswith(item, "@"))
            item2 = item[item_level+2:end]
        else
            item2 = item
        end

        # println("Item [$item]  Level [$item_level]")

        res = res * 
              """<div class="params-row">""" * 
              """<p class = "H$item_level">""" *
              item2 *
              """</p>""" *
              """</div>""" *
              """\n"""

        #Input: $(@bind nnnn MySlider(1:100, 10))
    end

    # closing the nav tags
    res = res *
          """</span>""" *
          """</section>""" *
          """</nav>""" 

    # add the style for this floating box
    res = res * HTMLFloatingBoxStyle(name; kwargs...)

    return HTMLFix(res)
end


function HTMLFloatingBoxStyle(name::String; right="1rem", top="20rem", width="25%", kwargs...)

    @info "HTMLFloatingBoxStyle: right=$right, top=$top, width=$width"

    res = """<style>
    @media screen and (min-width: 1081px) {
        .$name.aside {
            position: fixed;
            right: $right;
            top: $top;
            width: $width;
            padding: 10px;
            border: 3px solid rgba(0, 0, 0, 0.15);
            border-radius: 10px;
            box-shadow: 0 0 11px 0px #00000010;
            max-height: 500px;
            overflow: auto;
            z-index: 5;
            background: white;
        }
    }
    .$name header {
        display: block;
        font-size: 1.5em;
        margin-top: 0.67em;
        margin-bottom: 0.67em;
        margin-left: 0;
        margin-right: 0;
        font-weight: bold;
        border-bottom: 2px solid rgba(0, 0, 0, 0.15);
    }
    .$name section .params-row {
        white-space: nowrap;
        overflow: hidden;
        text-overflow: ellipsis;
        padding-bottom: 2px;
    }
    .highlight-pluto-cell-shoulder {
        background: rgba(0, 0, 0, 0.05);
        background-clip: padding-box;
    }
    .$name section a {
        text-decoration: none;
        font-weight: normal;
        color: gray;
    }
    /* hover */
    .$name section a:hover {
        color: black;
    }
    /* a-ref indentation */
    .$name.indent section a.H1 {
        font-weight: 700;
        line-height: 1em;
    }
    .$name.indent section a.H1 {
        padding-left: 0px;
    }
    .$name.indent section a.H2 {
        padding-left: 10px;
    }
    .$name.indent section a.H3 {
        padding-left: 20px;
    }
    .$name.indent section a.H4 {
        padding-left: 30px;
    }
    .$name.indent section a.H5 {
        padding-left: 40px;
    }
    .$name.indent section a.H6 {
        padding-left: 50px;
    }
    /* paragraph indentation */
    .$name.indent section p.H1 {
        font-weight: 700;
        line-height: 1em;
    }
    .$name.indent section p.H1 {
        padding-left: 0px;
    }
    .$name.indent section p.H2 {
        padding-left: 10px;
    }
    .$name.indent section p.H3 {
        padding-left: 20px;
    }
    .$name.indent section p.H4 {
        padding-left: 30px;
    }
    .$name.indent section p.H5 {
        padding-left: 40px;
    }
    .$name.indent section p.H6 {
        padding-left: 50px;
    }    
</style>"""

    return HTMLFix(res)
end

function HTMLNewDocLayout()

    res = """
    <style>
        body {
            display: block;
        }
        main {
            max-width: 73%;
            padding-left: 50px;
            width: 100%;
        }
    </style>
    """

    return HTMLFix(res)
end

function HTMLFixTOC()

    res = """
    <style>
    @media screen and (min-width: 1081px) {
        .plutoui-toc.aside {
            top: 4%;
            max-height: 40%;
        }
    }
	</style    
    """
    return HTMLFix(res)
end


end # module NotebooksUtils

export NotebooksUtils