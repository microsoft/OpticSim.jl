# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

# All of the objects can be displayed with Vis.draw. Example:
# Vis.draw(hex4RGB())

colornames(colors) = uppercase.(x[1]*string((i-1)÷3 + 1) for (i,x) in pairs(colors))
export colornames

function clustercolors(lattice)
    elements = Repeat.clusterelements(lattice)
    colors = Vector{String}(undef, length(elements))

    for (i, element) in pairs(elements)
        colors[i] = pointcolor(element, lattice)
    end
    return colors
end
export clustercolors

function hex3RGB()
    lattice = hex3()
    colors = clustercolors(lattice)
    names = colornames(colors)
    properties =  DataFrames.DataFrame(Color=colors, Name=names)
    return Repeat.ClusterWithProperties(hex3(), properties)
end
export hex3RGB

function hex4RGB()
    lattice = hex4()
    colors = clustercolors(lattice)
    names = colornames(colors)
    properties =  DataFrames.DataFrame(Color=colors, Name=names)
    return Repeat.ClusterWithProperties(hex4(), properties)
end
export hex4RGB

function hex7RGB()
    lattice = hex7()
    colors = clustercolors(lattice)
    names = colornames(colors)
    properties =  DataFrames.DataFrame(Color=colors, Name=names)
    return Repeat.ClusterWithProperties(hex7(), properties)
end
export hex7RGB

function hex9RGB()
    lattice = hex9()
    colors = clustercolors(lattice)
    names = colornames(colors)
    properties =  DataFrames.DataFrame(Color=colors, Name=names)
    return Repeat.ClusterWithProperties(hex9(), properties)
end
export hex9RGB

function hex12RGB()
   
    lattice = hex12()

    colors = clustercolors(lattice)

    names = colornames(colors)
    properties =  DataFrames.DataFrame(Color=colors, Name=names)
    return Repeat.ClusterWithProperties(hex12(), properties)
end
export hex12RGB

function hex19RGB()
    lattice = hex19()

    colors = clustercolors(lattice)
    
    names = colornames(colors)
    properties =  DataFrames.DataFrame(Color=colors, Name=names)
    # properties =  DataFrames.DataFrame(Color = colors)
    return Repeat.ClusterWithProperties(lattice, properties)
end
export hex19RGB