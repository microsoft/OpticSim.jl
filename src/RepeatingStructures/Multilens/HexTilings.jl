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

function hex3RGB(scale::T = 1.0) where{T<:Real}
    lattice = hex3(scale)
    colors = clustercolors(lattice)
    names = colornames(colors)
    properties =  DataFrames.DataFrame(Color=colors, Name=names)
    return Repeat.ClusterWithProperties{Repeat.RGBCluster}(hex3(scale), properties)
end
export hex3RGB

function hex4RGB(scale::T = 1.0) where{T<:Real}
    lattice = hex4(scale)
    colors = clustercolors(lattice)
    names = colornames(colors)
    properties =  DataFrames.DataFrame(Color=colors, Name=names)
    return Repeat.ClusterWithProperties{Repeat.RGBCluster}(hex4(scale), properties)
end
export hex4RGB

function hex7RGB(scale::T = 1.0) where{T<:Real}
    lattice = hex7(scale)
    colors = clustercolors(lattice)
    names = colornames(colors)
    properties =  DataFrames.DataFrame(Color=colors, Name=names)
    return Repeat.ClusterWithProperties{Repeat.RGBCluster}(hex7(scale), properties)
end
export hex7RGB

function hex9RGB(scale::T = 1.0) where{T<:Real}
    lattice = hex9(scale)
    colors = clustercolors(lattice)
    names = colornames(colors)
    properties =  DataFrames.DataFrame(Color=colors, Name=names)
    return Repeat.ClusterWithProperties{Repeat.RGBCluster}(hex9(scale), properties)
end
export hex9RGB

function hex12RGB(scale::T = 1.0) where{T<:Real}
   
    lattice = hex12(scale)

    colors = clustercolors(lattice)

    names = colornames(colors)
    properties =  DataFrames.DataFrame(Color=colors, Name=names)
    return Repeat.ClusterWithProperties{Repeat.RGBCluster}(lattice, properties)
end
export hex12RGB

function hex18RGB(scale::T = 1.0) where{T<:Real}
    lattice = hex18(scale)
    colors = clustercolors(lattice)
    names = colornames(colors)
    properties =  DataFrames.DataFrame(Color=colors, Name=names)
    # properties =  DataFrames.DataFrame(Color = colors)
    return Repeat.ClusterWithProperties{Repeat.RGBCluster}(lattice, properties)
end
export hex18RGB

function hex19RGB(scale::T = 1.0) where{T<:Real}
    lattice = hex19(scale)

    colors = clustercolors(lattice)
    
    names = colornames(colors)
    properties =  DataFrames.DataFrame(Color=colors, Name=names)
    # properties =  DataFrames.DataFrame(Color = colors)
    return Repeat.ClusterWithProperties{Repeat.RGBCluster}(lattice, properties)
end
export hex19RGB

function hex37RGB(scale::T = 1.0) where{T<:Real}
    lattice = hexn(3,scale)
    colors = clustercolors(lattice)
    names = colornames(colors)
    properties =  DataFrames.DataFrame(Color=colors, Name=names)
    # properties =  DataFrames.DataFrame(Color = colors)
    return Repeat.ClusterWithProperties{Repeat.RGBCluster}(lattice, properties)
end