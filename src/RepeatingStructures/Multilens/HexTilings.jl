# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

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
    names = ["R","G","B"]
    properties =  DataFrames.DataFrame(Color=colors, Name=names)
    return Repeat.ClusterWithProperties(hex3(), properties)
end
export hex3RGB

function hex4RGB()
    lattice = hex4()
    colors = clustercolors(lattice)
    names = ["R","G","B","W"]
    properties =  DataFrames.DataFrame(Color=colors, Name=names)
    return Repeat.ClusterWithProperties(hex4(), properties)
end
export hex4RGB

function hex7RGB()
    lattice = hex7()
    colors = clustercolors(lattice)
    names = ["R1","G1","B1","W","R2","G2","B2"]
    properties =  DataFrames.DataFrame(Color=colors, Name=names)
    return Repeat.ClusterWithProperties(hex7(), properties)
end
export hex7RGB

function hex9RGB()
    lattice = hex9()
    colors = clustercolors(lattice)
    names = ["R-1","G-1","B-1","R0","G0","B0","R1","G1","B1"]
    properties =  DataFrames.DataFrame(Color=colors, Name=names)
    return Repeat.ClusterWithProperties(hex9(), properties)
end
export hex9RGB

function hex12RGB()
   
    lattice = hex12()

    colors = clustercolors(lattice)

    names = [
        "G3","B0",
        "B2","R1","G2",
        "R0","G1","B1","R3",
        "B3","R2","G0"
    ]
    properties =  DataFrames.DataFrame(Color=colors, Name=names)
    return Repeat.ClusterWithProperties(hex12(), properties)
end
export hex12RGB

function hex19RGB()
    lattice = hex19()

    colors = clustercolors(lattice)
    
    names = [
        "W",
        "B0","G1","R2","B3","G4","R5",
        "G0","B1","R1","G2","B2","R3","G3","B4","R4","G5","B5","R0"
    ]
    properties =  DataFrames.DataFrame(Color=colors, Name=names)
    # properties =  DataFrames.DataFrame(Color = colors)
    return Repeat.ClusterWithProperties(lattice, properties)
end
export hex19RGB

function hex19fov1()
    lattice = hex19()

    offset = 1
    colornames = ["R","G","B"]
    colors = ["red","green","blue"]
    fov1 = [8,19,18,17,7,2]
    fov2 = [16,15,14,13,5,6] 
    fov3 = [12,11,10,9,3,4]

    allfov = vcat([0],fov1,fov2,fov3)
    names = [x -> string("fov",x) for x in 1:6]
    
    allnames = Vector{String}(undef,19)
    allcolors = Vector{String}(undef,19)

    for (index,fov) in pairs((fov1,fov2,fov3))
        for i in eachindex(fov)
            println("$index $i")
            latticeindex = (index-1)*length(fov1) + i + offset
            allcolors[allfov[latticeindex]] = colors[index]
            allnames[allfov[latticeindex]] = string("fov",string(index))
        end
    end
    allcolors[offset] = "white"
    allnames[offset] = "W"
    properties =  DataFrames.DataFrame(Color=allcolors, Name=allnames)
    # properties =  DataFrames.DataFrame(Color = colors)
    return Repeat.ClusterWithProperties(lattice, properties)
end
export hex19fov1


function hex19fov2()
    lattice = hex19()

    colors = clustercolors(lattice)
    names = [string("fov",x÷3) for x in 1:18]
    names = vcat("W",names)
    properties =  DataFrames.DataFrame(Color=colors, Name=names)
    # properties =  DataFrames.DataFrame(Color = colors)
    return Repeat.ClusterWithProperties(lattice, properties)
end
export hex19fov2