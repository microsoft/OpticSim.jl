

const defaultResolution = (1200,900)

##

mutable struct FigureView 
    figure::Makie.Figure
    name::AbstractString
    saveFileStub::AbstractString
    saveDir::AbstractString
    buttons::Union{Nothing, Vector{Makie.Button}}     #maybe: Union{Nothing, Vector{Button}}
    statusGrid::Any  #grid with status and buttons
    buttonGrid::Any  #buttons for figures
    userGrid::Any    #buttons for user of FigureView
    utilGrid::Any    #buttons for savescreen etc 
    workingGrid::Any #main area, as a grid
    saveRegion::Any  #region to save with SaveScreen
    function FigureView(figure, name, saveFileStub, saveDir)
        new(figure, name, saveFileStub, saveDir, nothing, nothing, nothing, nothing, nothing, nothing, nothing)
    end
end

global numActiveFig = 0
global numNextFigPos = 1
global logo = nothing
global figureViewList = []
global currentFigure = nothing
global currentUserGrid = nothing
global activeFigureView = false


function FigureView(; figureName = string(numNextFigPos), fileNameStub = nothing, dir = nothing,  activeButtonColor = Makie.RGBf0(0.8, 0.94, 0.8), resolution = defaultResolution)
    global logo
    global currentUserGrid
    global numActiveFig
    global currentFigure
    global numNextFigPos
    global activeFigureView
    activeFigureView = true   #direct drawing to managed figures

    if numActiveFig == 0  #first time through
        logo = load("logo.png")
        #Makie.set_window_config!(title = "OpticSim")
    else
        resolution = Makie.size(currentFigure.scene)
    end

    filename = fileNameStub === nothing ? figureName : fileNameStub
    dir = dir===nothing || !isdir(dir) ? pwd() : dir  #if dir isn't valid use current directory. 

    figc = Makie.Figure(resolution = resolution, backgroundcolor = Makie.RGBf0(0.98, 0.98, 0.98), )
    currentFigure = figc
    curFigView = FigureView(figc, figureName, filename, dir)
    push!(figureViewList, curFigView)
    numActiveFig += 1
    numNextFigPos += 1

    if numActiveFig > 0 #Add new buttons to exisiting figures 
        for figView in figureViewList[1:end-1]
            button = figView.buttonGrid[numActiveFig, 1] = Makie.Button(figView.figure, label = figureName)
            push!(figView.buttons, button)
            Makie.on(button.clicks) do n
                currentFigure = figc
                Makie.display(figc)
            end
        end
    end
    #  normally the items would be defined and then included in the constuctor. Some of the Makie objects have unexpected actions 
    #  when being assigned, so assign directly onto the FigureView fields.
    
    figc[1,1] = curFigView.statusGrid = Makie.GridLayout(tellwidth = true, tellheight=false)

    curFigView.statusGrid[1, 1] = curFigView.buttonGrid = Makie.GridLayout(tellwidth = true, tellheight=false)
    curFigView.statusGrid[2, 1] = curFigView.userGrid = Makie.GridLayout(tellwidth = true, tellheight=false)
    currentUserGrid=curFigView.userGrid
    curFigView.statusGrid[3, 1] = curFigView.utilGrid = Makie.GridLayout(tellwidth = true, tellheight=false)   

    figc[1,2] = curFigView.workingGrid = Makie.GridLayout(tellwidth = true, tellheight=true)
    curFigView.saveRegion = figc #workinggrid

    # there is a different "savebutton" on each Figure. inside this method, we define the new one.
    savebutton = curFigView.utilGrid[1,1] = Makie.Button(figc, label = "Save Screen")
    Makie.on(savebutton.clicks) do n
        #println("$(button.label[]) was clicked $n times.")
        saveFigure(dir, filename, figc)
    end

    #display the logo
    makieLogo, = Makie.image(curFigView.utilGrid[2,1], rotr90(logo), axis = (aspect = Makie.DataAspect(),title="OpticSim"))
    Makie.hidedecorations!(makieLogo)  #just the image
    Makie.hidespines!(makieLogo)

    # put figure selection buttons on new Figure

    curFigView.buttons = curFigView.buttonGrid[1:numActiveFig, 1] = [Makie.Button(figc, label = figView.name) for figView in figureViewList]
    curFigView.buttons[end].buttoncolor[] = activeButtonColor

    #link the new buttons to actions
    for (button, fig) in zip(curFigView.buttons,figureViewList)
        Makie.on(button.clicks) do n
            currentFigure = fig.figure
            display(currentFigure)
        end
    end
    Makie.display(figc)
    return curFigView
end

Makie.display(a::FigureView) = Makie.display(a.figure)

function clearFigureView()
    global numActiveFig = 0
    global numNextFigPos = 1
    global logo = nothing
    global figureViewList = []
    global currentFigure = nothing
    global currentUserGrid = nothing
end

function activateFigureView(stat = true)
    global activeFigureView
    activeFigureView = stat
end


global visPlaneX = Plane(SVector(1.0, 0., 0.), SVector(0.0, 0.0,0.0))
global visPlaneY = Plane(SVector(0.0, 1.0, 0.0), SVector(0.0, 0.0, 0.0))

mutable struct DrawingInR
    fV::FigureView
    lscene1::Makie.LScene
    lscene2::Makie.LScene
    lscene3::Makie.LScene
    oneButton::Makie.Button
    twoButton::Makie.Button
    oneswap::Bool
    twoswap::Bool
    xRegion
    yRegion
end

Makie.display(a::DrawingInR) = Makie.display(a.fV)

global currentDrawingInR= nothing

function DrawingInR(; name = "Draw-$numNextFigPos",  xRegion = visPlaneX, yRegion = visPlaneY, fileNameStub = nothing, dir = nothing, resolution = defaultResolution)
    global currentDrawingInR
    fV = FigureView(;figureName = name, fileNameStub, dir, resolution)
    lscene1 = Makie.LScene(fV.workingGrid[1,1], scenekw = (camera = Makie.cam3d_cad!, raw = false))
    lscene2 = Makie.LScene(fV.userGrid[1,1], scenekw = (camera = Makie.cam3d_cad!, raw = false))
    oneButton = Makie.Button(fV.userGrid[2,1], label = "Y")
    lscene3 = Makie.LScene(fV.userGrid[3,1], scenekw = (camera = Makie.cam3d_cad!, raw = false))
    twoButton = Makie.Button(fV.userGrid[4,1], label = "X")

    dInR = DrawingInR(fV, lscene1, lscene2, lscene3, oneButton, twoButton, false, false, xRegion, yRegion)

    Makie.on(dInR.oneButton.clicks) do n
        if dInR.oneswap
            dInR.fV.workingGrid[1,1] = lscene1
            dInR.fV.userGrid[1,1] = lscene2
            dInR.fV.userGrid[3,1] = lscene3
            dInR.oneswap = false
        else
            dInR.fV.workingGrid[1,1] = lscene2
            dInR.fV.userGrid[1,1] = lscene1
            dInR.fV.userGrid[3,1] = lscene3
            dInR.oneswap = true
        end
        dInR.twoswap = false
        Makie.display(dInR)
    end

    Makie.on(dInR.twoButton.clicks) do n
        if dInR.twoswap
            dInR.fV.workingGrid[1,1] = lscene1
            dInR.fV.userGrid[3,1] = lscene3
            dInR.fV.userGrid[1,1] = lscene2
            dInR.twoswap = false
        else
            dInR.fV.workingGrid[1,1] = lscene3
            dInR.fV.userGrid[3,1] = lscene1
            dInR.fV.userGrid[1,1] = lscene2
            dInR.twoswap = true
        end
        dInR.oneswap = false
        Makie.display(dInR)
    end

    currentDrawingInR = dInR
    return dInR
end




