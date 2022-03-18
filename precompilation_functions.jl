using Plots
using OpticSim

plot(rand(10),rand(10))

Vis.draw(OpticSim.Repeat.Multilens.hex16RGB())

