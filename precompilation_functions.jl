using Plots
using OpticSim

#execute functions from packages that normally take a long time to load and precompile

plot(rand(10),rand(10))

Vis.draw(OpticSim.Repeat.Multilens.hex18RGB())
Vis.draw(Examples.cooketriplet())

