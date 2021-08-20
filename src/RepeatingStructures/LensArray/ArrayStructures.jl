module ArrayStructures
export LensArray

using ...OpticSim

struct LensArray
    lenses::Vector{OpticSim.CSGOpticalSystem}
    emitters::Vector{Emitters.Sources.Source}
end

emitters(a::LensArray) = a.emitters
lenses(a::LensArray) = a.lenses

function trace(array::LensArray)
    lesvec, emitvec = lenses(array),emitters(array)
    Threads.@threads for i in 1:length(lensvec)
        lens = lensvec[i]
        emitter = emitvec[i]
        trace(lens,emitter,printprog=false)
    end
    return sum(detectorimage.(lenses))
end
export trace

moduletest() = println(CSGOpticalSystem)
export moduletest

end #module
export ArrayStructures
