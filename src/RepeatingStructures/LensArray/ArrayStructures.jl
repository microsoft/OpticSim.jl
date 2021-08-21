module ArrayStructures
export LensArray

using ...OpticSim,...OpticSim.Repeat

struct LensArray
    lenses::Vector{OpticSim.CSGOpticalSystem}
    emitters::Vector{Emitters.Sources.Source}
end

emitters(a::LensArray) = a.emitters
lenses(a::LensArray) = a.lenses

function trace(array::LensArray)
    lensvec, emitvec = lenses(array),emitters(array)

    Threads.@threads for i in 1:length(lensvec)
        println(Threads.threadid())
        lens = lensvec[i]
        emitter = emitvec[i]
        OpticSim.trace(lens,emitter,printprog=false)
    end
    return sum(OpticSim.detectorimage.(lensvec))
end
export trace

moduletest() = println(OpticSim.CSGOpticalSystem)
export moduletest

end #module
export ArrayStructures
