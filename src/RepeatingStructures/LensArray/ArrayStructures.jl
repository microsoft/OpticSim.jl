module ArrayStructures

using ...OpticSim

struct LensArray
    lenses::Vector{CSGOpticalSystem}
    emitters::Vector{Emitters.Source}
end

emitters(a::LensArray) = a.emitters
lenses(a::LensArray) = a.lenses

function trace(array::LensArray) 
    Threads.@thread for (lens,emitter) in zip(lenses(array),emitters(array))
        trace(lens,emitter,printprog=false)
    end
    return sum(detectorimage.(lenses))
end

end #module
