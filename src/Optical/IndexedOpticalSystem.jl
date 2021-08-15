import Colors
struct LensletProperties
    color::Colors.RGBX
    name::String
end

struct ArraySystem{O<:AbstractOpticalSystem}
    lenssystems::Vector{O}
    properties::Vector{LensletProperties}

    ArraySystem(lenssystems::Vector{O},properties::Vector{LensletProperties}) where{O<:AbstractOpticalSystem} = new{}
end

properties(system::ArraySystem,index::Int) = system.properties[i]

function indexedexample()
    
end
