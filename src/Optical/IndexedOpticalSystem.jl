struct IndexedOpticalSystem{O<:AbstractOpticalSystem}
    optics::O
    elements::Dictionary{String,DataFrame{}}
end

function indexedexample()
    
end
