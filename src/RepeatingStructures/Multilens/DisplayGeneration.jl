#display lenslet systems: emitter, paraxial lenslet, paraxial eye lens, retina detector.
#start with eyebox plane and hexagonal tiling. Project onto display surface, generate lenslets automatically.

function project(vertices::SMatrix{3,N,T},projectionvector::AbstractVector{T},surface::OpticSim.Surface{T}) where{N,T}
    for i in size(vertices)[2]
        origin = vertices[:,i]
        ray = OpticSim.Ray(origin,projectionvector)
        intsct = OpticSim.surfaceintersection(surface,ray)
        println(intsct)
    end
end

function testproject()
    normal = [0,0,1]
    hex = Repeat.HexBasis1()
    verts = vertices(hex)
    surf = OpticSim.Plane(normal,[0,0,10],)
end
