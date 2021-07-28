"""Encapsulates the lens assembly and the display for one lenslet."""

struct Display{T} <: Surface{T}
    surface::Rectangle{T}
    pixellattice::LatticeBasis{3,T}
    pixelpitch::Tuple{T,T} #x,y pitch of pixels
    xresolution::T
    yresolution::T
    transform::Geometry.Transform{T}

    function Display(xres::Int64,yres::Int64,xpitch::Unitful.Length,ypitch::Unitful.Length,transform::Geometry.Transform{T}) where{T<:Real}
        # convert xpitch,ypitch to unitless numbers representing μm
        pitchx = ustrip(μm,xpitch)
        pitchy = ustrip(μm,ypitch)
        size = (xres*pitchx,yres*pitchy)
        surface = Rectangle(size[1]/T(2),size[2]/T(2),SVector(T(0),T(0),T(1)),SVector(T(0),T(0),T(0))) #surface normal is +z axis
        pixellattice = Repeat.rectangularlattice(size[2],size[1]) #lattice takes ypitch,xpitch in that order
        return new{T}(surface,pixellattice,(pitchx,pitchy),xres,yres,transform)
    end
end

surfaceintersection(a::Display) = surfaceintersection(a.surface)
normal(a::Display) = normal(a.surface)
interface(a::Display) = interface(a.surface)
makemesh(a::Display) = makemesh(a.surface)
uv(a::Display,p::SVector{3,T}) where{T<:Real} = uv(a.surface,p)
onsurface(a::Display, point::SVector{3,T}) where {T<:Real} = onsurface(a.surface,point)
pixelposition(xindex::Int,yindex::Int) = transform*pixellattice[yindex,xindex]
struct LensletAssembly{T}
    lens::ParaxialLens{T}
    transform::Geometry.Transform{T}
    display::Display{T}
    worldtodisplay::Geometry.Transform{T}

    LensletAssembly(lens::ParaxialLens{T},transform::Geometry.Transform{T},display::Display{T}) where{T} = new{T}(lens,transform,display,display.transform*transform)
end

lens(a::LensletAssembly) = a.lens
display(a::LensletAssembly) = a.display
transform(a::LensletAssembly) = a.transform

worldtolens(a::LensletAssembly, pt::SVector{3,T}) where{T<:Real}= a.transform*pt
lenstodisplay(a::LensletAssembly,pt::SVector{3,T}) where{T<:Real} = a.display.transform*pt
worldtodisplay(a::LensletAssembly,pt::SVector{3,T}) where{T<:Real} = a.worldtodisplay*pt

"""Projects vertexpoints in world space onto the lens plane and converts them to two dimensional points represented in the local x,y coordinates of the lens coordinate frame. Used for projecting eye pupil onto lens plane."""
function project(lens::LensletAssembly{T},displaypoint::SVector{3,T},vertexpoints::SVector{N,SVector{3,T}}) where{T<:Real,N}
    #need local transform for lens. Not quite sure how to organize this yet
    projectedpoints = MVector{N,SVector{2,T}}(undef)
    locvertices = map(x->worldtolens(lens,x),vertexpoints) #transform pupil vertices into local coordinate frame of lens
    virtpoint = virtualpoint(lens, displaypoint) #compute virtual point corresponding to physical display point
    for (i,ppoint) in pairs(locvertices)
        vec = ppoint-virtpoint
        vecdist = distancefromplane(lens,ppoint)
        virtdist = distancefromplane(lens,virtpoint)
        scale =  virtdist/vecdist
        planepoint = scale*vec + virtpoint
        projectedpoints[i] = SVector{2,T}(planepoint[1],planepoint[2]) #local lens coordinate frame has z axis aligned with the positive normal to the lens plane.
    end
    return SVector{N,SVector{2,T}}(projectedpoints) #may need to make this a Vector{SVector} to be compatible with LazySets VPolygon constructors. Or maybe an MVector.
end

"""computes the intersection of the lens polygon with a projected polygon, which usually will represent the pupil of the eye. Used in beam energy calculations. Not something that an end user should ever need to access."""
function intersection(lens::ParaxialLens{T},projectedpoints::SVector{N,SVector{2,T}}) where{T<:Real,N}
    lpoly = LazySets.VPolygon(vertices(lens))
    ppoly = LazySets.VPolygon(projectedpoints)
    return lpoly ∩ ppoly
end

"""Returns a number between 0 and 1 representing the ratio of the lens radiance to the pupil radiance. Assume lᵣ is the radiance transmitted through the lens from the display point. Some of this radiance, pᵣ, passes through the pupil. The beam energy is the ratio pᵣ/lᵣ."""
function beamenergy(lens::ParaxialLens{T},displaypoint::AbstractVector{T},pupilpoints::AbstractVector{SVector{3,T}}) where{T<:Real}
    projectedpoints = project(lens,displaypoint,pupilpoints)
    virtpoint = virtualpoint(lens,displaypoint)
    beamlens = SphericalPolygon(vertices(lens),virtpoint,T(1)) #assumes lens vertices are represented in the local lens coordinate frame
    
    intsct = intersection(projectedpoints,vertices(lens))
    if isempty(intsct)
        return T(0)
    else
        beamintsct = SphericalPolygon(intsct,virtpoint,T(1))
        return area(beamintsct)/area(beamlens)
    end
end

"""Compute the bounding box in the display plane of the image of the worldpoly on the display plane. Pixels inside this area, conservatively, need to be turned on because some of their rays may pass through the worldpoly"""
function activebox(lenslet::LensletAssembly{T},worldpoly::SVector{N,SVector{3,T}}) where{T<:Real,N}
    maxx = typemin(T)
    maxy = typemin(T)
    minx = typemax(T)
    miny = typemax(T)
    lens = lens(lenslet)
    display = display(lenslet)

    for point in worldpoly
        locpoint = transform(lenslet)*point

        for lenspoint in vertices(lens)
            ray = Ray(locpoint,lenspoint-locpoint)
            refractedray,_,_ = processintersection(interface(lens),lenspoint,ray, T(OpticSim.GlassCat.TEMP_REF), T(OpticSim.GlassCat.PRESSURE_REF))
            displaypoint = surfaceintersection(display,refractedray)
            maxx = max(maxx,displaypoint[1])
            minx = min(minx,displaypoint[1])
            maxy = max(maxy,displaypoint[2])
            miny = min(miny,displaypoint[2])
        end
    end 
    return minx,miny,maxx,maxy
end

