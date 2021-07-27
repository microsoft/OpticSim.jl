"""Encapsulates the lens assembly and the display for one lenslet."""

struct Display{T} <: Surface{T}
    surface::Rectangle{T}
    pixellattice::LatticeBasis{3,T}
    pixelpitch::Tuple{T,T} #x,y pitch of pixels
    xresolution::T
    yresolution::T
    transform::Geometry.Transform{T}

    function Display(xres::Int64,yres::Int64,xpitch::T,ypitch::T,transform::Geometry.Transform{T}) where{T<:Real}
        size = (xres*xpitch,yres*ypitch)
        surface = Rectangle(size[1]/T(2),size[2]/T(2),SVector(T(0),T(0),T(1)),SVector(T(0),T(0),T(0))) #surface normal is +z axis
        pixellattice = Repeat.rectangularlattice(size[2],size[1]) #lattice takes ypitch,xpitch in that order
        return new{T}(surface,pixellattice,xres,yres,transform)
    end
end

surfaceintersection(a::Display) = surfaceintersection(a.surface)
normal(a::Display) = normal(a.surface)
interface(a::Display) = interface(a.surface)
makemesh(a::Display) = makemesh(a.surface)
uv(a::Display,p::SVector{3,T}) where{T<:Real} = uv(a.surface,p)
onsurface(a::Display, point::SVector{3,T}) where {T<:Real} = onsurface(a.surface,point)

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
