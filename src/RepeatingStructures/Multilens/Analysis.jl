
# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.


using Roots
import DataFrames
import OpticSim.Repeat



#Really should extend AbstractLattice to support cosets, then wouldn't need this ad hoc stuff. TODO
colorbasis(::Repeat.HexBasis1) = SMatrix{2,2}(2, -1, 1, 1)
colorbasis(::Repeat.HexBasis3) = SMatrix{2,2}(2, -1, 1, 1)
colororigins(::Repeat.HexBasis1) = ((0, 0), (-1, 0), (-1, 1))
colororigins(::Repeat.HexBasis3) = ((0, 0), (0, -1), (1, -1))

""" For lenslets arranged in a hexagonal pattern cross talk between lenslets can be reduced by arranging the color of the lenslets so that each lenslet is surrounded only by lenslets of a different color. This function computes the color to assign to any lattice point to ensure this property holds."""
function lensletcolor(point, cluster::Repeat.AbstractLatticeCluster)
    latticematrix = colorbasis(Repeat.elementbasis(cluster))
    origins = colororigins(Repeat.elementbasis(cluster))
    colors = zip(origins, ("red", "green", "blue"))
    for (origin, color) in colors
        if nothing !== Repeat.latticepoint(latticematrix, origin, point)
            return color
        end
    end
end
export lensletcolor

defaultclusterproperties() = (mtf = .2, minfnumber = 2.0, cyclesperdegree = 11, λ = 530nm, pixelpitch = .9μm)
export defaultclusterproperties

"""maximum allowable display for visibility reasons"""
maxlensletdisplaysize() = (250μm, 250μm)


const ρ_quartervalue = 2.21509 # value of ρ at which the airy disk function has magnitude .25
export ρ_quartervalue
const ρ_zerovalue = 3.832 # value of ρ at which the airy disk function has magnitude 0

"""compute focal length given maximum display size, fov, etc."""
function compute_focal_length(fov,lensletdiameter, minfnumber, maxdisplaysize)
    maxangle = deg2rad(max(fov...)) #leaving angles as degrees caused trouble for not obvious reason

    tempfl = uconvert(mm, maxdisplaysize/(2*tan(maxangle/2.0)))
    fnum = tempfl/lensletdiameter
    return tempfl*minfnumber/fnum
end
export compute_focal_length


"""compute focal length required for a given pixel pitch and angular size of the pixel"""
compute_focal_length(pixel_pitch::Unitful.Length,pixel_angle) = .5*pixel_pitch/tan(.5*pixel_angle)

"""minimum lenslet diameter given `ppd`, `fov`, and `pixel_pitch` assuming a Bayer pattern 2x2 pixel. `pixel_pitch` is the distance between individual pixels. The pitch for the Bayer pattern is 2*`pixel_pitch``."""
minimum_lenslet_diameter(ppd,fov::Tuple{T,T},pixel_pitch::Unitful.Length) where{T} = norm(ppd .* fov .* 2 .* pixel_pitch)
export minimum_lenslet_diameter


pixelsperdegree(focal_length,pixelpitch) = 1/(2.0*atand(uconvert(Unitful.NoUnits,pixelpitch/(2.0*focal_length))))
export pixelsperdegree

function lensletdisplaysize(fov,focal_length) 
    uconvert.(μm,2*focal_length .* tand.(fov./2))
end
export lenseletdisplaysize


"""diffraction limited response for a circular aperture, normalized by maximum cutoff frequency"""
function mtfcircular(freq, freqcutoff) 
    s =  freq / freqcutoff
    return 2 / π * (acos(s) - s * ((1 - s^2)^.5))
end
export mtfcircular

""" 
# Returns the diffraction limit frequency in cycles/deg. At this frequeny the response of the system is zero.

## Derivation from the more common formula relating f number and diffraction limit.

focal length = 𝒇𝒍
diffraction cutoff frequency,fc, in cycles/mm = 1/λF# = diameter/λ*𝒇𝒍
cutoff wavelength, Wc, = 1/cutoff frequency = λ*𝒇𝒍/diameter

angular wavelength, θc, radians/cycle, corresponding to cutoff wavelength:

θ ≈ tanθ for small θ
θc corresponding to Wc:  θc ≈ tanθ = Wc/𝒇𝒍
Wc = θc*𝒇𝒍

from equation for Wc:

λ*𝒇𝒍/diameter = Wc = θc*𝒇𝒍
θc = λ/diameter
cycles/rad = 1/θc = diameter/λ
"""
diffractionlimit(λ::Unitful.Length,diameter::Unitful.Length) = uconvert(Unitful.NoUnits, diameter / λ) / rad2deg(1)
export diffractionlimit

"""computes the diameter of a diffraction limited lens that has response mtf at frequency cyclesperdeg"""
function diameter_for_cycles_deg(mtf, cyclesperdeg, λ::Unitful.Length)
    cyclesperrad = rad2deg(1) * cyclesperdeg
    f(x) = mtfcircular(x, 1.0) - mtf
    normalizedfrequency = find_zero(f, (0.0, 1.0))
    fc = cyclesperrad / normalizedfrequency
    return uconvert(mm, λ * fc)
end
export diameter_for_cycles_deg

airydisk(ρ) = (2 * SpecialFunctions.besselj1(ρ) / ρ)^2

"""The f# required for the first zero of the airy diffraction disk to be at the next sample point"""
diffractionfnumber(λ::Unitful.Length,pixelpitch::Unitful.Length,indexofrefraction) = uconvert(Unitful.NoUnits, pixelpitch / (2.44λ / indexofrefraction))
export diffractionfnumber

"""returns ρ, the normalized distance, at which the airy disk will have the value airyvalue"""
function ρatairyvalue(airyvalue)
    ρ = 1e-5 # start at small angle
    while airydisk(ρ) > airyvalue
        ρ += 1e-5 # numerical precision not an issue here
    end
    return ρ
end
export ρatairyvalue

""" Spacing between lenslets which guarantess that for any pixel visible in all lenslets every point on the eyebox plane is covered. This is closest packing of circles."""
closestpackingdistance(pupildiameter::Unitful.Length) = pupildiameter * cosd(30)
export closestpackingdistance

"""Tries clusters of various sizes to choose the largest one which fits within the eye pupil. Larger clusters allow for greater reduction of the fov each lenslet must cover so it returns the largest feasible cluster"""
function choosecluster(pupildiameter::Unitful.Length, lensletdiameter::Unitful.Length, no_eyebox_subdivision::Bool = false)
   #for now don't use cluster sizes that are not divisible by 3 since this causes problems with RGB subdivision of the lenslets. Maybe add back in later.
    # clusters = (hex3RGB, hex4RGB, hex7RGB , hex9RGB, hex12RGB,hex19RGB) #, hex37RGB) # hex37RGB leave out for now. Leads to designs with thousands of small lenslets. May not be practical. 
    if no_eyebox_subdivision
        clusters = (hex3RGB,)
    else   
        clusters = (hex3RGB, hex9RGB, hex12RGB,hex18RGB) #, hex37RGB) # hex37RGB leave out for now. Leads to designs with thousands of small lenslets. May not be practical.
    end
    pupildiameter
    ratio = 0.0
    clusterindex = 0
    scale = 0

    for clusterfunc in clusters
        cluster = clusterfunc() #create an instance of the cluster type with default unit scale
        scale = ustrip(mm,lensletdiameter)/euclideandiameter(elementbasis(cluster)) #scale the element basis of the cluster to match lenslet diameter
        scaledcluster = clusterfunc(scale) #make a new instance of the cluster type scaled so the element basis has diameter equal to lenslet diameter

        temp = ustrip(mm,pupildiameter) / (euclideandiameter(scaledcluster))

        if temp >= 1.0
            ratio = temp
            clusterindex += 1
        else
            break
        end
    end

    maxcluster = clusters[clusterindex](ratio*scale)
    @assert isapprox(Repeat.euclideandiameter(maxcluster), ustrip(mm,pupildiameter))

    scaledlensletdiameter = lensletdiameter*ratio
    return (cluster = maxcluster, lensletdiameter = scaledlensletdiameter, packingdistance = pupildiameter * ustrip(mm, lensletdiameter))
end
export choosecluster

"""Computes the largest feasible cluster size meeting constraints"""
function choosecluster(pupildiameter::Unitful.Length, λ::Unitful.Length, mtf, cyclesperdeg::T, no_eyebox_subdivisions::Bool = false) where {T <: Real} 
    diam = diameter_for_cycles_deg(mtf, cyclesperdeg, λ)
    return choosecluster(pupildiameter, diam, no_eyebox_subdivisions) # use minimum diameter for now.
end

"""This computes the approximate size of the entire display, not the individual lenslet displays."""
sizeofdisplay(fov,eyerelief) = @. 2 * tand(fov / 2) * eyerelief
export sizeofdisplay

"""Computes the number of lenslets required to create a specified field of view at the given eyerelief"""
function numberoflenslets(fov, eyerelief::Unitful.Length, lensletdiameter::Unitful.Length)
    lensletarea = π * (lensletdiameter / 2)^2
    dispsize = sizeofdisplay(fov, eyerelief)
     return dispsize[1] * dispsize[2] / lensletarea
end
export numberoflenslets

"""angular size of the eyebox when viewed from distance eyerelief"""
eyeboxangles(eyebox,eyerelief) = @. atand(uconvert(Unitful.NoUnits, eyebox / eyerelief))
export eyeboxangles

"""This code is tightly linked to the cluster sizes returned by choosecluster. This is not goood design to link these two functions. Think about how to decouple these two functions. computes how the fov can be subdivided among lenslets based on cluster size. Assumes the horizontal size of the eyebox is larger than the vertical so the larger number of subdivisions will always be the first number in the returns Tuple."""
function anglesubdivisions(cluster, pupildiameter::Unitful.Length, λ::Unitful.Length, mtf, cyclesperdegree;RGB=true)
    numelements = Repeat.clustersize(cluster)
    if numelements == 37
        return RGB ? (4,3) : (6,6)
    elseif numelements == 19
        return RGB ? (3, 2) : (5, 3)
    elseif numelements == 18
        return RGB ? (3,2) : (6,3)
    elseif numelements == 12
        return RGB ? (2, 2) : (4, 3)
   elseif numelements == 9
        return RGB ? (3, 1) : (3, 3)
    elseif numelements == 7
        return RGB ? (3, 1) : (3, 2) 
    elseif numelements == 4 || numelements == 3
        return RGB ? (1, 1) : (3, 1)
    else throw(ErrorException("this should never happen"))
    end
end
export anglesubdivisions

"""Multilens displays tradeoff pixel redundancy for a reduction in total track of the display, by using many short focal length lenses to cover the eyebox. This function computes the ratio of pixels in the multilens display vs. a conventional display of the same nominal resolution"""
function pixelredundancy(fov, eyerelief::Unitful.Length, eyebox::NTuple{2,Unitful.Length}, pupildiameter::Unitful.Length, lensletdiameter,angles; RGB=true)
    nominalresolution = ustrip.(°, fov)  #remove degree units so pixel redundancy doesn't have units of °^-2
    numlenses = numberoflenslets(fov, eyerelief, lensletdiameter)
    return (numlenses * angles[1] * angles[2]) / ( nominalresolution[1] * nominalresolution[2])
end
export pixelredundancy

label(color) = color ? "RGB" : "Monochrome"

"""generates a contour plot of lenslet display size as a function of ppd and pupil diameter"""
function displaysize_ppdvspupildiameter()
    x = 20:2:45
    y = 3.0:.05:4
    RGB = true

    Plots.plot(Plots.contour(x, y, (x, y) -> maximum(ustrip.(μm, lensletdisplaysize((50, 35), 18mm, (10mm, 6mm), y * mm, x, RGB=RGB))), fill=true, xlabel="pixels per degree", ylabel="pupil diameter", legendtitle="display size μm", title="$(label(RGB)) lenslets"))
end
export displaysize_ppdvspupildiameter

""" Compute high level properties of a lenslet display system. Uses basic geometric and optical constraints to compute approximate values which should be within 10% or so of a real physical system.

`eyerelief` is in mm: `18mm` not `18`.

`eyebox` is a 2 tuple of lengths, also in mm, that specifies the size of the rectangular eyebox, assumed to lie on vertex of the cornea when the eye is looking directly forward.

`fov` is the field of view of the display as seen by the eye, in degrees (no unit type necessary here).

`pupildiameter` is the diameter of the eye pupil in mm.

`mtf` is the MTF response of the system at the angular frequency `cyclesperdegree`.

`pixelsperdegree` is the number of pixels per degree on the display, as seen by the eye through the lenslets. This is not a specification of the optical resolution of the system (use `mtf` and `cyclesperdegree` for that). It is a physical characteristic of the display.

Example:
```
julia> system_properties(18mm,(10mm,9mm),(55°,45°),4.0mm,.2,11.0)
Dict{Symbol, Any} with 14 entries:
  :lenslet_diameter     => 0.742781 mm
  :lenslet_display_size => (251.708 μm, 345.939 μm)
  :fnumber              => 2.0
  :number_lenslets      => 644.903
  :total_silicon_area   => 56.1555 mm^2
  :subdivisions         => (3, 2)
  :pixel_redundancy     => 33.5192
  :diffraction_limit    => 24.4603
  :eyebox_angles        => (29.0546, 26.5651)
  :pixels_per_degree    => 28.8088
  :lenslet_fov          => (9.68487, 13.2825)
  :display_size         => (18.7404 mm, 14.9117 mm)
  :cluster_data         => (cluster = ClusterWithProperties{19, 2, Float64}(Lat…
  :focal_length         => 1.48556 mm

```
"""
function system_properties(eyerelief::Unitful.Length, eyebox::NTuple{2,Unitful.Length}, fov, pupildiameter::Unitful.Length, mtf, cyclesperdegree; minfnumber=2.0,RGB=true,λ=530nm,pixelpitch=.9μm, maxdisplaysize = 250μm, no_eyebox_subdivision::Bool = false)::Dict{Symbol,Any}
    diameter = diameter_for_cycles_deg(mtf, cyclesperdegree, λ)
    clusterdata = choosecluster(pupildiameter, diameter,no_eyebox_subdivision)
    difflimit = diffractionlimit(λ, clusterdata.lensletdiameter)
    numlenses = numberoflenslets(fov, eyerelief, clusterdata.lensletdiameter)

    subdivisions = anglesubdivisions(clusterdata.cluster, pupildiameter, λ, mtf, cyclesperdegree, RGB=RGB)
    
    eyebox_angles = eyeboxangles(eyebox,eyerelief)
    angles = eyebox_angles ./ subdivisions
    redundancy = pixelredundancy(fov, eyerelief, eyebox, pupildiameter,clusterdata.lensletdiameter, angles,  RGB=RGB)

    focal_length = compute_focal_length(angles,clusterdata.lensletdiameter,minfnumber,maxdisplaysize)
    dispsize = lensletdisplaysize(angles,focal_length)
    fnumber = focal_length/clusterdata.lensletdiameter
    siliconarea = uconvert(mm^2,numlenses  * dispsize[1]*dispsize[2])
    fulldisplaysize = sizeofdisplay(fov,eyerelief)

    pixels_per_degree = pixelsperdegree(focal_length,pixelpitch)
    return Dict(:cluster_data => clusterdata, :lenslet_diameter => clusterdata.lensletdiameter,:pixels_per_degree => pixels_per_degree, :diffraction_limit => difflimit, :fnumber => fnumber, :focal_length => focal_length, :display_size => fulldisplaysize, :lenslet_display_size => dispsize, :total_silicon_area => siliconarea, :number_lenslets => numlenses, :pixel_redundancy => redundancy, :eyebox_angles => eyebox_angles, :lenslet_fov => angles, :subdivisions => subdivisions)
end
export system_properties

"""prints system properties nicely"""
function printsystemproperties(eyerelief::Unitful.Length, eyebox::NTuple{2,Unitful.Length}, fov, pupildiameter::Unitful.Length, mtf, cyclesperdegree; minfnumber=2.0,RGB=true,λ=530nm,pixelpitch=.9μm,maxdisplaysize = 350μm) 
    println("Input parameters")
    println()
    println("pixel pitch $pixelpitch")
    println("eye relief = $eyerelief")
    println("eye box = $eyebox")
    println("fov = $(fov)")
    println("pupil diameter = $pupildiameter")
    println("mtf = $mtf @ $cyclesperdegree cycles/°")
    println()
    printsystemproperties(system_properties(eyerelief, eyebox, fov, pupildiameter, mtf, cyclesperdegree, minfnumber = minfnumber,RGB=RGB,λ=λ,pixelpitch=pixelpitch,maxdisplaysize = maxdisplaysize))
end
export printsystemproperties

function printsystemproperties(props)
    println("Output values")
    println()

    for (key,value) in pairs(props)
        if key == :cluster_data
            println("cluster = $(typeof(props[:cluster_data][:cluster]))")
        elseif key == :diffraction_limit
                println("$key = $value cycles/°")
        else
                println("$key = $value")
        end
    end

    println("occlusion ratio  $(π*uconvert(Unitful.NoUnits,*(props[:lenslet_display_size]...)/(props[:lenslet_diameter])^2))")
    cluster = props[:cluster_data][:cluster]
    clusterdiameter = Repeat.euclideandiameter(cluster)
    println("cluster diameter (approx): $(clusterdiameter)")
end



