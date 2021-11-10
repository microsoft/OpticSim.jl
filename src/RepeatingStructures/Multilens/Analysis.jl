
# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.


using Roots
import DataFrames


# luminance (cd/m2)	Multiple	Value	Item
# 10−6	µcd/m2	1 µcd/m2	Absolute threshold of vision[1]
# 10−5			
# 10−4		400 µcd/m2	Darkest sky[2]
# 10−3	mcd/m2	1 mcd/m2	Night sky[3]
# 1.4   mcd/m2	Typical photographic scene lit by full moon[4]
# 5     mcd/m2	Approximate scotopic/mesopic threshold[5]
# 10−2		40 mcd/m2	Phosphorescent markings on a watch dial after 1 h in the dark[6][7]
# 10−1			
# 100	cd/m2	2 cd/m2	Floodlit buildings, monuments, and fountains[8]
# 5     cd/m2	Approximate mesopic/photopic threshold[5]
# 101		25 cd/m2	Typical photographic scene at sunrise or sunset[4]
# 30    cd/m2	Green electroluminescent source[2]
# 102		250 cd/m2	Peak luminance of a typical LCD monitor[10][11]
# 700   cd/m2	Typical photographic scene on overcast day[4][8][11]
# 103	kcd/m2	2 kcd/m2	Average cloudy sky[2]
# 5     kcd/m2	Typical photographic scene in full sunlight[4][8]
"""
# Pupil diameter as a function of scene luminance
https://jov.arvojournals.org/article.aspx?articleid=2279420
https://en.wikipedia.org/wiki/Orders_of_magnitude_(luminance)

Pupil diameter is approximately 2.8mm at 100cd/m^2. A typical overcast day is 700cd/m^2 
"""

"""computes pupil diameter as a function of scene luminance, in cd/m², and the angular area, a, over which this luminance is presented to the eye."""
𝐃sd(L,a) = 7.75 - 5.75 * ((L * a / 846)^.41) / ((L * a / 846)^.41 + 2) # the first letter of this function name is \bfD not D.


"""This function returns the radius of the longest basis vector of the lattice cluster. Most lattices defined in this project have symmetric basis vectors so the radii of all basis vectors will be identical."""
function latticediameter(basismatrix::SMatrix)
    maximum(norm.([basismatrix[:,i] for i in 1:size(basismatrix)[2]]))
end
export latticediameter

latticediameter(a::Repeat.AbstractLatticeCluster) =   latticediameter(Repeat.basismatrix(Repeat.clusterbasis(a)))
latticediameter(a::Repeat.AbstractBasis) =  latticediameter(Repeat.basismatrix(a))

const hex3latticeclusterbasis = [2 // 1 -1 // 1;-1 // 1 2 // 1]
export hex3latticeclusterbasis

"""returns the integer lattice coords of point in the given basis if the point is in the span of latticebasis. Otherwise returns nothing"""
function latticepoint(latticebasis::AbstractMatrix, origin, point) where {R <: Rational,I <: Integer}
    Ainv = inv(Rational.(latticebasis))
    b = [(point .- origin)...]
    x = Ainv * b
    println(" X $x b $b Ainv $Ainv")
    if reduce(&, (1, 1) .== denominator.(x))
        return Integer.(x)
    else
        return nothing
    end
end
export latticepoint

colorbasis(::Repeat.HexBasis1) = SMatrix{2,2}(2, 0, 0, 2)
colorbasis(::Repeat.HexBasis3) = SMatrix{2,2}(2, -1, 1, 1)
colororigins(::Repeat.HexBasis1) = ((0, 0), (-1, 0), (-1, 1))
colororigins(::Repeat.HexBasis3) = ((0, 0), (0, -1), (1, -1))

"""computes the color associated with a lattice point in the lattice"""
function pointcolor(point, cluster::Repeat.AbstractLatticeCluster) where {T <: Integer}
    latticematrix = colorbasis(Repeat.elementbasis(cluster))
    origins = colororigins(Repeat.elementbasis(cluster))
    colors = zip(origins, ("red", "green", "blue"))
    for (origin, color) in colors
        if nothing !== latticepoint(latticematrix, origin, point)
            return color
        end
    end
end
export pointcolor

""" mtf is the desired response at cycles per degree"""
struct LensletClusterProperties
    mtf
    minfnumber
    cyclesperdegree
    λ
    pixelpitch
end

defaultclusterproperties() = LensletClusterProperties(.2, 2.0, 11, 530nm, .9μm)
export defaultclusterproperties

"""maximum allowable display for visibility reasons"""
maxlensletdisplaysize() = (150μm, 150μm)


const ρ_quartervalue = 2.21509 # value of ρ at which the airy disk function has magnitude .25
export ρ_quartervalue
const ρ_zerovalue = 3.832 # value of ρ at which the airy disk function has magnitude 0

"""given pixelpitch and angular subtense of pixel returns focal length"""
focallength(pixelpitch,θ) = uconvert(mm, .5 * pixelpitch / tand(θ / 2))
export focallength

"""returns the diffraction limit frequency in cycles/deg

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
cyclesperdegree(diameter,λ) = uconvert(Unitful.NoUnits, diameter / (rad2deg(1) * λ))
export cyclesperdegree

"""diffraction limited response for a circular aperture, normalized by maximum cutoff frequency"""
function mtfcircular(freq, freqcutoff) 
    s =  freq / freqcutoff
    return 2 / π * (acos(s) - s * ((1 - s^2)^.5))
end
export mtfcircular


diffractionlimit(λ,diameter) = uconvert(Unitful.NoUnits, diameter / λ) / rad2deg(1)
export diffractionlimit

"""computes the diameter of a diffraction limited lens that has respone mtf at frequency cyclesperdeg"""
function diameter_for_cycles_deg(mtf, cyclesperdeg, λ)
    cyclesperrad = rad2deg(1) * cyclesperdeg
    f(x) = mtfcircular(x, 1.0) - mtf
    normalizedfrequency = find_zero(f, (0.0, 1.0))
    fc = cyclesperrad / normalizedfrequency
    return uconvert(mm, λ * fc)
end
export diameter_for_cycles_deg

airydisk(ρ) = (2 * SpecialFunctions.besselj1(ρ) / ρ)^2

"""The f# required for the first zero of the airy diffraction disk to be at the next sample point"""
diffractionfnumber(λ,pixelpitch,indexofrefraction) = uconvert(Unitful.NoUnits, pixelpitch / (2.44λ / indexofrefraction))
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
closestpackingdistance(pupildiameter) = pupildiameter * cosd(30)
export closestpackingdistance

function choosecluster(pupildiameter, lensletdiameter)
    clusters = (hex3RGB(), hex4RGB(), hex7RGB(), hex9RGB(), hex12RGB(), hex19RGB())
    cdist = closestpackingdistance(pupildiameter)
    maxcluster = clusters[1]
    ratio = 0.0

    for cluster in clusters
        temp = cdist / (lensletdiameter * latticediameter(cluster))
            if temp >= 1.0
            ratio = temp
            maxcluster = cluster
        else
            break
        end
    end

    @assert ratio ≥ 1.0 "ratio $ratio cdist $cdist lensletdiameter $lensletdiameter latticediameter $(latticediameter(maxcluster)) scaled=$(lensletdiameter * latticediameter(maxcluster))"

    return (cluster = maxcluster, lensletdiameter = lensletdiameter, diameteroflattice = latticediameter(maxcluster) / ratio, packingdistance = cdist * ustrip(mm, lensletdiameter), scalefactor = ratio)
end
export choosecluster

function choosecluster(pupildiameter, λ, mtf, cyclesperdeg::T) where {T <: Real} 
    diam = diameter_for_cycles_deg(mtf, cyclesperdeg, λ)
    return choosecluster(pupildiameter, diam) # use minimum diameter for now.
end

"""Computes display size assuming lenslet normal to eyebox plane passes through lenslet. This is an approximation but for the narrow fov we are considering it is accurate enough to estimate pixel redundandcy, etc."""
function sizeoflensletdisplay(eyerelief::T, eyebox::AbstractVector{T}, ppd, pixelpitch::S) where {T <: Unitful.Length,S <: Unitful.Length}
    θ = eyeboxangles(eyebox, eyerelief)
    displaysize = @.  θ * ppd * pixelpitch
    return uconvert.(mm, displaysize)
end
export sizeoflensletdisplay

"""This computes the approximate size of the entire display, not the individual lenslet displays."""
sizeofdisplay(fov,eyerelief) = @. 2 * tand(fov / 2) * eyerelief
export sizeofdisplay

function numberoflenslets(fov, eyerelief, lensletdiameter)
    lensletarea = π * (lensletdiameter / 2)^2
    dispsize = sizeofdisplay(fov, eyerelief)
     return dispsize[1] * dispsize[2] / lensletarea
end
export numberoflenslets

lensletpixels(angles,ppd) = angles .* ppd

"""given the angles each lenslet has to cover compute the corresponding display size"""
sizeoflensletdisplay(angles,ppd,pixelpitch) = @. angles * ppd * pixelpitch

"""angular size of the eyebox when viewed from distance eyerelief"""
eyeboxangles(eyebox,eyerelief) = @. atand(uconvert(Unitful.NoUnits, eyebox / eyerelief))
export eyeboxangles

"""computes how the fov can be subdivided among lenslets based on cluster size"""
function anglesubdivisions(pupildiameter, λ, mtf, cyclesperdegree;RGB=true)
    cluster, _ = choosecluster(pupildiameter, λ, mtf, cyclesperdegree)
    numelements = Repeat.clustersize(cluster)
    if numelements == 19
        return RGB ? (3, 2) : (5, 3)
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


"""computes the approximate fov required of each lenslet for the given constraints. This is strictly correct only for a lenslet centered in front of the eyebox, but the approximation is good enough for high level analysis"""
function lensletangles(eyerelief, eyebox, pupildiameter, ppd; clusterproperties=defaultclusterproperties(), RGB=true)
    cyclesperdegree = ppd / 2.0
    return eyeboxangles(eyebox, eyerelief) ./ anglesubdivisions(pupildiameter, clusterproperties.λ, clusterproperties.mtf, cyclesperdegree, RGB=RGB)
end
export lensletangles

testangles() = lensletangles(18mm, (10mm, 6mm), 4mm, 45)
export testangles

lensletresolution(angles,ppd) = angles .* ppd

"""computes how many more pixels the multilens display will use than a conventional display of the same nominal resolution"""
function pixelredundancy(fov, eyerelief, eyebox, pupildiameter, ppd; RGB=true)
    lensprops = defaultclusterproperties()
    clusterdata = choosecluster(pupildiameter, lensprops.λ, lensprops.mtf, lensprops.cyclesperdegree)
    nominalresolution = fov .* ppd
    angles = lensletangles(eyerelief, eyebox, pupildiameter, ppd, RGB=RGB)
    pixelsperlenslet = lensletpixels(angles, ppd)
    numlenses = numberoflenslets(fov, eyerelief, clusterdata.lensletdiameter * clusterdata.scalefactor)
    return (numlenses * pixelsperlenslet[1] * pixelsperlenslet[2]) / ( nominalresolution[1] * nominalresolution[2])
end
export pixelredundancy

testpixelredundancy() = pixelredundancy((55, 35), 18mm, (10mm, 6mm), 4mm, 45, RGB=false)
export testpixelredundancy

label(color) = color ? "RGB" : "Monochrome"

"""generates a contour plot showing pixel redundancy as a function of ppd and pupil diameter"""
function redundancy_ppdvspupildiameter()
    x = 20:2:45
    y = 3.0:.05:4

    RGB = true
    Plots.plot(Plots.contour(x, y, (x, y) -> pixelredundancy((50, 35), 18mm, (10mm, 6mm), y * mm, x, RGB=RGB), fill=true, xlabel="pixels per degree", ylabel="pupil diameter", legendtitle="pixel redundancy", title="$(label(RGB)) lenslets"))
end
export redundancy_ppdvspupildiameter

"""computes lenslet display size to match the design constraints"""
function lensletdisplaysize(fov, eyerelief, eyebox, pupildiameter, ppd; RGB=true)
    cyclesperdegree = ppd / 2.0
    lensprops = defaultclusterproperties()
    clusterdata = choosecluster(pupildiameter, lensprops.λ, lensprops.mtf, cyclesperdegree)
    nominalresolution = fov .* ppd
    angles = lensletangles(eyerelief, eyebox, pupildiameter, ppd, RGB=RGB)
    pixelsperlenslet = lensletpixels(angles, ppd)
    return @. angles * ppd * lensprops.pixelpitch
end
export lensletdisplaysize

testlensletdisplaysize() = lensletdisplaysize((55, 35), 18mm, (10mm, 6mm), 4mm, 30, RGB=true)
export testlensletdisplaysize

"""generates a contour plot of lenslet display size as a function of ppd and pupil diameter"""
function displaysize_ppdvspupildiameter()
    x = 20:2:45
    y = 3.0:.05:4
    RGB = true

    Plots.plot(Plots.contour(x, y, (x, y) -> maximum(ustrip.(μm, lensletdisplaysize((50, 35), 18mm, (10mm, 6mm), y * mm, x, RGB=RGB))), fill=true, xlabel="pixels per degree", ylabel="pupil diameter", legendtitle="display size μm", title="$(label(RGB)) lenslets"))
end
export displaysize_ppdvspupildiameter
