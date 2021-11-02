
using Roots

"""This function returns the radius of the longest basis vector of the lattice cluster. Most lattices defined in this project have symmetric basis vectors so the radii of all basis vectors will be identical."""
function latticediameter(basismatrix::SMatrix)
    maximum(norm.([basismatrix[:,i] for i in 1:size(basismatrix)[2]]))
end
export latticediameter

latticediameter(a::Repeat.AbstractLatticeCluster) =   latticediameter(Repeat.basismatrix(Repeat.clusterbasis(a)))
latticediameter(a::Repeat.AbstractBasis) =  latticediameter(Repeat.basismatrix(a))

struct LensletClusterProperties
    indexofrefraction
    minfnumber
    psfvalue
    λ
    pixelpitch
end

defaultclusterproperties() = LensletClusterProperties(1.5,2.0,.3,530nm,.9μm)
export defaultclusterproperties

"""maximum allowable display for visibility reasons"""
maxlensletdisplaysize() = (150μm,150μm)


const ρ_quartervalue = 2.21509 #value of ρ at which the airy disk function has magnitude .25
export ρ_quartervalue
const ρ_zerovalue = 3.832 #value of ρ at which the airy disk function has magnitude 0

"""given pixelpitch and angular subtense of pixel returns focal length"""
focallength(pixelpitch,θ) = uconvert(mm,.5*pixelpitch/tand(θ/2))
export focallength

"""cutoff frequency for diffraction limited lens is 1/(λ*fnumber) = diameter/(λ*focallength). This has units cycles/distance. For small angles tan(θ) ≈ θ and θ ≈ distance/focallength. distance/rad = focallength. Combine with the first equation to get cutoff frequency in cycles/radian: diameter/(λ*focallength) * focallength = cycles/radian."""
cyclesperdegree(diameter,λ) = uconvert(Unitful.NoUnits,diameter/(rad2deg(1)*λ))
export cyclesperdegree

function mtfcircular(freq,freqcutoff) 
    s =  freq/freqcutoff
    return 2/π * (acos(s) - s*((1-s^2)^.5))
end
export mtfcircular

"""returns the diffraction limit frequency in cycles/deg"""
diffractionlimit(λ,diameter) = uconvert(Unitful.NoUnits,diameter/λ)/rad2deg(1)
export diffractionlimit

function diameter_for_cycles_deg(mtf,cyclesperdeg,λ)
    cyclesperrad = rad2deg(1)*cyclesperdeg
    f(x) = mtfcircular(x,1.0) - mtf
    normalizedfrequency = find_zero(f,(0.0,1.0))
    fc = cyclesperrad/normalizedfrequency
    return uconvert(mm,λ*fc)
end
export diameter_for_cycles_deg

airydisk(ρ) = (2*SpecialFunctions.besselj1(ρ)/ρ)^2

"""The f# required for the first zero of the airy diffraction disk to be at the next sample point"""
diffractionfnumber(λ,pixelpitch,indexofrefraction) = uconvert(Unitful.NoUnits,pixelpitch/(2.44λ/indexofrefraction))
export diffractionfnumber

"""returns ρ, the normalized distance, at which the airy disk will have the value airyvalue"""
function ρatairyvalue(airyvalue)
    ρ = 1e-5 #start at small angle
    while airydisk(ρ) > airyvalue
        ρ += 1e-5 #numerical precision not an issue here
    end
    return ρ
end
export ρatairyvalue

""" Spacing between lenslets which guarantess that for any pixel visible in all lenslets every point on the eyebox plane is covered. This is closest packing of circles."""
closestpackingdistance(pupildiameter) = pupildiameter*cosd(30)
export closestpackingdistance


"""computes the maximum lenslet diameter given the smallest allowed f#. Larger lenslets would require lower f#"""
function maxlensletdiameter(ppd,pixelpitch,minfnumber)
    focallength = pixelpitch/(2*tand(.5/ppd))
    maxdiameter = focallength/minfnumber
    return uconvert(mm,maxdiameter)
end
export maxlensletdiameter

"""Computes the range of lenslet diameters that satisfy two constraints: mtf response at a given frequency, assuming diffraction limited response, and f# is no lower than minfnumber. Returns nothing if constraints can't be satisfied"""
function diameterrange(ppd,pixelpitch,minfnumber,λ,mtf,cyclesperdeg) 
    low = diameter_for_cycles_deg(mtf,cyclesperdeg,λ)
    high = maxlensletdiameter(ppd,pixelpitch,minfnumber)
    if low <= high 
        return (low,high)
    else
        return nothing
    end
end
export diameterrange

function choosecluster(pupildiameter,lensletdiameter)
    clusters = (hex3RGB(),hex4RGB(),hex7RGB(),hex9RGB(),hex12RGB(),hex19RGB())
    cdist = closestpackingdistance(pupildiameter)
    maxcluster = clusters[1]
    ratio = 0.0

    for cluster in clusters
        temp = cdist/(lensletdiameter*latticediameter(cluster))
        if temp >= 1.0
            ratio = temp
            maxcluster = cluster
        else
            break
        end
    end

    @assert ratio ≥ 1.0 "ratio $ratio cdist $cdist lensletdiameter $lensletdiameter latticediameter $(latticediameter(maxcluster)) scaled=$(lensletdiameter*latticediameter(maxcluster))"

    return (cluster = maxcluster,lensletdiameter = lensletdiameter, diameteroflattice = latticediameter(maxcluster)/ratio, packingdistance = cdist*ustrip(mm,lensletdiameter), scalefactor = ratio)
end
export choosecluster

function choosecluster(pupildiameter,ppd,pixelpitch,minfnumber,λ,mtf,cyclesperdeg) 
    diams = diameterrange(ppd,pixelpitch,minfnumber,λ,mtf,cyclesperdeg)
    if diams === nothing
        return nothing
    else
        clusterdata = choosecluster(pupildiameter,diams[1]) #use minimum diameter for now.
    end
    return clusterdata
end

"""Computes display size assuming lenslet normal to eyebox plane passes through lenslet. This is an approximation but for the narrow fov we are considering it is accurate enough to estimate pixel redundandcy, etc."""
function sizeoflensletdisplay(eyerelief::T,eyebox::AbstractVector{T},ppd,pixelpitch::S) where{T<:Unitful.Length,S<:Unitful.Length}
    θ = eyeboxangles(eyebox,eyerelief)
    displaysize = @.  θ*ppd*pixelpitch
    return uconvert.(mm,displaysize)
end
export sizeoflensletdisplay

"""This computes the approximate size of the entire display, not the individual lenslet displays."""
sizeofdisplay(fov,eyerelief) = @. 2*tand(fov/2)*eyerelief
export sizeofdisplay

function numberoflenslets(fov,eyerelief,lensletdiameter)
    lensletarea = π*(lensletdiameter/2)^2
    dispsize = sizeofdisplay(fov,eyerelief)
    return dispsize[1]*dispsize[2]/lensletarea
end
export numberoflenslets

lensletpixels(angles,ppd) = angles .* ppd

"""given the angles each lenslet has to cover compute the corresponding display size"""
sizeoflensletdisplay(angles,ppd,pixelpitch) = @. angles*ppd*pixelpitch

eyeboxangles(eyebox,eyerelief) = @. atand(uconvert(Unitful.NoUnits,eyebox / eyerelief))
export eyeboxangles

function anglesubdivisions(pupildiameter,ppd,pixelpitch,minfnumber,λ,indexofrefraction,psfvalue;RGB = true)
    cluster,_ = choosecluster(pupildiameter,ppd,pixelpitch,minfnumber,λ,indexofrefraction,psfvalue)
    numelements = Repeat.clustersize(cluster)
    if numelements == 19
        return RGB ? (3,2) : (5,3)
    elseif numelements == 12
        return RGB ? (2,2) : (4,3)
   elseif numelements == 9
        return RGB ? (3,1) : (3,3)
    elseif numelements == 7
        return RGB ? (3,1) : (3,2) 
    elseif numelements == 4 || numelements == 3
        return RGB ? (1,1) : (3,1)
    else throw(ErrorException("this should never happen"))
    end
end
export anglesubdivisions


"""computes the approximate fov required of each lenslet for the given constraints. This is strictly correct only for a lenslet centered in front of the eyebox, but the approximation is good enough for high level analysis"""
function lensletangles(eyerelief,eyebox,pupildiameter,ppd; clusterproperties = defaultclusterproperties(), RGB = true)

    return eyeboxangles(eyebox,eyerelief) ./ anglesubdivisions(pupildiameter,ppd,clusterproperties.pixelpitch,clusterproperties.minfnumber,clusterproperties.λ,clusterproperties.indexofrefraction,clusterproperties.psfvalue,RGB = RGB)
end
export lensletangles

testangles() = lensletangles(18mm,(10mm,6mm),4mm,45)
export testangles

lensletresolution(angles,ppd) = angles .* ppd

function pixelredundancy(fov,eyerelief,eyebox,pupildiameter,ppd; RGB = true)
    lensprops = defaultclusterproperties()
    clusterdata = choosecluster(pupildiameter,ppd,lensprops.pixelpitch,lensprops.minfnumber,lensprops.λ,lensprops.indexofrefraction,lensprops.psfvalue)
    nominalresolution = fov .* ppd
    angles = lensletangles(eyerelief,eyebox,pupildiameter,ppd,RGB = RGB)
    pixelsperlenslet = lensletpixels(angles,ppd)
    numlenses = numberoflenslets(fov,eyerelief,clusterdata.lensletdiameter*clusterdata.scalefactor)

    return (numlenses*pixelsperlenslet[1]*pixelsperlenslet[2]) /( nominalresolution[1]*nominalresolution[2])
end
export pixelredundancy

testpixelredundancy() = pixelredundancy((55,35),18mm,(10mm,6mm),4mm,45,RGB = false)
export testpixelredundancy

using Plots
function redundancy_ppdvspupildiameter()
    x = 20:2:45
    y = 3.5:.05:4
   
    plot(Plots.contour(x,y,(x,y) -> pixelredundancy((50,35),18mm,(10mm,6mm),y*mm,x,RGB = false),fill = true))
end
export ppdvspupildiameter

function lensletdisplaysize(fov,eyerelief,eyebox,pupildiameter,ppd; RGB = true)
    lensprops = defaultclusterproperties()
    clusterdata = choosecluster(pupildiameter,ppd,lensprops.pixelpitch,lensprops.minfnumber,lensprops.λ,lensprops.indexofrefraction,lensprops.psfvalue)
    # println("cluster size $(Repeat.clustersize(clusterdata.cluster))")
    nominalresolution = fov .* ppd
    angles = lensletangles(eyerelief,eyebox,pupildiameter,ppd,RGB = RGB)
    pixelsperlenslet = lensletpixels(angles,ppd)
    return @. angles*ppd*lensprops.pixelpitch
end
export lensletdisplaysize

testlensletdisplaysize() = lensletdisplaysize((55,35),18mm,(10mm,6mm),4mm,30,RGB = true)
export testlensletdisplaysize

function displaysize_ppdvspupildiameter()
    x = 20:2:45
    y = 3.5:.05:4
   
    plot(Plots.contour(x,y,(x,y) -> maximum(ustrip.(μm, lensletdisplaysize((50,35),18mm,(10mm,6mm),y*mm,x,RGB = true))),fill = true))
end
export displaysize_ppdvspupildiameter
