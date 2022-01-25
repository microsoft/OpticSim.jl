# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.


# This file contains average properties of the human eye and optical human eye models

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

"""computes pupil diameter as a function of scene luminance `L`, in cd/m², and the angular area, `a`, over which this luminance is presented to the eye."""
𝐃sd(L,a) = 7.75 - 5.75 * ((L * a / 846)^.41) / ((L * a / 846)^.41 + 2) # the first letter of this function name is \bfD not D.
export 𝐃sd

eyeradius() = 12mm
export eyeradius

"""Posterior focal length, i.e., optical distance from entrance pupil to the retina. Focal length will change depending on accomodation. This value is for focus at ∞. When the eye is focused at 25cm focal length will be ≈ 22mm. Because the index of refraction of the vitreous humor is approximately 1.33 the physical distance from the entrance pupil to the retina will be 24mm/1.33 = 18mm."""
eyefocallength() = 24mm
export eyefocallength


vc_epupil() = 3mm #distance from vertex of cornea to entrance pupil

""" distance from vertex of cornea to center of rotation"""
cornea_to_eyecenter() = 13.5mm
export cornea_to_eyecenter

entrancepupil_to_retina() = 22.9mm

"""distance from entrance pupil to center of rotation."""
entrancepupil_to_eyecenter() = entrancepupil_to_retina() - eyeradius()
export entrancepupil_to_eyecenter

"""average angle, in degrees, the eye will rotate before users will turn their head"""
comfortable_eye_rotation_angle() = 20°

"""Average one sided translation of the entrance pupil associated with comfortable eye rotation. If you are using this to define an eyebox multiply this value by 2"""
comfortable_entrance_pupil_translation() = sin(comfortable_eye_rotation_angle())*entrancepupil_to_eyecenter()
export comfortable_entrance_pupil_translation

