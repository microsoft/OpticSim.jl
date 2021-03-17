# MIT License

# Copyright (c) Microsoft Corporation.

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE

# small offset along ray direction to prevent ray from reintersecting the surface it just intersected. Optical pathlength is measured in mm so this is much smaller than a wavelength of light.
const RAY_OFFSET = 1e-6
# if the power of a ray drops below this value return nothing and do not continue the ray tracing recursion
const POWER_THRESHOLD = 1e-6
# if the ray is traced recursively more than this many times then the trace will be terminated and the ray ignored (mostly happens with TIR)
const TRACE_RECURSION_LIMIT = 256

# triangulated surfaces are imprecise mostly due to floating point precision, particularly if using trig functions
# expanding every triangle by a fraction resolves missed intersections
const TRIANGULATION_EXPANSION = 2 * eps(Float64)
# we also need to extend the bounds of the surface to ensure we don't miss anythgin at the edges of the surface
const TRIANGULATION_EXTENSION = 0.02

# for accelerated surface the triangulated surface can often be very imprecise, this will mean rays which start close
# to the surface might miss intersection because the ray origin lies above the triangulated surface but below the true surface
# (where for this example the ray direction is straight upwards), offsetting the ray fixes this, and as we check the alpha
# of the intersection point afterwards there are no problems
const ACCEL_SURF_RAY_OFFSET = 0.1

# settings for visualization
const MIN_VIS_TRI_SIZE = 1e-3
const MESH_PRECISION = 1e-3
const VIS_RECURSION_LIMIT = 64

# number of M/N to precompute for various variables involved in the QType calculations
const QTYPE_PRECOMP = 16

# maximum number of orders for a thin grating interface
const GRATING_MAX_ORDERS = 10

# LensAssemblyN methods will be generated at compile time (faster on first run) for N up to this number
const PREGENERATED_LENS_ASSEMBLY_SIZE = 10
