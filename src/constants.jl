# MIT license
# Copyright (c) Microsoft Corporation. All rights reserved.
# See LICENSE in the project root for full license information.

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
