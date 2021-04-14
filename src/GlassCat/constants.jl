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

using Unitful

const TEMP_REF = 20.0
const PRESSURE_REF = 1.0
const TEMP_REF_UNITFUL = TEMP_REF * u"°C"
export TEMP_REF, PRESSURE_REF

const DISPFORM_NAMES = ["Schott", "Sellmeier1", "Herzberger", "Sellmeier2", "Conrady", "Sellmeier3", "HandbookOfOptics1", "HandbookOfOptics2", "Sellmeier4", "Extended1", "Sellmeier5", "Extended2", "Extended3"]
const STATUS = ["Standard", "Preferred", "Obsolete", "Special", "Melt"]

# paths for GlassCat source file builds
const GLASSCAT_DIR = @__DIR__ # contains GlassCat.jl (pre-existing)
const AGF_DIR = joinpath(GLASSCAT_DIR, "data", "agf") # contains SCHOTT.agf, Sumita.agf, etc.
const JL_DIR = joinpath(GLASSCAT_DIR, "data", "jl") # contains AGFGlasscat.jl, SCHOTT.jl, etc.

const SOURCES_PATH = joinpath(GLASSCAT_DIR, "data", "sources.txt")
const AGFGLASSCAT_NAME = "AGFGlassCat.jl"
