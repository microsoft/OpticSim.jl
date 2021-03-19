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

# all other glasses should follow the format below, new glasses must be added to OTHER_GLASSES and OTHER_GLASS_NAMES where the index in the array matches the numeric part of the GlassID

module CARGILLE
using ..GlassCat: Glass, GlassID, OTHER
using StaticArrays: SVector

"""    CARGILLE.OG0608
```
ID:                      OTHER:1
RI @ 587nm:              1.457518
Abbe Number:             57.18978
ΔPgF:                    0.008
TCE (÷1e-6):             800.0
Density:                 0.878g/m³
Valid wavelengths:       0.32μm to 1.55μm
Reference Temp:          25°C
```
"""
const OG0608 = Glass(GlassID(OTHER, 1), -2, 1.4451400, 0.0043176, -1.80659e-5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.32, 1.55, -0.0009083144750540808, 0.0, 0.0, 0.0, 0.0, 0.0, 25.0, 0.008, -1.0, -1.0, 800.0, -1.0, 0, -1.0, [SVector(0.32, 0.03, 10.0), SVector(0.365, 0.16, 100.0), SVector(0.4047, 0.40, 100.0), SVector(0.480, 0.71, 100.0), SVector(0.4861, 0.72, 100.0), SVector(0.5461, 0.80, 100.0), SVector(0.5893, 0.90, 100.0), SVector(0.6328, 0.92, 100.0), SVector(0.6439, 0.95, 100.0), SVector(0.6563, 0.96, 100.0), SVector(0.6943, 0.99, 100.0), SVector(0.840, 0.99, 100.0), SVector(0.10648, 0.74, 100.0), SVector(0.1300, 0.39, 100.0), SVector(0.1550, 0.16, 100.0)], 1.457518, -1.0, -1.0, 0, 57.18978, 0, 0.878, -1)

"""    CARGILLE.OG0607
```
ID:                      OTHER:2
RI @ 587nm:              1.457587
Abbe Number:             57.19833
ΔPgF:                    0.008
TCE (÷1e-6):             700.0
Density:                 0.878g/m³
Valid wavelengths:       0.32μm to 1.55μm
Reference Temp:          25°C
```
"""
const OG0607 = Glass(GlassID(OTHER, 2), -2, 1.44503, 0.0044096, -2.85878e-5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.32, 1.55, -0.0009083144750540808, 0.0, 0.0, 0.0, 0.0, 0.0, 25.0, 0.008, -1.0, -1.0, 700.0, -1.0, 0, -1.0, [SVector(0.32, 0.15, 10.0), SVector(0.365, 0.12, 100.0), SVector(0.4047, 0.42, 100.0), SVector(0.480, 0.78, 100.0), SVector(0.4861, 0.79, 100.0), SVector(0.5461, 0.86, 100.0), SVector(0.5893, 0.90, 100.0), SVector(0.6328, 0.92, 100.0), SVector(0.6439, 0.90, 100.0), SVector(0.6563, 0.92, 100.0), SVector(0.6943, 0.98, 100.0), SVector(0.840, 0.99, 100.0), SVector(0.10648, 0.61, 100.0), SVector(0.1300, 0.39, 100.0), SVector(0.1550, 0.11, 100.0)], 1.457587, -1.0, -1.0, 0, 57.19833, 0, 0.878, -1)

"""    CARGILLE.OG081160
```
ID:                      OTHER:3
RI @ 587nm:              1.515549
Abbe Number:             36.82493
ΔPgF:                    0.014
TCE (÷1e-6):             700.0
Density:                 1.11g/m³
Valid wavelengths:       0.32μm to 1.55μm
Reference Temp:          25°C
```
"""
const OG081160 = Glass(GlassID(OTHER, 3), -2, 1.49614, 0.00692199, -8.07052e-5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.32, 1.55, -0.000885983052189022, 0.0, 0.0, 0.0, 0.0, 0.0, 25.0, 0.014, -1.0, -1.0, 700.0, -1.0, 0, -1.0, [SVector(0.32, 0.04, 100.0), SVector(0.365, 0.13, 100.0), SVector(0.4047, 0.26, 100.0), SVector(0.480, 0.48, 100.0), SVector(0.4861, 0.49, 100.0), SVector(0.5461, 0.60, 100.0), SVector(0.5893, 0.68, 100.0), SVector(0.6328, 0.71, 100.0), SVector(0.6439, 0.73, 100.0), SVector(0.6563, 0.74, 100.0), SVector(0.6943, 0.76, 100.0), SVector(0.840, 0.83, 100.0), SVector(0.10648, 0.86, 100.0), SVector(0.1300, 0.89, 100.0), SVector(0.1550, 0.90, 100.0)], 1.515549, -1.0, -1.0, 0, 36.82493, 0, 1.11, -1)

end

const OTHER_GLASSES = [CARGILLE.OG0607, CARGILLE.OG0608, CARGILLE.OG081160]
const OTHER_GLASS_NAMES = ["CARGILLE.OG0607", "CARGILLE.OG0608", "CARGILLE.OG081160"]
