using OpticSim
using OpticSim.Emitters
using OpticSim.Geometry

using Unitful: mm
using StaticArrays
using FileIO
using Images

focaldistance = 3.
eye_radius = 12.

image = Gray.(load("test.png"))
pixel_size = (.1, .1)
pixel_pitch = (.15, .15)
pixel_pos = (0., 0.)
transform = translation(2., 2., 1 / (1/focaldistance - .5/eye_radius)) * rotationY(π) * rotationZ(π/2)
raygenerator = Emitters.imageemitter(image, pixel_size, pixel_pitch, pixel_pos; transform)

detpixels = 400
eye = OpticSim.ParaxialEye(; focaldistance, eye_radius, detpixels)

trackallrays = test = colorbysourcenum = verbose = true; resolution = (852, 480)
Vis.drawtracerays(eye; raygenerator, trackallrays, test, colorbysourcenum, resolution, verbose)
Vis.drawtraceimage(eye; raygenerator, test, verbose);
