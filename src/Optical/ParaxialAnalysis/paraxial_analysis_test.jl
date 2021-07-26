function testvirtualpoint()
    lens = ParaxialLensRect(10.0,5.0,5.0,[0.0,0.0,1.0],[0.0,0.0,0.0])
    displaypoint = [0.0,0.0,-10.0]
    println(virtualpoint(lens,displaypoint))
end
export testvirtualpoint