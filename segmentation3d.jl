using FileIO
using Images
using Statistics

"""
Version Comment
v0.1	hf, use yen-threshold select threshold for 3D stack at each time point
v0.2    hf, add multithreads 
"""

function create3dmask(zstack)
    mask = zeros(size(zstack))
    #thresholds_z = [real(yen_threshold(zstack[:, :, i])) for i in 1:20]
    threshold_3d = real(yen_threshold(zstack))
    #threshold_3d = median(thresholds_z)
    mask = opening(zstack .> threshold_3d)
    #mask = zstack .> threshold_3d
    mask, threshold_3d
end

function extract3dnucleus(stack)
    z_depth = 20
    t_len = size(stack)[3] ÷z_depth
    nucleus = zeros(size(stack))
    thresholds = zeros(t_len)
    #nucleus_3dmask = zeros(size(stack)[1], size(stack)[2], z_depth)
	println("Extracting nucleus")
    Threads.@threads for i in 1:t_len
        nucleus_3dmask, thresholds[i] = create3dmask(stack[:, :,(i-1)*20+1:20*i])
        nucleus[:,:,(i-1)*20+1:20*i] = nucleus_3dmask .* stack[:, :,(i-1)*20+1:20*i]
    end
    nucleus, thresholds
end

#s3c2 = load("../mRNA_confocal_hamamatsu-60X-TIRF/20200316_result/s3_c2.tiff");
#@time nucleus_all, threshold_all = extract3dnucleus(s3c2);
#@time save(File(format"TIFF", "s5-c2_clear.ome.tiff"), N0f16.(nucleus_all))
#x, y, z_all = size(nucleus_all)
#@time embedxml(x, y, 20, z_all÷20, "s5-c2_clear.ome.tiff")
