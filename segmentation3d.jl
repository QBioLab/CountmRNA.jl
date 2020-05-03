using FileIO
using Images
using Statistics

"""
Version Comment
v0.1	hf, use yen-threshold select threshold for 3D stack at each time point
v0.2    hf, add multithreads 
"""

function create3dmask(zstack)
    #zstack = imfilter(zstack, Kernel.gaussian((2,2,2)))
	z_depth = 20
    mask = zeros(size(zstack))
    #thresholds_z = [real(yen_threshold(zstack[:, :, i])) for i in 1:20]
	thresholds_z = zeros(Float64, z_depth)
	@inbounds Threads.@threads for z in 1:z_depth
		zstack[:,:,z] = imfilter(zstack[:, :, z], Kernel.gaussian(2))
		thresholds_z[z] = real(otsu_threshold( 
			zstack[81:end-81,81:end-81,z][(zstack[81:end-81, 81:end-81,  z] .> 0.001)] ))
	end
    #threshold_3d = real(yen_threshold(imfilter(zstack, Kernel.gaussian((2,2,1)))))
    #threshold_3d = real(yen_threshold(zstack))
    threshold_3d = maximum(thresholds_z) # choose clearest one
    #mask = opening(zstack .> threshold_3d)
    mask = zstack .> threshold_3d
    remove_small_area!(mask)
    mask, threshold_3d
end

function extract3dnucleus(stack)
    z_depth = 20
    t_len = size(stack)[3] ÷z_depth
    nucleus = zeros(Float64, size(stack))
    thresholds = zeros(Float64, t_len)
    #nucleus_3dmask = zeros(size(stack)[1], size(stack)[2], z_depth)
	println("Extracting nucleus")
    @inbounds Threads.@threads for i in 1:t_len
		z = (i-1)*20+1 : 20*i
        nucleus_3dmask, thresholds[i] = create3dmask(stack[:, :, z])
        nucleus[:, :, z] = nucleus_3dmask .* stack[:, :, z]
    end
    nucleus, thresholds
end

"""
Just remove small regions in binary images
"""
function remove_small_area!(mask)
    mask_con = label_components(mask);
    #mask_res = BitArray(undef, size(mask));
    #mask_res .= false;
    con_size = component_lengths(mask_con)
    selected_con = (1:length(con_size))[con_size .> 5e4]
    """
    for i in 1:maximum(mask_con)
        if sum(mask_con .== i) > 5e3
            mask_res .+= (mask_con.==i);
        end
    end
    """
    @inbounds Threads.@threads for i in 1:length(mask)
        if mask_con[i] ∉ selected_con
            mask[i] = false
        end
    end
    mask
end

#s3c2 = load("../mRNA_confocal_hamamatsu-60X-TIRF/20200316_result/s3_c2.tiff");
#@time nucleus_all, threshold_all = extract3dnucleus(s3c2);
#@time save(File(format"TIFF", "s5-c2_clear.ome.tiff"), N0f16.(nucleus_all))
#x, y, z_all = size(nucleus_all)
#@time embedxml(x, y, 20, z_all÷20, "s5-c2_clear.ome.tiff")
