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
    local mask = zeros(size(zstack))
    #thresholds_z = [real(yen_threshold(zstack[:, :, i])) for i in 1:20]
	local thresholds_z = zeros(Float64, z_depth)
	#@inbounds Threads.@threads for z in 1:z_depth
	@inbounds for z in 1:z_depth
		zstack[:,:,z] = imfilter(zstack[:, :, z], Kernel.gaussian(2))
		local select_pixels = zstack[81:end-81, 81:end-81,  z] .> 0.001
		thresholds_z[z]
		if sum(select_pixels) > 2e2
			thresholds_z[z] = real(otsu_threshold( 
				zstack[81:end-81,81:end-81,z][select_pixels]))
		end
	end
    #threshold_3d = real(yen_threshold(imfilter(zstack, Kernel.gaussian((2,2,1)))))
    #threshold_3d = real(yen_threshold(zstack))
    threshold_3d = maximum(thresholds_z) # choose clearest one
    #mask = opening(zstack .> threshold_3d)
    mask = zstack .> threshold_3d
    #remove_small_area!(mask)
    mask, mask_size = remove_small_area(mask)
    mask, mask_size, threshold_3d
end

function extract3dnucleus(stack)
    z_depth = 20
    t_len = size(stack)[3] ÷z_depth
    local nucleus = zeros(Gray{Normed{UInt16,16}}, size(stack))
	local nucleus_len = length(stack)
    local thresholds = zeros(Float64, t_len)
    local mask_size = zeros(Float64, t_len)
    #nucleus_3dmask = zeros(size(stack)[1], size(stack)[2], z_depth)
	println("Extracting nucleus")
    @inbounds Threads.@threads for t in 1:t_len
    #@inbounds for i in 1:t_len
		local z = (t-1)*20+1 : 20*t
		nucleus_3dmask, mask_size[t], thresholds[t] = create3dmask(stack[:, :, z])
		nucleus_t = view(nucleus, :,:,z)
		stack_t = view(stack, :,:,z)
        #nucleus[:, :, z] = nucleus_3dmask .* stack[:, :, z]
		#for i in 1:nucleus_len
		for i in eachindex(nucleus_3dmask)
			if nucleus_3dmask[i]
				nucleus_t[i] = stack_t[i]
			end
		end
    end
    nucleus, mask_size, thresholds
end

"""
Just remove small regions in binary images
"""
#function remove_small_area!(mask)
function remove_small_area(mask)
    mask_con = label_components(mask);
    #mask_res = BitArray(undef, size(mask));
    #mask_res .= false;
    con_size = component_lengths(mask_con)
    selected_con = (0:length(con_size)-1)[con_size .> 2e4]
    #@inbounds Threads.@threads for i in 1:length(mask)
	local mask_size = 0
    @inbounds for i in 1:length(mask)
        if mask_con[i] ∉ selected_con
            mask[i] = false
		else
			mask_size += 1
        end
    end
    mask, mask_size
end

"""
Computes threshold for grayscale image using Otsu's method. (modify from Imagas.jl)
Parameters:
-    img         = Grayscale input image
-    bins        = Number of bins used to compute the histogram. Needed for floating-point images.
"""
function otsu_threshold2(img::AbstractArray{T, N}, min, bins::Int = 256) where {T<:Union{Gray,Real}, N}

    #min, max = extrema(img)
    max = maximum(img)
	if max == 0 # avoid image without cell
		return T(0)
	end
    edges, counts = imhist(img, range(gray(min), stop=gray(max), length=bins))
	if sum(counts)<2e2 # don't count small region
		return T(0)
	end
    histogram = counts./sum(counts)

    ω0 = 0
    μ0 = 0
    μt = 0
    μT = sum((1:(bins+1)).*histogram)
    max_σb=0.0
    thres=1

    for t in 1:bins
        ω0 += histogram[t]
        ω1 = 1 - ω0
        μt += t*histogram[t]

        σb = (μT*ω0-μt)^2/(ω0*ω1)

        if(σb > max_σb)
            max_σb = σb
            thres = t
        end
    end

    return T((edges[thres-1]+edges[thres])/2)
end

#s3c2 = load("../mRNA_confocal_hamamatsu-60X-TIRF/20200316_result/s3_c2.tiff");
#@time nucleus_all, threshold_all = extract3dnucleus(s3c2);
#@time save(File(format"TIFF", "s5-c2_clear.ome.tiff"), N0f16.(nucleus_all))
#x, y, z_all = size(nucleus_all)
#@time embedxml(x, y, 20, z_all÷20, "s5-c2_clear.ome.tiff")
