using FileIO
using Printf
include("julia2ims.jl")

"""
20201216: Package Spinning disk tif stack to imaris
"""

src_dir = "/datahub/rawdata/tandeng/mRNA_imaging/mRNA_confocal_again/20201215/sCMOS"
out_dir = "/datahub/rawdata/tandeng/mRNA_imaging/mRNA_confocal_again/20201215"
#src_dir = ARGS[1]
t_len = 42 # 79-120 time point
z_len = 20
pos_len = 40

buffer = zeros(N0f16, 1900, 1300, z_len*t_len*pos_len)
img_num = 866

raw = load(File(format"TIFF", "$src_dir/HE7-1NLS-uw-PWM_1_MMStack_72.ome.tif"))
buffer[:, :, 1:818] = raw[:, :, 49:end] # push 364 slice

Threads.@threads for index in 73:109 # 18 file here, remain 412 slice
    println(index)
    filename = "HE7-1NLS-uw-PWM_1_MMStack_$index.ome.tif"
    raw = load(File(format"TIFF", "$src_dir/$filename"))
	GC.gc()
    id = (index-73)*img_num + 1 + 818
    buffer[:, :, id:id+img_num-1] = raw
end

raw = load(File(format"TIFF", "$src_dir/HE7-1NLS-uw-PWM_1_MMStack_110.ome.tif"))
buffer[:, :, end-740+1:end] = raw[:, :, :] # push 48 slice

buffer = reshape(buffer, (1900, 1300, z_len, t_len*pos_len))
Threads.@threads for pos in 1:pos_len
    imsname = @sprintf("%s/pos%02d-4.tiff", out_dir, pos)
    #save2ims(reinterpret(UInt16, buffer[:, :, :, :, pos]),
    #         imsname, dx=0.108, compression=1)
    save(imsname, reshape(buffer[:, :, :, pos:40:end],
                          (1900, 1300, z_len*t_len)))
end
