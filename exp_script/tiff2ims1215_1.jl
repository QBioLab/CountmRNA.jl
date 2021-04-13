using TIFF
using FileIO
using Printf
include("julia2ims.jl")

"""
20201216: Package Spinning disk tif stack to imaris
"""

src_dir = "/datahub/rawdata/tandeng/mRNA_imaging/mRNA_confocal_again/20201215/sCMOS"
out_dir = "/datahub/rawdata/tandeng/mRNA_imaging/mRNA_confocal_again/20201215"
t_len = 28 #1-28
z_len = 20
pos_len = 40

buffer = zeros(N0f16, 1900, 1300, z_len*t_len*pos_len)
img_num = 866

Threads.@threads for index in 0:24
    println(index)
    filename = "HE7-1NLS-uw-PWM_1_MMStack_$index.ome.tif"
    raw = TIFF.load("$src_dir/$filename") 
    id = index*img_num + 1
    buffer[:, :, id:id+img_num-1] = raw.data
    GC.gc()
end

#Threads.@threads for index in 25:54
#    filename = "HE7-1NLS-uw-PWM_1_MMStack_$index.ome.tif"
#    raw = load(File(format"TIFF", "$src_dir/$filename")) 
#    id = index*img_num + 1
#    buffer[:, :, id:id+img_num-1] = raw.data
#    GC.gc()
#end

t0 = 25*866
remain = t_len*z_len*pos_len - t0 # 750
raw = load(File(format"TIFF", "$src_dir/HE7-1NLS-uw-PWM_1_MMStack_25.ome.tif"))
buffer[:, :, t0+1:end] = raw[:, :, 1:remain]

buffer = reshape(buffer, (1900, 1300, z_len, t_len*pos_len))
Threads.@threads for pos in 1:pos_len
    imsname = @sprintf("%s/pos%02d-1.tiff", out_dir, pos)
    #save2ims(reinterpret(UInt16, buffer[:, :, :, :, pos]),
    #         imsname, dx=0.108, compression=1)
    save(imsname, reshape(buffer[:, :, :, pos:40:end], 
        (1900, 1300, z_len*t_len)))
end
