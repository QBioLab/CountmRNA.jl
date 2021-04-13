using FileIO
using Printf
include("julia2ims.jl")

"""
20201216: Package Spinning disk tif stack to imaris
"""

#src_dir = "/datahub/rawdata/tandeng/mRNA_imaging/mRNA_confocal_again/20201214/sCMOS/HE7-1NLS-uw-PWM_1/Default"
src_dir = "/home/hf/cf2/20201217/HE7-1NLS-uw-PWM-4_1/Default"
out_dir = "/datahub/rawdata/tandeng/mRNA_imaging/mRNA_confocal_again/20201217"
#src_dir = ARGS[1]
t_len = 120
z_len = 20
pos_len = 40

for pos in 0:pos_len-1
#for pos in pos_len-1:pos_len-1
    println(pos)
    buffer = zeros(N0f16, 1900, 1300, z_len, t_len)
    Threads.@threads for t in 0:t_len-1
        id0 = t*z_len*pos_len + pos*z_len
        for z in 0:z_len-1
            if id0+z > 1
                id = @sprintf("img_channel000_position000_time%09d_z000.tif", id0+z-2)
            else # because I miss 2 frame at the begin
                id = @sprintf("img_channel000_position000_time%09d_z000.tif", 0)
            end
            buffer[:, :, z+1, t+1] = load("$src_dir/$id");
        end
    end
    imsname = @sprintf("%s/pos%02d.ims", out_dir, pos)
    save2ims(reinterpret(UInt16, buffer), imsname, dx=0.108, compression=0)
end

