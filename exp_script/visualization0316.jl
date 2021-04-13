
using FileIO
using Images

field = [ "s3", "s8", "s14", "s15", "s20","s37", "s40"]
data_dir = "/datahub/rawdata/tandeng/mRNA_imaging/mRNA_confocal_hamamatsu-60X-TIRF/20200316";
output_dir = "$(data_dir)_visualization"
try
    mkdir("$output_dir")
catch
end

function max_projection(pos)
    println("Processing $pos")
    raw_img = load(File(format"TIFF", "$data_dir/HE7-11-1-80uw-PWM_1_$pos.ome.tiff"));
    img_size = size(raw_img);
    img_max = zeros(N0f16, img_size[1], img_size[2], Int(img_size[3]/20));
    for i in 1:Int(img_size[3]/20)
        img_max[:, :, i] = maximum(raw_img[:, :, (i-1)*20+1:i*20], dims=3);
    end
    save("$output_dir/$pos.tiff", img_max);
end

for pos in field
    max_projection(pos)
end
