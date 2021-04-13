using FileIO
using Images

field = ["s1", "s2", "s9", "s13", "s14", "s16", "s21"];
data_dir = "/datahub/rawdata/tandeng/mRNA_imaging/mRNA_confocal_hamamatsu-60X-TIRF/20200315";
output_dir = "$(data_dir)_visualization"
try
    mkdir("$output_dir")
catch
end

function max_projection(pos)
    println("Processing $pos")
    raw_img = load(File(format"TIFF", "$data_dir/HE7-11-1-80uw-const_1_$pos.ome.btf.tiff"));
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
