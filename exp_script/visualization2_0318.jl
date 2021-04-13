using FileIO
include("julia2ims.jl")
field = ["s4-2" "s9-1" "s9-2" "s10-5" "s24-2" "s6-1" "s6-2" "s7-5" "s11-1" "s15-1" "s28-2" "s32-1" "s33-1"]
data_dir = "/datahub/rawdata/tandeng/mRNA_imaging/mRNA_confocal_hamamatsu-60X-TIRF/20200318";
output_dir = "$(data_dir)_visualization"

try
    mkdir("$output_dir")
catch
end

function max_projection(cell)
    println("Processing $cell")
    #raw_img = load(File(format"TIFF", "$data_dir/HE7-11-1-80uw-const_3_$pos.ome.tiff"));
    raw_img = loadims2("$(data_dir)_result/output/d18$cell.ims");
    img_size = size(raw_img);
    img_max = zeros(N0f16, img_size[1], img_size[2], Int(img_size[3]/20));
    for i in 1:Int(img_size[3]/20)
        img_max[:, :, i] = maximum(raw_img[:, :, (i-1)*20+1:i*20], dims=3);
    end
    save("$output_dir/$cell.tiff", img_max);
end

for cell in field
    max_projection(cell)
end
