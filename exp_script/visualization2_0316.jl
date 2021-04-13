using FileIO
include("julia2ims.jl")

field = ["s1-6"  "s3-4"  "s7-12" "s8-3"  "s8-4"  "s12-1" "s17-9" "s14-7" "s14-5" "s15-11"    "s19-9" "s20-7" "s21-1" "s26-2" "s31-1" "s33-8" "s37-7" "s37-8" "s30-7" "s40-3" ]
data_dir = "/datahub/rawdata/tandeng/mRNA_imaging/mRNA_confocal_hamamatsu-60X-TIRF/20200316";
output_dir = "$(data_dir)_visualization"

try
    mkdir("$output_dir")
catch
end

function max_projection(cell)
    println("Processing $cell")
    #raw_img = load(File(format"TIFF", "$data_dir/HE7-11-1-80uw-const_3_$pos.ome.tiff"));
    raw_img = loadims2("$(data_dir)_result/output/d16$cell.ims");
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
