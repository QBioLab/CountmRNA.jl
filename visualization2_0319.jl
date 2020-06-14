using FileIO
include("julia2ims.jl")
field = ["s1-1" "s11-3" "s12-1" "s23-2" "s27-3" "s29-1" "s32-2" "s34-4" "s36-1" "s3-1" "s11-2" "s23-3" "s27-2" "s32-1" "s33-3" ]
data_dir = "/datahub/rawdata/tandeng/mRNA_imaging/mRNA_confocal_hamamatsu-60X-TIRF/20200319";
output_dir = "$(data_dir)_visualization"

try
    mkdir("$output_dir")
catch
end

function max_projection(cell)
    println("Processing $cell")
    #raw_img = load(File(format"TIFF", "$data_dir/HE7-11-1-80uw-const_3_$pos.ome.tiff"));
    raw_img = loadims2("$(data_dir)_result/output/d19$cell.ims");
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
