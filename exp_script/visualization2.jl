using FileIO
include("julia2ims.jl")

field = ["s3-2" "s3-5" "s3-2" "s3-5" "s3-6"	"s4-9"	"s4-6"	"s4-4"	"s5-4"	"s5-2"	"s5-6" "s7-4"	"s8-9"	"s8-8"	"s11-3"	"s11-7"	"s13-3"	"s14-2"	"s16-5"	"s16-6" "s17-2"	"s20-1"	"s21-2"	"s23-4"	"s25-1"	"s30-4"	"s31-6"	"s32-2"	"s34-4" "s36-1"	"s36-2"	"s36-4" ];

data_dir = "/datahub/rawdata/tandeng/mRNA_imaging/mRNA_confocal_hamamatsu-60X-TIRF/20200314";
output_dir = "$(data_dir)_visualization"

try
    mkdir("$output_dir")
catch
end

function max_projection(cell)
    println("Processing $cell")
    #raw_img = load(File(format"TIFF", "$data_dir/HE7-11-1-80uw-const_3_$pos.ome.tiff"));
    raw_img = loadims2("$(data_dir)_result/output/d14$cell.ims");
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
