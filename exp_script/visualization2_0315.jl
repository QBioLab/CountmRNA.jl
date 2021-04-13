using FileIO
include("julia2ims.jl")

field = ["s1-1"  "s4-3"  "s4-6"  "s5-6"  "s9-5"  "s9-4"  "s10-2" "s11-4" "s12-1" "s13-3" "s13-1" "s14-1" "s14-6" "s14-4" "s14-8" "s16-5" "s16-1" "s21-1" "s22-8" "s23-1" "s27-5" "s34-3"  ]
data_dir = "/datahub/rawdata/tandeng/mRNA_imaging/mRNA_confocal_hamamatsu-60X-TIRF/20200315";
output_dir = "$(data_dir)_visualization"

try
    mkdir("$output_dir")
catch
end

function max_projection(cell)
    println("Processing $cell")
    #raw_img = load(File(format"TIFF", "$data_dir/HE7-11-1-80uw-const_3_$pos.ome.tiff"));
    raw_img = loadims2("$(data_dir)_result/output/d15$cell.ims");
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
