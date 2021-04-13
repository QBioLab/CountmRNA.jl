include("trackingcell.jl");
include("lineage.jl");

data_dir = "/datahub/rawdata/tandeng/mRNA_imaging/mRNA_confocal_hamamatsu-60X-TIRF";
ret_dir = "/datahub/rawdata/tandeng/mRNA_imaging/mRNA_confocal_hamamatsu-60X-TIRF/20200316_result";
img = load(File(format"TIFF","$data_dir/20200316/HE7-11-1-80uw-PWM_1_s3.ome.tiff"));
img_edge, mask_edge, mask_clear = split_cell_LoG(img, size(img)[3]÷20);
save("$ret_dir/img_edge.tiff", img_edge);

nucleus_pos = [ component_centroids(label_components(mask_clear[:,:,i])) for i in 1:size(mask_clear)[3] ];

paths = generate_path(nucleus_pos);

for i in 1:length(paths)
	save("$ret_dir/s3_c$i.tiff", N0f16.(export_cell(paths[i], img, mask_edge)));
end

