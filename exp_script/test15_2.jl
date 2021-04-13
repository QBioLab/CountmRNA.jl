include("trackingcell.jl");
include("lineage.jl");

data_dir = "/datahub/rawdata/tandeng/mRNA_imaging/mRNA_confocal_hamamatsu-60X-TIRF";
ret_dir = "/datahub/rawdata/tandeng/mRNA_imaging/mRNA_confocal_hamamatsu-60X-TIRF/20200315_result";
img = load(File(format"TIFF","$data_dir/20200315/HE7-11-1-80uw-const_1_s2.ome.btf.tiff"));
img_edge, mask_edge, mask_clear = split_cell_LoG(img, size(img)[3]÷20);
save("$ret_dir/img_edge.tiff", img_edge);

nucleus_pos = [ component_centroids(label_components(mask_clear[:,:,i])) for i in 1:size(mask_clear)[3] ];

paths = generate_path(nucleus_pos);

for i in 1:length(paths)
	save("$ret_dir/s2_c$cell.tiff", N0f16.(export_cell(paths[i], img, mask_edge)));
end

