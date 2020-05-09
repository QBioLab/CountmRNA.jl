include("splitcell.jl")
include("lineage.jl");
include("segmentation3d.jl")
include("normalization3d.jl")
#include("tiffxml.jl");
include("julia2ims.jl")

data_dir = "/datahub/rawdata/tandeng/mRNA_imaging/mRNA_confocal_hamamatsu-60X-TIRF";
ret_dir = "/datahub/rawdata/tandeng/mRNA_imaging/mRNA_confocal_hamamatsu-60X-TIRF/20200316_result/output";

#for s in 3:19
#for s in 2:20
function playing(s)
	println("Loading $s")
	@time img = load(File(format"TIFF","$data_dir/20200316/HE7-11-1-80uw-PWM_1_s$s.ome.tiff"));
	GC.gc()
	println("Load done")
    t_len = size(img)[3]÷20
	@time markers = split_cell_LoG(img, t_len);
    GC.gc()
	@time time_line, longlived_labels, livingtime, time_line_whole = 
		find_time_line(markers);
	@time split_contacted_cell!(time_line, longlived_labels, livingtime, time_line_whole);
	tracks = walking(time_line, longlived_labels, livingtime);
	@time longlived_maps, watershed_maps = grant_domain(img, time_line, longlived_labels, livingtime, time_line_whole);
	for index in 1:length(longlived_labels)
		println("tracking cell $index")
		@time local cell = pick_cell(img, longlived_maps, longlived_labels[index], tracks[:,:,index],
								livingtime[index]);
		@time cell_nu, cell_th = extract3dnucleus(cell);
		@time cell_nu_nor, nor_para = normalize(cell_nu);
        save("$ret_dir/d16s$(s)-$index.jld", "threshold", cell_th, "normal", nor_para,
             "livingtime", livingtime[index])
        save2ims(reinterpret.(reshape(cell_nu_nor, (512,512,20,t_len))), 
                 "$ret_dir/d16s$(s)-$index.ims")
	end
	GC.gc()
end

for i in 22:25
    playing(i)
end
