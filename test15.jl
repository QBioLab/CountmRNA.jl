include("splitcell.jl")
include("lineage.jl");
include("segmentation3d.jl")
include("normalization3d.jl")
include("tiffxml.jl");
include("julia2ims.jl")
#
const data_dir = "/datahub/rawdata/tandeng/mRNA_imaging/mRNA_confocal_hamamatsu-60X-TIRF/20200315";
const ret_dir = "/datahub/rawdata/tandeng/mRNA_imaging/mRNA_confocal_hamamatsu-60X-TIRF/20200315_result/output";

#for s in 3:19
#function playing(s::Int)
#
#cell = zeros(Gray{Normed{UInt16,16}}, 512, 512, *20)

#for in 19:20
#for s in 19
function playing(s)
	println("Loading $s")
	@time img = load(File(format"TIFF","$data_dir/HE7-11-1-80uw-const_1_s$s.ome.btf.tiff"));
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
	#cell = zeros(Gray{Normed{UInt16,16}}, 512, 512, t_len*20)
	for index in 1:length(longlived_labels)
		println("tracking cell $index")
		@time cell = pick_cell(img, longlived_maps, longlived_labels[index], tracks[:,:,index],
								livingtime[index]);
		@time cell_nu, cell_th = extract3dnucleus(cell);
		@time cell_nu_nor, nor_para = normalize(cell_nu);
        save("$ret_dir/d15s$(s)-$index.jld", "threshold", cell_th, "normal", nor_para, 
             "livingtime", livingtime[index])
		save2ims( reinterpret.(reshape(cell_nu_nor,(512,512,20,t_len))), 
				 "$ret_dir/d15s$(s)-$index.ims")
		#GC.gc()
	end
	nothing
end

for img in 26:40
	@time playing(img)
end
