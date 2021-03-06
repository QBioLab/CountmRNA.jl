include("splitcell.jl")
include("lineage.jl");
include("segmentation3d.jl")
include("normalization3d.jl")
include("julia2ims.jl")
using Dates

data_dir = "/datahub/rawdata/tandeng/mRNA_imaging/mRNA_confocal_hamamatsu-60X-TIRF/20200320";
ret_dir = "/datahub/rawdata/tandeng/mRNA_imaging/mRNA_confocal_hamamatsu-60X-TIRF/20200320_result/output";

#for s in 3:19
#for s in 2:20
function playing(s)
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println("$(Dates.now()) Loading $s")
	@time img = load(File(format"TIFF","$data_dir/HE7-1NLS-22uw-const_1_s$s.ome.tiff"));
	GC.gc()
	println("Load done")

    t_len = size(img)[3]÷20
	@time markers = split_cell_LoG(img, t_len);
    GC.gc()
	@time time_line, longlived_labels, livingtime, time_line_whole = 
		find_time_line(markers);
	if length(longlived_labels) == 0
		println("No longlived cells found, skipping")
		return nothing
	end
	@time split_contacted_cell!(time_line, longlived_labels, livingtime, time_line_whole);
	tracks = walking(time_line, longlived_labels, livingtime);
	@time longlived_maps, watershed_maps = grant_domain(img, time_line, longlived_labels, livingtime, time_line_whole);
    #save("$ret_dir/d21s$(s)-watershed.jld", "maps", watershed_maps)
	for index in 1:length(longlived_labels)
		println("tracking cell $index")
		@time local cell = pick_cell(img, longlived_maps, longlived_labels[index], tracks[:,:,index],
								livingtime[index]);
		@time cell_nu, cell_th = extract3dnucleus(cell);
		@time cell_nu_nor, nor_para = normalize(cell_nu);
        save("$ret_dir/d20s$(s)-$index.jld", "threshold", cell_th, "normal", nor_para,
             "livingtime", livingtime[index], "track", tracks[:,:,index])
        save2ims(reinterpret.(reshape(cell_nu_nor, (512,512,20,t_len))), 
                 "$ret_dir/d20s$(s)-$index.ims")
	end
	GC.gc()
end

println("Processing $data_dir, output $ret_dir")
for i in parse(Int, ARGS[1]):parse(Int, ARGS[2])
    playing(i)
end
