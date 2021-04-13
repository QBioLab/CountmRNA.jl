include("splitcell.jl")
include("lineage.jl");
include("segmentation3d.jl")
include("normalization3d.jl")
include("julia2ims.jl")
using Printf
using Images
using Dates

data_dir = "/datahub/rawdata/tandeng/mRNA_imaging/mRNA_confocal_III/20210119";
ret_dir = "/datahub/rawdata/tandeng/mRNA_imaging/mRNA_confocal_III/2021019_result";

function playing(s)
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println("$(Dates.now()) Loading $s")
    local id = @sprintf("%02d", s)
    @time img = loadims("$data_dir/pos$id.ims");
    GC.gc()
	println("Load done")

    t_len = size(img)[4]
	@time markers = split_cell_LoG(img, LoG=40, thres=-2e-7);
    GC.gc()
	@time time_line, longlived_labels, livingtime, time_line_whole = 
		find_time_line(markers, shortest=70);
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
		@time cell_nu, nucleus_size, nucleus_th = extract3dnucleus(cell);
		@time cell_nu_nor, nor_para = normalize(cell_nu);
        save("$ret_dir/d19s$(s)-$index.jld", "threshold", nucleus_th, "normal", nor_para,
             "livingtime", livingtime[index], "track", tracks[:,:,index], "size", nucleus_size)
        save2ims(reinterpret.(cell_nu_nor), "$ret_dir/d19s$(s)-$index.ims", compression=2)
	end
	GC.gc()
end

println("Processing $data_dir, output $ret_dir")
for i in parse(Int, ARGS[1]):parse(Int, ARGS[2])
#for i in parse(Int, ARGS[1]):4:parse(Int, ARGS[2])
    playing(i)
end
