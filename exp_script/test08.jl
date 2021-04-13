include("splitcell.jl")
include("lineage.jl");
include("segmentation3d.jl")
include("normalization3d.jl")
include("julia2ims.jl")
using HDF5
using Images
using Dates

data_dir = "/datahub/rawdata/tandeng/mRNA_imaging/mRNA_confocal/20200308-11-2-piezo";
ret_dir = "/datahub/rawdata/tandeng/mRNA_imaging/mRNA_confocal/20200308_result";

function playing(s)
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println("$(Dates.now()) Loading $s")
	#@time img = load(File(format"TIFF","$data_dir/HE7-11-1-80uw-const_3_s$s.ome.tiff"));
    @time img = loadims("$data_dir/HED-1NLS-11-2-11uW-const-200ms_s$(s)_.ims")[12:500, 12:500,:,:];
	GC.gc()
	println("Load done")

    t_len = size(img)[4]
	@time markers = split_cell_LoG(img, LoG=40, thres=-1e-6);
    GC.gc()
	@time time_line, longlived_labels, livingtime, time_line_whole = 
		find_time_line(markers, shortest=60);
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
        save("$ret_dir/d08$(s)-$index.jld", "threshold", nucleus_th, "normal", nor_para,
             "livingtime", livingtime[index], "track", tracks[:,:,index], "size", nucleus_size)
        save2ims(reinterpret.(cell_nu_nor), "$ret_dir/d08s$(s)-$index.ims", dx=0.16, compression=0)
	end
	GC.gc()
end

println("Processing $data_dir, output $ret_dir")
for i in parse(Int, ARGS[1]):parse(Int, ARGS[2])
    playing(i)
end
