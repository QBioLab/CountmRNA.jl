include("splitcell.jl")
include("lineage.jl");
include("segmentation3d.jl")
include("normalization3d.jl")
include("julia2ims.jl")
using Printf
using Images
using Dates

data_dir = "/datahub/rawdata/tandeng/mRNA_imaging/mRNA_confocal_IV/20210305";
ret_dir = "/datahub/rawdata/tandeng/mRNA_imaging/mRNA_confocal_IV/20210305_result";
day = "05"


function playing(s)
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println("$(Dates.now()) Loading $s")
    local id = @sprintf("%02d", s)
    @time img = Gray.(reinterpret(N0f16, loadims("$data_dir/pos$id.ims")));
    #p_img = img .- df
    GC.gc()
	println("Load done")

    t_len = size(img)[4]
	@time markers = split_cell_LoG(img, LoG=40, thres=-1e-7);
    GC.gc()
	@time time_line, longlived_labels, livingtime, time_line_whole = 
		find_time_line(markers, shortest=80);
	if length(longlived_labels) == 0
		println("No longlived cells found, skipping")
		return nothing
	end

	@time split_contacted_cell!(time_line, longlived_labels, livingtime, time_line_whole);
	tracks = walking(time_line, longlived_labels, livingtime);

	@time longlived_maps, watershed_maps = grant_domain(img, time_line, 
		longlived_labels, livingtime, time_line_whole);

	for index in 1:length(longlived_labels)
		println("tracking cell $index")

		@time local cell = pick_cell(img , longlived_maps, longlived_labels[index], 
			tracks[:,:,index], livingtime[index]);
        #TiffImages.save("$ret_dir/d$(day)s$id-$index_id.tiff", cell)
		@time cell_nu, nucleus_size, nucleus_th = extract3dnucleus(cell);

        local index_id = @sprintf("%02d", index)
		@time cell_nu_nor, nor_para = normalize(cell_nu);
        save("$ret_dir/d$(day)s$(id)-$index_id.jld", "threshold", nucleus_th, "normal", nor_para,
             "livingtime", livingtime[index], "track", tracks[:,:,index], "size", nucleus_size)
        save2ims(reinterpret.(cell_nu_nor), "$ret_dir/d$(day)s$id-$index_id.ims", compression=2)
	end
	GC.gc()
end

println("Processing $data_dir, output $ret_dir")
for i in parse(Int, ARGS[1]):parse(Int, ARGS[2])
#for i in parse(Int, ARGS[1]):4:parse(Int, ARGS[2])
    playing(i)
end
