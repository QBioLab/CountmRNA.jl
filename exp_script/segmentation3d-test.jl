using FileIO

include("tiffxml.jl")
include("Segmentation3D.jl")

data_dir = "/datahub/rawdata/tandeng/mRNA_imaging/mRNA_confocal_hamamatsu-60X-TIRF/20200316_result"
#cell_name = filter(x->occursin(r"^s",x), readdir(data_dir) )
cell_name = filter(x->occursin(r"^s[0-9][0-9]?-c[0-9][0-9]?.tiff$",x), readdir(data_dir) )


for i in cell_name
	print("processing $i")
	@time cell_raw = load(File(format"TIFF", "$data_dir/$i"))
	@time nucleus_all, threshold_all = extract3dnucleus(cell_raw)
    save(File(format"TIFF", "$data_dir/$i-clear.ome.tiff"), N0f16.(nucleus_all))
	save("$data_dir/$i.jld", "threshold", threshold_all)
	x, y, z_all = size(nucleus_all)
	@time embedxml(x, y, 20, z_all÷20, "$data_dir/$i-clear.ome.tiff")
end
